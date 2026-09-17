"""Cell-free DNA resuspension and Gator plate setup on a Hamilton STAR.

Reads a Twist-style platemap (CSV/Excel) and a Gator layout Excel sheet, resuspends
dried DNA to a normalized concentration, dispenses cell-free mastermix, and transfers
DNA into the named Gator wells.
"""

from __future__ import annotations

import asyncio
import csv
import io
import re
import zipfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import BinaryIO, Callable, Iterable, Optional, Union
import sys 

sys.path.insert(0, r"C:\Users\16122\Documents\Github\pylabrobot")

ROWS = "ABCDEFGH"
ROW_TO_CHANNEL = {row: idx for idx, row in enumerate(ROWS)}
WELL_RE = re.compile(r"^([A-Ha-h])\s*0*([1-9]|1[0-2])$")
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}

LogFn = Callable[[str], None]


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

TIP_MAX_UL = {10: 10.0, 50: 50.0, 300: 300.0}


@dataclass
class Settings:
    dna_vol_ul: float = 4.0
    mastermix_vol_ul: float = 16.0
    target_concentration: float = 40.0
    concentration_unit: str = "ng/uL"  # "ng/uL" or "nM"
    min_resuspend_ul: float = 5.0
    max_resuspend_ul: float = 120.0
    resuspend_mix_cycles: int = 5
    dest_mix_cycles: int = 3
    water_flow_rate: float = 100.0
    resuspend_dispense_flow_rate: float = 30.0
    resuspend_mix_flow_rate: float = 30.0
    resuspend_dispense_height_mm: float = 4.0
    resuspend_settle_s: float = 2.0
    mastermix_flow_rate: float = 50.0
    mastermix_lld: bool = True
    mastermix_immersion_mm: float = 2.0
    mastermix_lld_sensitivity: int = 2  # 1 = high, 4 = low
    # 2 mL tubes sit higher in the 1.5 mL 32-position insert than the carrier model.
    mastermix_tube_z_offset_mm: float = 18.0
    mastermix_min_height_mm: float = 5.0
    dna_flow_rate: float = 40.0
    dna_lld: bool = True
    dna_immersion_mm: float = 1.0
    dna_lld_sensitivity: int = 1  # 1 = high; small volumes
    dna_min_height_mm: float = 0.3
    dna_aspirate_xy_offset_mm: float = 1.2
    do_resuspend: bool = True
    do_mastermix: bool = True
    do_dna_transfer: bool = True
    simulation: bool = False
    tip_carrier_rail: int = 7
    plate_carrier_rail: int = 1
    trough_carrier_rail: int = 43
    tube_carrier_rail: int = 35
    dna_plate_site: int = 0
    gator_plate_1_site: int = 1
    gator_plate_2_site: int = 2
    water_trough_site: int = 0
    mastermix_tube_site: int = 0
    load_10ul_tips: bool = True
    load_50ul_tips: bool = False
    load_300ul_tips: bool = True
    tips_10ul_site: int = 0
    tips_10ul_extra_site: int = 1
    tips_50ul_site: int = 2
    tips_50ul_extra_site: int = -1
    tips_300ul_site: int = 3
    tips_300ul_extra_site: int = 4
    resuspend_tip_ul: int = 300
    mastermix_tip_ul: int = 300
    dna_tip_ul: int = 10


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

@dataclass
class DNASample:
    name: str
    well: str
    yield_ng: float
    length_bp: Optional[int] = None


@dataclass
class DestWell:
    plate_title: str
    plate_resource: str
    well: str
    name: str


@dataclass
class Transfer:
    name: str
    source_well: str
    dest_plate: str
    dest_plate_title: str
    dest_well: str
    yield_ng: float
    length_bp: Optional[int]
    resuspend_ul: float
    actual_conc: float
    dna_vol_ul: float
    mastermix_vol_ul: float
    warnings: list[str] = field(default_factory=list)

    @property
    def source_row(self) -> str:
        return self.source_well[0]

    @property
    def source_col(self) -> int:
        return int(self.source_well[1:])

    @property
    def dest_row(self) -> str:
        return self.dest_well[0]

    @property
    def dest_col(self) -> int:
        return int(self.dest_well[1:])


@dataclass
class RunPlan:
    transfers: list[Transfer]
    gator_titles: list[str]
    unmatched_dest_names: list[str]
    unused_dna_names: list[str]
    water_ul: float
    mastermix_ul: float
    settings: Settings


# ---------------------------------------------------------------------------
# File helpers
# ---------------------------------------------------------------------------

FileInput = Union[str, Path, bytes, BinaryIO]


def normalize_well(well: str) -> str:
    text = str(well).strip().upper().replace(" ", "")
    match = WELL_RE.match(text)
    if not match:
        raise ValueError(f"Invalid well location: {well!r}")
    return f"{match.group(1)}{int(match.group(2))}"


def resolve_file(upload_value, fallback_path: str = "") -> tuple[bytes, str]:
    """Return (bytes, filename) from an ipywidgets FileUpload value or a filesystem path."""
    if upload_value:
        if isinstance(upload_value, dict):
            filename, info = next(iter(upload_value.items()))
            content = info["content"] if isinstance(info, dict) else info.content
            return bytes(content), filename
        first = upload_value[0]
        if isinstance(first, dict):
            return bytes(first["content"]), first.get("name", "uploaded")
        return bytes(first.content), getattr(first, "name", "uploaded")
    if fallback_path:
        path = Path(fallback_path)
        if not path.exists():
            raise FileNotFoundError(f"File not found: {fallback_path}")
        return path.read_bytes(), path.name
    raise ValueError("Select a file or enter a path.")


def loaded_tip_racks(settings: Settings) -> list[tuple[str, int, int]]:
    """Return (resource_name, tip_size_ul, carrier_site) for racks on the deck."""
    racks: list[tuple[str, int, int]] = []

    def add(loaded: bool, size: int, primary: int, extra: int) -> None:
        if not loaded:
            return
        racks.append((f"tips_{size}uL_0", size, primary))
        if extra >= 0:
            racks.append((f"tips_{size}uL_1", size, extra))

    add(settings.load_10ul_tips, 10, settings.tips_10ul_site, settings.tips_10ul_extra_site)
    add(settings.load_50ul_tips, 50, settings.tips_50ul_site, settings.tips_50ul_extra_site)
    add(settings.load_300ul_tips, 300, settings.tips_300ul_site, settings.tips_300ul_extra_site)
    return racks


def racks_for_tip_size(settings: Settings, size: int) -> list[str]:
    return [name for name, tip_size, _ in loaded_tip_racks(settings) if tip_size == size]


def validate_settings(settings: Settings) -> None:
    plate_sites = [settings.dna_plate_site, settings.gator_plate_1_site, settings.gator_plate_2_site]
    if len(set(plate_sites)) != 3:
        raise ValueError("DNA plate and both Gator plates must occupy different carrier sites.")
    racks = loaded_tip_racks(settings)
    if not racks:
        raise ValueError("Load at least one tip rack.")
    sites = [site for _, _, site in racks]
    if len(set(sites)) != len(sites):
        raise ValueError("Loaded tip racks must occupy different carrier sites.")
    if not 0 <= settings.mastermix_tube_site <= 31:
        raise ValueError("Mastermix 2 mL tube site must be between 0 and 31.")
    if settings.dna_vol_ul <= 0 or settings.mastermix_vol_ul <= 0:
        raise ValueError("DNA and mastermix volumes must be greater than 0.")
    if settings.concentration_unit not in ("ng/uL", "nM"):
        raise ValueError("Concentration unit must be ng/uL or nM.")
    for step, size in (
        ("resuspend", settings.resuspend_tip_ul),
        ("mastermix", settings.mastermix_tip_ul),
        ("DNA transfer", settings.dna_tip_ul),
    ):
        if size not in TIP_MAX_UL:
            raise ValueError(f"{step} tip size must be 10, 50, or 300 µL.")
        if not racks_for_tip_size(settings, size):
            raise ValueError(f"{step} uses {size} µL tips, but that rack is not loaded.")
    if settings.dna_vol_ul > TIP_MAX_UL[settings.dna_tip_ul]:
        raise ValueError("DNA volume exceeds the selected DNA tip size.")
    if settings.mastermix_vol_ul > TIP_MAX_UL[settings.mastermix_tip_ul]:
        raise ValueError("Mastermix volume exceeds the selected mastermix tip size.")
    if not 1 <= settings.mastermix_lld_sensitivity <= 4:
        raise ValueError("Mastermix cLLD sensitivity must be 1 (high) through 4 (low).")
    if settings.mastermix_immersion_mm <= 0:
        raise ValueError("Mastermix immersion depth must be greater than 0 mm.")
    if not 0 <= settings.mastermix_tube_z_offset_mm <= 40:
        raise ValueError("Mastermix tube Z offset must be between 0 and 40 mm.")
    if settings.mastermix_min_height_mm < 2:
        raise ValueError("Mastermix minimum height above tube bottom must be at least 2 mm.")
    if settings.resuspend_dispense_flow_rate <= 0 or settings.resuspend_mix_flow_rate <= 0:
        raise ValueError("Resuspend dispense and mix flow rates must be greater than 0.")
    if settings.resuspend_dispense_height_mm < 0:
        raise ValueError("Resuspend dispense height cannot be negative.")
    if settings.resuspend_settle_s < 0:
        raise ValueError("Resuspend settle time cannot be negative.")
    if not 1 <= settings.dna_lld_sensitivity <= 4:
        raise ValueError("DNA cLLD sensitivity must be 1 (high) through 4 (low).")
    if settings.dna_immersion_mm <= 0:
        raise ValueError("DNA immersion depth must be greater than 0 mm.")
    if settings.dna_min_height_mm < 0:
        raise ValueError("DNA minimum height above well bottom cannot be negative.")


def _as_bytes(source: FileInput) -> bytes:
    if isinstance(source, (str, Path)):
        return Path(source).read_bytes()
    if hasattr(source, "read"):
        data = source.read()
        return data if isinstance(data, bytes) else data.encode("utf-8")
    if isinstance(source, bytes):
        return source
    raise TypeError(f"Unsupported file input: {type(source)}")


def _col_to_index(col: str) -> int:
    n = 0
    for char in col:
        n = n * 26 + (ord(char) - 64)
    return n


def _read_xlsx_cells(data: bytes) -> dict[tuple[int, int], str]:
    """Return {(row, col_index_1based): string_value} for the first sheet."""
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        wb = ET.fromstring(zf.read("xl/workbook.xml"))
        sheets = wb.findall("m:sheets/m:sheet", NS)
        rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        rid = sheets[0].attrib[
            "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
        ]
        target = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels}[rid]
        if not target.startswith("xl/"):
            target = "xl/" + target.lstrip("/")

        shared: list[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in root.findall("m:si", NS):
                shared.append(
                    "".join(t.text or "" for t in si.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t"))
                )

        cells: dict[tuple[int, int], str] = {}
        sheet = ET.fromstring(zf.read(target))
        for cell in sheet.findall("m:sheetData/m:row/m:c", NS):
            ref = cell.attrib.get("r", "")
            match = re.match(r"^([A-Z]+)(\d+)$", ref)
            if not match:
                continue
            col_idx = _col_to_index(match.group(1))
            row_idx = int(match.group(2))
            kind = cell.attrib.get("t")
            if kind == "s":
                value_el = cell.find("m:v", NS)
                value = shared[int(value_el.text)] if value_el is not None and value_el.text else ""
            elif kind == "inlineStr":
                value = "".join(
                    t.text or ""
                    for t in cell.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t")
                )
            else:
                value_el = cell.find("m:v", NS)
                value = value_el.text if value_el is not None and value_el.text else ""
            if value != "":
                cells[(row_idx, col_idx)] = value
        return cells


def _read_tabular(data: bytes, filename: str = "") -> list[dict[str, str]]:
    name = filename.lower()
    if name.endswith(".xlsx") or data[:2] == b"PK":
        cells = _read_xlsx_cells(data)
        if not cells:
            return []
        min_row = min(r for r, _ in cells)
        max_row = max(r for r, _ in cells)
        max_col = max(c for _, c in cells)
        headers = [
            str(cells.get((min_row, col), "")).strip()
            for col in range(1, max_col + 1)
        ]
        rows = []
        for row in range(min_row + 1, max_row + 1):
            record = {
                headers[col - 1]: str(cells.get((row, col), "")).strip()
                for col in range(1, max_col + 1)
                if headers[col - 1]
            }
            if any(record.values()):
                rows.append(record)
        return rows

    text = data.decode("utf-8-sig")
    return list(csv.DictReader(io.StringIO(text)))


def _pick_column(row: dict[str, str], candidates: Iterable[str]) -> Optional[str]:
    lookup = {re.sub(r"\s+", " ", key.strip().lower()): key for key in row}
    for candidate in candidates:
        key = lookup.get(candidate.lower())
        if key and str(row.get(key, "")).strip():
            return str(row[key]).strip()
    return None


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_platemap(source: FileInput, filename: str = "") -> list[DNASample]:
    data = _as_bytes(source)
    if not filename and isinstance(source, (str, Path)):
        filename = str(source)
    rows = _read_tabular(data, filename)
    if not rows:
        raise ValueError("Platemap is empty.")

    samples: list[DNASample] = []
    for row in rows:
        name = _pick_column(row, ["Name", "Custom Label", "Sample", "Sample Name"])
        well = _pick_column(row, ["Well Location", "Well", "Well Position", "Position"])
        yield_text = _pick_column(row, ["Yield (ng)", "Yield", "Yield ng", "ng"])
        length_text = _pick_column(
            row,
            [
                "Construct Length (Insert + Adapters)",
                "Insert Length",
                "Length",
                "bp",
            ],
        )
        if not name or not well:
            continue
        if yield_text is None:
            raise ValueError(f"Missing yield for {name} in {well}.")
        length_bp = int(float(length_text)) if length_text else None
        samples.append(
            DNASample(
                name=name,
                well=normalize_well(well),
                yield_ng=float(yield_text),
                length_bp=length_bp,
            )
        )
    if not samples:
        raise ValueError("No DNA samples found in platemap.")
    return samples


def parse_gatorsetup(source: FileInput) -> list[DestWell]:
    data = _as_bytes(source)
    if data[:2] != b"PK":
        return _parse_gator_table(data)

    cells = _read_xlsx_cells(data)
    grids = _find_plate_grids(cells)
    if grids:
        dests: list[DestWell] = []
        for idx, (title, well_map) in enumerate(grids):
            resource = f"gator_plate_{idx + 1}"
            for well, name in well_map.items():
                dests.append(
                    DestWell(
                        plate_title=title,
                        plate_resource=resource,
                        well=well,
                        name=name,
                    )
                )
        return dests
    return _parse_gator_table(data)


def _parse_gator_table(data: bytes) -> list[DestWell]:
    rows = _read_tabular(data)
    dests: list[DestWell] = []
    plate_resources: dict[str, str] = {}
    for row in rows:
        name = _pick_column(row, ["Name", "Sample", "Sample Name"])
        well = _pick_column(row, ["Well", "Well Location", "Well Position"])
        plate = _pick_column(row, ["Plate", "Plate Name", "Gator Plate"])
        if not name or not well:
            continue
        plate = plate or "gator_plate_1"
        if plate not in plate_resources:
            plate_resources[plate] = f"gator_plate_{len(plate_resources) + 1}"
        dests.append(
            DestWell(
                plate_title=plate,
                plate_resource=plate_resources[plate],
                well=normalize_well(well),
                name=name,
            )
        )
    if not dests:
        raise ValueError("No Gator sample wells found.")
    return dests


def _find_plate_grids(cells: dict[tuple[int, int], str]) -> list[tuple[str, dict[str, str]]]:
    """Find 8x12 blocks labeled A-H / 1-12 and return (title, {well: name})."""
    grids: list[tuple[str, dict[str, str]]] = []
    used_headers: set[tuple[int, int]] = set()

    for (row, col), value in sorted(cells.items()):
        if (row, col) in used_headers:
            continue
        if _cell_number(value) != 1:
            continue
        if any(_cell_number(cells.get((row, col + offset), "")) != offset + 1 for offset in range(12)):
            continue
        label_col = col - 1
        label_row = row + 1
        if not all(cells.get((label_row + i, label_col), "").strip().upper() == ROWS[i] for i in range(8)):
            continue

        used_headers.add((row, col))
        title = (
            cells.get((row - 1, label_col), "")
            or cells.get((row - 1, col), "")
            or f"Plate {len(grids) + 1}"
        ).strip()
        well_map: dict[str, str] = {}
        for r_i, plate_row in enumerate(ROWS):
            for c_i in range(12):
                name = str(cells.get((label_row + r_i, col + c_i), "")).strip()
                if name and name.upper() not in ROWS and _cell_number(name) is None:
                    well_map[f"{plate_row}{c_i + 1}"] = name
        grids.append((title, well_map))
    return grids


def _cell_number(value: object) -> Optional[int]:
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    if number.is_integer():
        return int(number)
    return None


# ---------------------------------------------------------------------------
# Volumes and worklist
# ---------------------------------------------------------------------------

def resuspension_volume_ul(
    sample: DNASample,
    settings: Settings,
) -> tuple[float, float, list[str]]:
    """Return (volume_ul, actual_concentration, warnings)."""
    warnings: list[str] = []
    unit = settings.concentration_unit
    if unit == "nM":
        if not sample.length_bp:
            raise ValueError(f"{sample.name} is missing construct length needed for nM resuspension.")
        # ng = nM * µL * bp * 650 / 1e6
        raw = sample.yield_ng / (settings.target_concentration * sample.length_bp * 650 / 1e6)
        conc_from_vol = lambda vol: sample.yield_ng / (vol * sample.length_bp * 650 / 1e6)
    else:
        raw = sample.yield_ng / settings.target_concentration
        conc_from_vol = lambda vol: sample.yield_ng / vol

    volume = raw
    if volume < settings.min_resuspend_ul:
        volume = settings.min_resuspend_ul
        warnings.append(
            f"{sample.name}: calculated {raw:.1f} µL < min {settings.min_resuspend_ul} µL"
        )
    if volume > settings.max_resuspend_ul:
        volume = settings.max_resuspend_ul
        warnings.append(
            f"{sample.name}: calculated {raw:.1f} µL > max {settings.max_resuspend_ul} µL"
        )
    volume = round(volume, 2)
    actual = round(conc_from_vol(volume), 3)
    return volume, actual, warnings


def build_plan(
    samples: list[DNASample],
    dests: list[DestWell],
    settings: Optional[Settings] = None,
) -> RunPlan:
    """Plan only names that appear on both the platemap and the Gator sheet."""
    settings = settings or Settings()
    by_name = {sample.name: sample for sample in samples}
    dest_names = {dest.name for dest in dests}
    transfers: list[Transfer] = []
    unmatched: list[str] = []

    for dest in dests:
        sample = by_name.get(dest.name)
        if sample is None:
            unmatched.append(f"{dest.name} ({dest.plate_title} {dest.well})")
            continue
        resuspend_ul, actual_conc, warnings = resuspension_volume_ul(sample, settings)
        if resuspend_ul + 1e-6 < settings.dna_vol_ul:
            warnings.append(f"{sample.name}: resuspend volume is smaller than the DNA transfer volume")
        transfers.append(
            Transfer(
                name=sample.name,
                source_well=sample.well,
                dest_plate=dest.plate_resource,
                dest_plate_title=dest.plate_title,
                dest_well=dest.well,
                yield_ng=sample.yield_ng,
                length_bp=sample.length_bp,
                resuspend_ul=resuspend_ul,
                actual_conc=actual_conc,
                dna_vol_ul=settings.dna_vol_ul,
                mastermix_vol_ul=settings.mastermix_vol_ul,
                warnings=warnings,
            )
        )

    transfers.sort(key=lambda t: (t.dest_plate, t.dest_col, t.dest_row))

    used = {t.name for t in transfers}
    unused = [sample.name for sample in samples if sample.name not in dest_names]
    titles = []
    for dest in dests:
        if dest.plate_title not in titles:
            titles.append(dest.plate_title)

    unique_sources = {t.source_well: t.resuspend_ul for t in transfers}
    water_ul = round(sum(unique_sources.values()) * 1.15, 1)
    mastermix_ul = round(len(transfers) * settings.mastermix_vol_ul * 1.15, 1)
    if mastermix_ul > 1800 and transfers:
        transfers[0].warnings.append(
            "Total mastermix exceeds ~1.8 mL; use a second 2 mL tube or reduce reactions."
        )
    return RunPlan(
        transfers=transfers,
        gator_titles=titles,
        unmatched_dest_names=unmatched,
        unused_dna_names=unused,
        water_ul=water_ul,
        mastermix_ul=mastermix_ul,
        settings=settings,
    )


def format_plan(plan: RunPlan) -> str:
    shared = sorted({t.name for t in plan.transfers})
    lines = [
        f"Shared names: {len(shared)}",
        f"Reactions: {len(plan.transfers)} (dilution + cell-free only for shared names)",
        f"Gator plates: {', '.join(plan.gator_titles) or '(none)'}",
        f"DNA transfer: {plan.settings.dna_vol_ul} µL"
        + (" with cLLD (skip bottom bubbles)" if plan.settings.dna_lld else ""),
        f"Resuspend: slow dispense from {plan.settings.resuspend_dispense_height_mm:.0f} mm, "
        f"gentle mix, {plan.settings.resuspend_settle_s:.0f} s settle",
        f"Mastermix: {plan.settings.mastermix_vol_ul} µL",
        f"Target: {plan.settings.target_concentration} {plan.settings.concentration_unit}",
        f"Water needed (with 15% extra): {plan.water_ul:.0f} µL in the water trough",
        f"Mastermix needed (with 15% extra): {plan.mastermix_ul:.0f} µL in a 2 mL tube",
        (
            f"Mastermix aspiration: cLLD, {plan.settings.mastermix_immersion_mm:.1f} mm below surface; "
            f"tube Z +{plan.settings.mastermix_tube_z_offset_mm:.0f} mm, "
            f"min height {plan.settings.mastermix_min_height_mm:.0f} mm above bottom"
            if plan.settings.mastermix_lld
            else (
                f"Mastermix aspiration: fixed height; "
                f"tube Z +{plan.settings.mastermix_tube_z_offset_mm:.0f} mm, "
                f"min height {plan.settings.mastermix_min_height_mm:.0f} mm above bottom"
            )
        ),
        f"Tips: resuspend {plan.settings.resuspend_tip_ul} µL, "
        f"mastermix {plan.settings.mastermix_tip_ul} µL, "
        f"DNA {plan.settings.dna_tip_ul} µL",
    ]
    if plan.unmatched_dest_names:
        lines.append(
            "Skipped Gator wells (name not in platemap): " + ", ".join(plan.unmatched_dest_names)
        )
    if plan.unused_dna_names:
        lines.append(
            "Skipped platemap DNA (name not in gatorsetup): " + ", ".join(plan.unused_dna_names)
        )
    lines.append("")
    lines.append(
        f"{'Name':<16} {'Src':<5} {'Dest plate':<16} {'Dest':<5} "
        f"{'Yield ng':>8} {'Resuspend':>10} {'Conc':>8}"
    )
    for t in plan.transfers:
        lines.append(
            f"{t.name:<16} {t.source_well:<5} {t.dest_plate_title:<16} {t.dest_well:<5} "
            f"{t.yield_ng:8.0f} {t.resuspend_ul:10.1f} {t.actual_conc:8.2f}"
        )
        for warning in t.warnings:
            lines.append(f"  warning: {warning}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 8-channel grouping
# ---------------------------------------------------------------------------

def _row_index(well: str) -> int:
    return ROW_TO_CHANNEL[well[0]]


def can_parallelize(group: list[Transfer]) -> bool:
    if not group or len(group) > 8:
        return False
    ordered = sorted(group, key=lambda t: _row_index(t.dest_well))
    dest_plate = {t.dest_plate for t in ordered}
    dest_col = {t.dest_col for t in ordered}
    src_col = {t.source_col for t in ordered}
    if len(dest_plate) != 1 or len(dest_col) != 1 or len(src_col) != 1:
        return False
    dest_idx = [_row_index(t.dest_well) for t in ordered]
    src_idx = [_row_index(t.source_well) for t in ordered]
    dest_gaps = [b - a for a, b in zip(dest_idx, dest_idx[1:])]
    src_gaps = [b - a for a, b in zip(src_idx, src_idx[1:])]
    return dest_gaps == src_gaps and all(gap >= 1 for gap in dest_gaps)


def split_parallel_transfers(transfers: list[Transfer]) -> list[list[Transfer]]:
    batches: list[list[Transfer]] = []
    remaining = list(transfers)
    remaining.sort(key=lambda t: (t.dest_plate, t.dest_col, _row_index(t.dest_well)))
    while remaining:
        seed = remaining.pop(0)
        batch = [seed]
        keep: list[Transfer] = []
        for item in remaining:
            if can_parallelize(batch + [item]):
                batch.append(item)
            else:
                keep.append(item)
        remaining = keep
        batches.append(sorted(batch, key=lambda t: _row_index(t.dest_well)))
    return batches


def group_by_source_column(transfers: list[Transfer]) -> list[list[Transfer]]:
    columns: dict[int, list[Transfer]] = {}
    seen = set()
    unique: list[Transfer] = []
    for transfer in transfers:
        if transfer.source_well in seen:
            continue
        seen.add(transfer.source_well)
        unique.append(transfer)
        columns.setdefault(transfer.source_col, []).append(transfer)
    return [sorted(group, key=lambda t: _row_index(t.source_well)) for _, group in sorted(columns.items())]


def group_dest_columns(transfers: list[Transfer]) -> list[list[Transfer]]:
    columns: dict[tuple[str, int], list[Transfer]] = {}
    for transfer in transfers:
        columns.setdefault((transfer.dest_plate, transfer.dest_col), []).append(transfer)
    return [
        sorted(group, key=lambda t: _row_index(t.dest_well))
        for _, group in sorted(columns.items())
    ]


# ---------------------------------------------------------------------------
# Deck + protocol (pylabrobot imported lazily)
# ---------------------------------------------------------------------------

class TipCursor:
    """Allocate unused tips without skipping leftover rows in a column."""

    def __init__(self, rack_names: list[str]):
        self.rack_names = rack_names
        self.used: set[tuple[str, int, str]] = set()

    def allocate(self, channels: list[int]) -> tuple[str, int, list[str]]:
        if not channels:
            raise ValueError("Need at least one channel to allocate tips.")
        offsets = [channel - channels[0] for channel in channels]
        for rack_name in self.rack_names:
            for col in range(1, 13):
                for start in range(8):
                    rows: list[str] = []
                    fits = True
                    for offset in offsets:
                        idx = start + offset
                        if idx >= 8:
                            fits = False
                            break
                        row = ROWS[idx]
                        if (rack_name, col, row) in self.used:
                            fits = False
                            break
                        rows.append(row)
                    if fits and len(rows) == len(channels):
                        for row in rows:
                            self.used.add((rack_name, col, row))
                        return rack_name, col, rows
        raise RuntimeError("Out of tips. Add another tip rack on the tip carrier.")


def setup_deck(settings: Settings):
    from pylabrobot.liquid_handling import LiquidHandler
    from pylabrobot.liquid_handling.backends import STARBackend
    from pylabrobot.liquid_handling.backends.chatterbox import LiquidHandlerChatterboxBackend
    from pylabrobot.resources import (
        PLT_CAR_L5AC_A00,
        TIP_CAR_480_A00,
        Coordinate,
        Cor_96_wellplate_360ul_Fb,
        hamilton_96_tiprack_10uL_filter,
        hamilton_96_tiprack_50uL_filter,
        hamilton_96_tiprack_300uL_filter,
        hamilton_1_trough_60ml_Vb,
    )
    from pylabrobot.resources.eppendorf.tubes import Eppendorf_DNA_LoBind_2ml_Ub
    from pylabrobot.resources.hamilton import STARDeck
    from pylabrobot.resources.hamilton.trough_carriers import Trough_CAR_5R60_A00
    from pylabrobot.resources.hamilton.tube_carriers import (
        hamilton_tube_carrier_32_a00_insert_eppendorf_1_5mL,
    )

    validate_settings(settings)

    if settings.simulation:
        backend = LiquidHandlerChatterboxBackend(num_channels=8)
    else:
        backend = STARBackend(read_timeout=2400)

    lh = LiquidHandler(
        backend=backend,
        deck=STARDeck(core_grippers="1000uL-5mL-on-waste"),
    )

    tip_builders = {
        10: hamilton_96_tiprack_10uL_filter,
        50: hamilton_96_tiprack_50uL_filter,
        300: hamilton_96_tiprack_300uL_filter,
    }
    tip_car = TIP_CAR_480_A00(name="tip_carrier")
    for name, size, site in loaded_tip_racks(settings):
        tip_car[site] = tip_builders[size](name=name)
    lh.deck.assign_child_resource(tip_car, rails=settings.tip_carrier_rail)

    plt_car = PLT_CAR_L5AC_A00(name="plate_carrier")
    plt_car[settings.dna_plate_site] = Cor_96_wellplate_360ul_Fb(name="dna_plate")
    plt_car[settings.gator_plate_1_site] = Cor_96_wellplate_360ul_Fb(name="gator_plate_1")
    plt_car[settings.gator_plate_2_site] = Cor_96_wellplate_360ul_Fb(name="gator_plate_2")
    lh.deck.assign_child_resource(plt_car, rails=settings.plate_carrier_rail)

    trough_car = Trough_CAR_5R60_A00(name="trough_carrier")
    trough_car[settings.water_trough_site] = hamilton_1_trough_60ml_Vb(name="water_trough")
    lh.deck.assign_child_resource(trough_car, rails=settings.trough_carrier_rail)

    tube_car = hamilton_tube_carrier_32_a00_insert_eppendorf_1_5mL(name="tube_carrier")
    tube_car[settings.mastermix_tube_site].assign_child_resource(
        Eppendorf_DNA_LoBind_2ml_Ub(name="mastermix_tube"),
        location=Coordinate(0, 0, settings.mastermix_tube_z_offset_mm),
    )
    lh.deck.assign_child_resource(tube_car, rails=settings.tube_carrier_rail)

    return lh


async def connect_liquid_handler(lh, settings: Settings) -> str:
    from pylabrobot.resources import set_tip_tracking, set_volume_tracking

    await lh.setup(skip_iswap=True) if not settings.simulation else await lh.setup()
    set_volume_tracking(enabled=False)
    set_tip_tracking(enabled=False)
    summary = lh.summary()
    return summary if isinstance(summary, str) else "Deck setup complete."


def _wells(plate, names: list[str]):
    return [plate.get_item(name) for name in names]


def _spaced_channels(rows: list[str]) -> list[int]:
    """Map plate rows onto 8-channel indices while preserving 9 mm gaps."""
    indexes = [_row_index(f"{row}1") for row in rows]
    base = min(indexes)
    return [idx - base for idx in indexes]


async def _pick_tips(lh, cursor: TipCursor, channels: list[int]):
    rack_name, col, tip_rows = cursor.allocate(channels)
    rack = lh.deck.get_resource(rack_name)
    spots = [rack.get_item(f"{row}{col}") for row in tip_rows]
    await lh.pick_up_tips(spots, use_channels=channels)
    return rack_name, col, tip_rows


async def _mix(
    lh,
    wells,
    vols: list[float],
    channels: list[int],
    cycles: int,
    flow_rate: float,
    tip_max_ul: float,
    mix_fraction: float = 0.4,
    leave_ul: float = 2.0,
):
    """Mix without emptying the well, which traps air bubbles at the bottom."""
    if cycles <= 0:
        return
    mix_vols = [
        max(1.0, min(v * mix_fraction, tip_max_ul * 0.8, max(1.0, v - leave_ul)))
        for v in vols
    ]
    rates = [flow_rate] * len(wells)
    settle = [0.4] * len(wells)
    for _ in range(cycles):
        await lh.aspirate(
            wells,
            vols=mix_vols,
            use_channels=channels,
            flow_rates=rates,
            settling_time=settle,
        )
        await lh.dispense(
            wells,
            vols=mix_vols,
            use_channels=channels,
            flow_rates=rates,
            settling_time=settle,
        )


def _mastermix_aspirate_kwargs(lh, tube, settings: Settings) -> dict:
    """Keep the tip above the physical tube bottom, then cLLD to the surface.

    The 32-position carrier model is for 1.5 mL inserts. A 2 mL tube sits higher, so
    the modeled cavity bottom is raised by mastermix_tube_z_offset_mm. minimum_height
    is an additional floor so LLD cannot search into the plastic if detection misses.
    """
    cavity_bottom = tube.get_location_wrt(lh.deck).z + tube.material_z_thickness
    min_z = cavity_bottom + settings.mastermix_min_height_mm
    kwargs: dict = {
        "liquid_height": [settings.mastermix_min_height_mm],
        "minimum_height": [min_z],
    }
    if not settings.mastermix_lld:
        return kwargs
    try:
        from pylabrobot.liquid_handling.backends.hamilton.STAR_backend import STARBackend
    except ImportError:
        return kwargs
    if not isinstance(lh.backend, STARBackend):
        return kwargs
    kwargs.update(
        {
            "lld_mode": [STARBackend.LLDMode.GAMMA],
            "immersion_depth": [settings.mastermix_immersion_mm],
            "gamma_lld_sensitivity": [settings.mastermix_lld_sensitivity],
        }
    )
    return kwargs


def _dna_aspirate_kwargs(lh, wells, settings: Settings) -> dict:
    """Find liquid instead of a bottom-center air bubble."""
    from pylabrobot.resources import Coordinate

    n = len(wells)
    kwargs: dict = {
        "liquid_height": [settings.dna_min_height_mm] * n,
        "offsets": [Coordinate(settings.dna_aspirate_xy_offset_mm, 0, 0) for _ in wells],
    }
    try:
        from pylabrobot.liquid_handling.backends.hamilton.STAR_backend import STARBackend
    except ImportError:
        return kwargs
    if not isinstance(lh.backend, STARBackend):
        return kwargs
    bottoms = [
        well.get_location_wrt(lh.deck).z + well.material_z_thickness for well in wells
    ]
    kwargs["minimum_height"] = [bottom + settings.dna_min_height_mm for bottom in bottoms]
    if settings.dna_lld:
        kwargs.update(
            {
                "lld_mode": [STARBackend.LLDMode.GAMMA] * n,
                "immersion_depth": [settings.dna_immersion_mm] * n,
                "gamma_lld_sensitivity": [settings.dna_lld_sensitivity] * n,
            }
        )
    return kwargs


async def _aspirate_tube_sequential(
    lh,
    tube,
    vols: list[float],
    channels: list[int],
    flow_rate: float,
    settings: Optional[Settings] = None,
):
    """Aspirate from one 2 mL tube one channel at a time, using cLLD when available."""
    settings = settings or Settings()
    lld_kwargs = _mastermix_aspirate_kwargs(lh, tube, settings)
    for vol, channel in zip(vols, channels):
        await lh.aspirate(
            [tube],
            vols=[vol],
            use_channels=[channel],
            flow_rates=[flow_rate],
            **lld_kwargs,
        )


async def run_protocol(lh, plan: RunPlan, log: LogFn = print) -> None:
    settings = plan.settings
    dna_plate = lh.deck.get_resource("dna_plate")
    water = lh.deck.get_resource("water_trough")
    mastermix = lh.deck.get_resource("mastermix_tube")
    cursors: dict[int, TipCursor] = {}
    for size in (settings.resuspend_tip_ul, settings.mastermix_tip_ul, settings.dna_tip_ul):
        if size not in cursors:
            cursors[size] = TipCursor(racks_for_tip_size(settings, size))
    resuspend_cursor = cursors[settings.resuspend_tip_ul]
    mastermix_cursor = cursors[settings.mastermix_tip_ul]
    dna_cursor = cursors[settings.dna_tip_ul]
    resuspend_max = TIP_MAX_UL[settings.resuspend_tip_ul]
    dna_max = TIP_MAX_UL[settings.dna_tip_ul]
    shared = {t.name for t in plan.transfers}
    log(
        f"Running {len(plan.transfers)} reactions for {len(shared)} names shared by both sheets."
    )
    if plan.unused_dna_names:
        log(
            f"Skipping {len(plan.unused_dna_names)} platemap names not in gatorsetup "
            "(no dilution or cell-free)."
        )
    if plan.unmatched_dest_names:
        log(
            f"Skipping {len(plan.unmatched_dest_names)} gatorsetup wells not in platemap "
            "(no dilution or cell-free)."
        )

    if settings.do_resuspend:
        max_resuspend = max(t.resuspend_ul for t in plan.transfers)
        if max_resuspend > resuspend_max:
            raise ValueError(
                f"Resuspend volume {max_resuspend:.1f} µL exceeds {settings.resuspend_tip_ul} µL tips."
            )
        log("Resuspending DNA to normalized concentration (slow dispense, gentle mix)...")
        for group in group_by_source_column(plan.transfers):
            rows = [t.source_row for t in group]
            channels = _spaced_channels(rows)
            vols = [t.resuspend_ul for t in group]
            wells = _wells(dna_plate, [t.source_well for t in group])
            log(
                f"  water -> DNA {group[0].source_well}-{group[-1].source_well} "
                f"({min(vols):.1f}-{max(vols):.1f} µL, {settings.resuspend_tip_ul} µL tips)"
            )
            await _pick_tips(lh, resuspend_cursor, channels)
            await lh.aspirate(
                [water] * len(channels),
                vols=vols,
                use_channels=channels,
                flow_rates=[settings.water_flow_rate] * len(group),
                spread="wide",
            )
            await lh.dispense(
                wells,
                vols=vols,
                use_channels=channels,
                flow_rates=[settings.resuspend_dispense_flow_rate] * len(group),
                liquid_height=[settings.resuspend_dispense_height_mm] * len(group),
            )
            await _mix(
                lh,
                wells,
                vols,
                channels,
                settings.resuspend_mix_cycles,
                settings.resuspend_mix_flow_rate,
                resuspend_max,
            )
            if settings.resuspend_settle_s > 0:
                await asyncio.sleep(settings.resuspend_settle_s)
            await lh.discard_tips()

    if settings.do_mastermix:
        if settings.mastermix_lld:
            log(
                "Dispensing cell-free mastermix from 2 mL tube "
                f"(cLLD, {settings.mastermix_immersion_mm:.1f} mm below surface, "
                f"Z floor +{settings.mastermix_tube_z_offset_mm:.0f} mm / "
                f"{settings.mastermix_min_height_mm:.0f} mm above bottom)..."
            )
        else:
            log(
                "Dispensing cell-free mastermix from 2 mL tube "
                f"(no LLD, Z floor +{settings.mastermix_tube_z_offset_mm:.0f} mm / "
                f"{settings.mastermix_min_height_mm:.0f} mm above bottom)..."
            )
        mm_tips_loaded = False
        mm_rows: list[str] = []
        for group in group_dest_columns(plan.transfers):
            rows = [t.dest_row for t in group]
            channels = _spaced_channels(rows)
            plate = lh.deck.get_resource(group[0].dest_plate)
            wells = _wells(plate, [t.dest_well for t in group])
            vols = [t.mastermix_vol_ul for t in group]
            if not mm_tips_loaded:
                await _pick_tips(lh, mastermix_cursor, channels)
                mm_tips_loaded = True
                mm_rows = rows
            elif rows != mm_rows:
                await lh.discard_tips()
                await _pick_tips(lh, mastermix_cursor, channels)
                mm_rows = rows
            log(
                f"  mastermix -> {group[0].dest_plate_title} "
                f"{group[0].dest_well}-{group[-1].dest_well} ({vols[0]} µL, sequential)"
            )
            await _aspirate_tube_sequential(
                lh, mastermix, vols, channels, settings.mastermix_flow_rate, settings
            )
            await lh.dispense(
                wells,
                vols=vols,
                use_channels=channels,
                flow_rates=[settings.mastermix_flow_rate] * len(group),
            )
        if mm_tips_loaded:
            await lh.discard_tips()

    if settings.do_dna_transfer:
        log("Transferring resuspended DNA into Gator wells...")
        for group in split_parallel_transfers(plan.transfers):
            src_rows = [t.source_row for t in group]
            channels = _spaced_channels(src_rows)
            await _pick_tips(lh, dna_cursor, channels)
            src_wells = _wells(dna_plate, [t.source_well for t in group])
            dest_plate = lh.deck.get_resource(group[0].dest_plate)
            dest_wells = _wells(dest_plate, [t.dest_well for t in group])
            vols = [t.dna_vol_ul for t in group]
            rates = [settings.dna_flow_rate] * len(group)
            log(
                f"  DNA {group[0].source_well}-{group[-1].source_well} -> "
                f"{group[0].dest_plate_title} {group[0].dest_well}-{group[-1].dest_well} "
                f"({vols[0]} µL, {settings.dna_tip_ul} µL tips"
                f"{', cLLD' if settings.dna_lld else ''})"
            )
            await lh.aspirate(
                src_wells,
                vols=vols,
                use_channels=channels,
                flow_rates=rates,
                **_dna_aspirate_kwargs(lh, src_wells, settings),
            )
            blowout = min(3.0, max(1.0, dna_max - vols[0] - 1.0))
            await lh.dispense(
                dest_wells,
                vols=vols,
                use_channels=channels,
                flow_rates=rates,
                blow_out_air_volume=[blowout] * len(group),
            )
            dest_total = [t.dna_vol_ul + t.mastermix_vol_ul for t in group]
            await _mix(
                lh,
                dest_wells,
                dest_total,
                channels,
                settings.dest_mix_cycles,
                settings.dna_flow_rate,
                dna_max,
            )
            await lh.discard_tips()

    log("Protocol complete.")


def load_run(platemap: FileInput, gatorsetup: FileInput, settings: Optional[Settings] = None, platemap_name: str = "") -> RunPlan:
    settings = settings or Settings()
    validate_settings(settings)
    samples = parse_platemap(platemap, filename=platemap_name)
    dests = parse_gatorsetup(gatorsetup)
    plan = build_plan(samples, dests, settings)
    if not plan.transfers:
        raise ValueError(
            "No shared names between the platemap and gatorsetup sheets."
        )
    return plan


def settings_from_dict(values: dict) -> Settings:
    base = Settings()
    for key, value in values.items():
        if hasattr(base, key):
            setattr(base, key, value)
    return base


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Preview a cell-free Gator worklist.")
    parser.add_argument("platemap")
    parser.add_argument("gatorsetup")
    parser.add_argument("--ng-ul", type=float, default=40.0)
    parser.add_argument("--nM", type=float, default=None)
    args = parser.parse_args()
    cfg = Settings(target_concentration=args.ng_ul)
    if args.nM is not None:
        cfg.concentration_unit = "nM"
        cfg.target_concentration = args.nM
    print(format_plan(load_run(args.platemap, args.gatorsetup, cfg, args.platemap)))
