#!/usr/bin/env python3
"""Replicate-aware inventory and downloader for the ENCODE/JOLI benchmark.

The inventory command only retrieves JSON metadata. The download command first
rebuilds the inventory and then downloads the selected data files.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable


API_ROOT = "https://www.encodeproject.org"
DEFAULT_CONFIG = Path(__file__).with_name("encode_2026_experiments.csv")
DEFAULT_OUTPUT = Path(
    "/gpfs/commons/home/atalukder/RNA_Splicing/data/ENCODE_2026"
)
USER_AGENT = "JOLI-ENCODE-downloader/1.0"

LR_OUTPUT_TYPES = {
    "reads": "reads",
    "transcript quantifications": "transcript_quantifications",
    "transcriptome annotations": "transcriptome_annotations",
}
SR_OUTPUT_TYPES = {
    "reads": "reads",
    "transcript quantifications": "transcript_quantifications",
}
MANIFEST_FIELDS = [
    "cell_line",
    "arm",
    "experiment_accession",
    "experiment_assay",
    "experiment_description",
    "experiment_status",
    "file_accession",
    "file_status",
    "biological_replicates",
    "technical_replicates",
    "biosample_accessions",
    "library_accessions",
    "file_format",
    "file_format_type",
    "output_type",
    "read_mate",
    "paired_with",
    "assembly",
    "genome_annotation",
    "platform",
    "md5sum",
    "file_size",
    "download_url",
    "selected",
    "selection_reason",
    "assignment_status",
    "local_path",
]


@dataclass(frozen=True)
class ExperimentSpec:
    cell_line: str
    arm: str
    accession: str


def request_json(path_or_url: str, attempts: int = 4) -> dict[str, Any]:
    if path_or_url.startswith("http"):
        url = path_or_url
    else:
        url = f"{API_ROOT}{path_or_url}"
    separator = "&" if "?" in url else "?"
    if "format=json" not in url:
        url += f"{separator}format=json"
    request = urllib.request.Request(
        url, headers={"Accept": "application/json", "User-Agent": USER_AGENT}
    )
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                return json.load(response)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            if attempt == attempts:
                raise
            time.sleep(2 ** (attempt - 1))
    raise RuntimeError("unreachable")


def embedded(value: Any) -> dict[str, Any]:
    """Return an embedded object, fetching an ENCODE @id when necessary."""
    if isinstance(value, dict):
        return value
    if isinstance(value, str) and value.startswith("/"):
        return request_json(value)
    return {}


def clean_accession(value: str) -> str:
    match = re.search(r"(ENCSR[A-Z0-9]+)", value.upper())
    if not match:
        raise ValueError(f"Cannot find an ENCSR accession in {value!r}")
    return match.group(1)


def load_specs(path: Path) -> list[ExperimentSpec]:
    specs: list[ExperimentSpec] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"cell_line", "arm", "experiment_accession"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain columns: {', '.join(sorted(required))}")
        for line_number, row in enumerate(reader, start=2):
            arm = row["arm"].strip().upper()
            if arm not in {"LR", "SR"}:
                raise ValueError(f"{path}:{line_number}: arm must be LR or SR")
            specs.append(
                ExperimentSpec(
                    cell_line=row["cell_line"].strip(),
                    arm=arm,
                    accession=clean_accession(row["experiment_accession"]),
                )
            )
    duplicates = [x for x in {s.accession for s in specs} if sum(s.accession == x for s in specs) > 1]
    if duplicates:
        raise ValueError(f"Duplicate experiment accessions: {', '.join(sorted(duplicates))}")
    return specs


def listify(value: Any) -> list[Any]:
    if value in (None, ""):
        return []
    return value if isinstance(value, list) else [value]


def join_values(values: Iterable[Any]) -> str:
    cleaned = sorted({str(v) for v in values if v not in (None, "")})
    return ";".join(cleaned)


def accession_from_object(value: Any) -> str:
    if isinstance(value, dict):
        return str(value.get("accession", ""))
    if isinstance(value, str):
        match = re.search(r"/(ENC[A-Z0-9]+)/", value)
        return match.group(1) if match else value
    return ""


def replicate_map(experiment: dict[str, Any]) -> dict[int, dict[str, set[str]]]:
    result: dict[int, dict[str, set[str]]] = {}
    for raw_replicate in experiment.get("replicates", []):
        replicate = embedded(raw_replicate)
        bio = replicate.get("biological_replicate_number")
        if not isinstance(bio, int):
            continue
        entry = result.setdefault(
            bio,
            {
                "biosamples": set(),
                "biosample_summaries": set(),
                "libraries": set(),
                "technical_replicates": set(),
            },
        )
        tech = replicate.get("technical_replicate_number")
        if isinstance(tech, int):
            entry["technical_replicates"].add(str(tech))
        library = embedded(replicate.get("library"))
        library_accession = accession_from_object(library or replicate.get("library"))
        if library_accession:
            entry["libraries"].add(library_accession)
        biosample_value = library.get("biosample") if library else None
        biosample = embedded(biosample_value)
        biosample_accession = accession_from_object(biosample or biosample_value)
        if biosample_accession:
            entry["biosamples"].add(biosample_accession)
        biosample_summary = biosample.get("summary") if biosample else None
        if biosample_summary:
            entry["biosample_summaries"].add(str(biosample_summary))
    return result


def object_label(value: Any) -> str:
    if isinstance(value, dict):
        obj = value
        return str(
            obj.get("title")
            or obj.get("name")
            or obj.get("summary")
            or obj.get("accession")
            or ""
        )
    return accession_from_object(value)


def display_date(value: Any) -> str:
    text = str(value or "")
    if not text:
        return ""
    candidate = text.replace("Z", "+00:00")
    try:
        parsed = datetime.fromisoformat(candidate)
        return parsed.strftime("%B %-d, %Y")
    except (ValueError, TypeError):
        return text


def experiment_libraries(experiment: dict[str, Any]) -> list[dict[str, Any]]:
    libraries: dict[str, dict[str, Any]] = {}
    for raw_replicate in experiment.get("replicates", []):
        replicate = embedded(raw_replicate)
        library = embedded(replicate.get("library"))
        if library:
            key = str(library.get("accession") or library.get("@id") or id(library))
            libraries[key] = library
    return list(libraries.values())


def aggregate_fields(objects: Iterable[dict[str, Any]], *field_names: str) -> str:
    values: list[str] = []
    for obj in objects:
        for field_name in field_names:
            for value in listify(obj.get(field_name)):
                if isinstance(value, dict) or (isinstance(value, str) and value.startswith("/")):
                    label = object_label(value)
                else:
                    label = str(value)
                if label:
                    values.append(label)
    return join_values(values)


def experiment_platforms(experiment: dict[str, Any]) -> str:
    platforms: list[str] = []
    for raw_file in experiment.get("files", []):
        file_object = embedded(raw_file)
        for raw_platform in listify(file_object.get("platform")):
            label = object_label(raw_platform)
            if label:
                platforms.append(label)
    return join_values(platforms)


def audit_counts(experiment: dict[str, Any]) -> str:
    audit = experiment.get("audit", {})
    if not isinstance(audit, dict):
        return ""
    labels = []
    for level in ("ERROR", "WARNING", "NOT_COMPLIANT", "INTERNAL_ACTION"):
        count = len(listify(audit.get(level)))
        if count:
            labels.append(f"{level.lower().replace('_', ' ')}={count}")
    return ", ".join(labels)


def file_replicates(file_object: dict[str, Any]) -> tuple[list[int], list[int]]:
    biological = [x for x in listify(file_object.get("biological_replicates")) if isinstance(x, int)]
    technical = [x for x in listify(file_object.get("technical_replicates")) if isinstance(x, int)]
    replicate_value = file_object.get("replicate")
    if replicate_value and not biological:
        replicate = embedded(replicate_value)
        bio = replicate.get("biological_replicate_number")
        tech = replicate.get("technical_replicate_number")
        if isinstance(bio, int):
            biological = [bio]
        if isinstance(tech, int):
            technical = [tech]
    return sorted(set(biological)), sorted(set(technical))


def normalize_output_type(value: Any) -> str:
    return " ".join(str(value or "").strip().lower().replace("_", " ").split())


def select_file(file_object: dict[str, Any], arm: str) -> tuple[bool, str, str]:
    status = str(file_object.get("status", "")).lower()
    if status != "released":
        return False, "", f"file status is {status or 'unknown'}, not released"
    output_type = normalize_output_type(file_object.get("output_type"))
    allowed = LR_OUTPUT_TYPES if arm == "LR" else SR_OUTPUT_TYPES
    if output_type not in allowed:
        return False, "", f"output type {output_type or 'unknown'} is not requested for {arm}"
    file_format = str(file_object.get("file_format", "")).lower()
    if output_type == "reads":
        allowed_formats = {"fastq", "fasta"} if arm == "LR" else {"fastq"}
        if file_format not in allowed_formats:
            return False, "", f"read format {file_format or 'unknown'} is not accepted for {arm}"
    elif output_type == "transcript quantifications" and file_format not in {"tsv", "txt"}:
        return False, "", f"quantification format {file_format or 'unknown'} is not TSV/text"
    elif output_type == "transcriptome annotations" and file_format not in {"gtf", "gff", "gff3"}:
        return False, "", f"annotation format {file_format or 'unknown'} is not GTF/GFF"
    return True, allowed[output_type], "requested output type and format"


def safe_component(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return value.strip("._") or "unknown"


def file_name(file_object: dict[str, Any]) -> str:
    href = str(file_object.get("href", ""))
    name = Path(href.rstrip("/")).name if href else ""
    if name:
        return safe_component(name)
    accession = str(file_object.get("accession", "file"))
    file_format = str(file_object.get("file_format", "dat"))
    return f"{safe_component(accession)}.{safe_component(file_format)}"


def file_url(file_object: dict[str, Any]) -> str:
    href = str(file_object.get("href", ""))
    if href.startswith("http"):
        return href
    return f"{API_ROOT}{href}" if href else ""


def platform_name(file_object: dict[str, Any], experiment: dict[str, Any]) -> str:
    values: list[str] = []
    for raw in listify(file_object.get("platform")):
        obj = embedded(raw)
        values.append(str(obj.get("title") or obj.get("accession") or accession_from_object(raw)))
    if not values:
        for raw in experiment.get("replicates", []):
            replicate = embedded(raw)
            library = embedded(replicate.get("library"))
            for platform in listify(library.get("platform")):
                obj = embedded(platform)
                values.append(str(obj.get("title") or obj.get("accession") or accession_from_object(platform)))
    return join_values(values)


def inventory_experiment(spec: ExperimentSpec, output_root: Path) -> list[dict[str, str]]:
    experiment = request_json(f"/experiments/{spec.accession}/?frame=embedded")
    actual_accession = experiment.get("accession")
    if actual_accession != spec.accession:
        raise RuntimeError(f"Requested {spec.accession}, received {actual_accession!r}")
    experiment_dir = output_root / spec.arm / spec.accession
    experiment_dir.mkdir(parents=True, exist_ok=True)
    with (experiment_dir / "experiment_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(experiment, handle, indent=2, sort_keys=True)
        handle.write("\n")

    rep_info = replicate_map(experiment)
    rows: list[dict[str, str]] = []
    for raw_file in experiment.get("files", []):
        file_object = embedded(raw_file)
        if not file_object:
            continue
        biological, technical = file_replicates(file_object)
        selected, category, reason = select_file(file_object, spec.arm)
        biosamples = set()
        libraries = set()
        for bio in biological:
            biosamples.update(rep_info.get(bio, {}).get("biosamples", set()))
            libraries.update(rep_info.get(bio, {}).get("libraries", set()))

        if len(biological) == 1:
            bio = biological[0]
            biosample_label = join_values(biosamples) or "NO_BIOSAMPLE"
            replicate_dir = f"replicate_{bio}_{safe_component(biosample_label)}"
            assignment = "replicate_resolved"
        elif len(biological) > 1:
            replicate_dir = "shared_or_unresolved"
            assignment = "pooled_replicates"
        else:
            replicate_dir = "shared_or_unresolved"
            assignment = "no_biological_replicate"

        destination = experiment_dir / replicate_dir
        if selected:
            destination = destination / category / file_name(file_object)
        else:
            destination = Path("")
        paired_with = accession_from_object(file_object.get("paired_with"))
        read_mate = str(file_object.get("paired_end", ""))
        if read_mate in {"1", "2"}:
            read_mate = f"R{read_mate}"
        row = {
            "cell_line": spec.cell_line,
            "arm": spec.arm,
            "experiment_accession": spec.accession,
            "experiment_assay": str(experiment.get("assay_title") or experiment.get("assay_term_name") or ""),
            "experiment_description": str(experiment.get("description", "")),
            "experiment_status": str(experiment.get("status", "")),
            "file_accession": str(file_object.get("accession", "")),
            "file_status": str(file_object.get("status", "")),
            "biological_replicates": join_values(biological),
            "technical_replicates": join_values(technical),
            "biosample_accessions": join_values(biosamples),
            "library_accessions": join_values(libraries),
            "file_format": str(file_object.get("file_format", "")),
            "file_format_type": join_values(listify(file_object.get("file_format_type"))),
            "output_type": str(file_object.get("output_type", "")),
            "read_mate": read_mate,
            "paired_with": paired_with,
            "assembly": accession_from_object(file_object.get("assembly")),
            "genome_annotation": str(file_object.get("genome_annotation", "")),
            "platform": platform_name(file_object, experiment),
            "md5sum": str(file_object.get("md5sum", "")),
            "file_size": str(file_object.get("file_size", "")),
            "download_url": file_url(file_object),
            "selected": "yes" if selected else "no",
            "selection_reason": reason,
            "assignment_status": assignment,
            "local_path": str(destination) if selected else "",
        }
        rows.append(row)

    rows.sort(key=lambda row: (row["biological_replicates"], row["output_type"], row["file_accession"]))
    write_csv(experiment_dir / "experiment_file_manifest.csv", rows)
    write_summary(experiment_dir / "experiment_summary.txt", experiment, spec, rep_info, rows)
    return rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def write_summary(
    path: Path,
    experiment: dict[str, Any],
    spec: ExperimentSpec,
    rep_info: dict[int, dict[str, set[str]]],
    rows: list[dict[str, str]],
) -> None:
    selected = sum(row["selected"] == "yes" for row in rows)
    libraries = experiment_libraries(experiment)
    biosample_ontology = embedded(experiment.get("biosample_ontology"))
    lab = object_label(experiment.get("lab"))
    award = embedded(experiment.get("award"))
    project = str(award.get("project") or experiment.get("project") or "")

    starting_quantities = []
    for library in libraries:
        amount = library.get("nucleic_acid_starting_quantity", library.get("starting_amount"))
        units = library.get(
            "nucleic_acid_starting_quantity_units", library.get("starting_amount_units")
        )
        if amount not in (None, ""):
            starting_quantities.append(f"{amount} {units or ''}".strip())

    strand_specificity = aggregate_fields(libraries, "strand_specificity")
    if strand_specificity in {"forward", "reverse"}:
        strand_specificity = f"Strand-specific ({strand_specificity})"

    spike_ins = aggregate_fields(
        [experiment, *libraries], "spikeins_used", "spike_ins", "spikeins"
    )
    aliases = join_values(listify(experiment.get("aliases")))
    tags = join_values(
        [
            *listify(experiment.get("internal_tags")),
            *listify(experiment.get("tags")),
        ]
    )
    attribution = "Community submission" if project.lower() == "community" else ""
    audits = audit_counts(experiment)
    lines = [
        f"Experiment summary for {spec.accession}",
        "",
        f"DOI: 10.17989/{spec.accession}",
        f"ENCODE URL: {API_ROOT}/experiments/{spec.accession}/",
    ]
    if audits:
        lines.extend([f"Audits: {audits}"])
    lines.extend(
        [
        "",
        "Summary",
        "-------",
        f"Status: {experiment.get('status', '')}",
        f"Assay: {experiment.get('assay_title') or experiment.get('assay_term_name') or ''}",
        f"Biosample summary: {experiment.get('biosample_summary', '')}",
        f"Biosample type: {biosample_ontology.get('classification', '')}",
        f"Replication type: {experiment.get('replication_type', '')}",
        f"Description: {experiment.get('description', '')}",
        f"Nucleic acid type: {aggregate_fields(libraries, 'nucleic_acid_term_name')}",
        f"Library starting quantity: {join_values(starting_quantities)}",
        f"Library construction method: {aggregate_fields(libraries, 'construction_method')}",
        f"Strand specificity: {strand_specificity}",
        f"Spike-ins datasets: {spike_ins}",
        f"Platform: {experiment_platforms(experiment)}",
        "",
        "Attribution",
        "-----------",
        f"Submission category: {attribution}",
        f"Lab: {lab}",
        f"Project: {project}",
        f"Aliases: {aliases}",
        f"Date submitted: {display_date(experiment.get('date_created') or experiment.get('date_submitted'))}",
        f"Date released: {display_date(experiment.get('date_released'))}",
        f"Tags: {tags}",
        "",
        f"{str(experiment.get('replication_type', '')).capitalize()} replicates",
        "-" * (len(str(experiment.get('replication_type', ''))) + 11),
        "Biological replicate | Technical replicate | Biosample summary | Biosample accession | Library accession",
        ]
    )
    for bio, details in sorted(rep_info.items()):
        lines.append(
            f"{bio} | {join_values(details['technical_replicates']) or 'unknown'} | "
            f"{join_values(details['biosample_summaries']) or 'unknown'} | "
            f"{join_values(details['biosamples']) or 'unknown'} | "
            f"{join_values(details['libraries']) or 'unknown'}"
        )
    lines.extend([
        "",
        "Local inventory",
        "---------------",
        f"Cell line: {spec.cell_line}",
        f"Arm: {spec.arm}",
        f"Biological replicates found: {len(rep_info)}",
        f"Associated files found: {len(rows)}",
        f"Files selected for download: {selected}",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_inventory(specs: list[ExperimentSpec], output_root: Path) -> list[dict[str, str]]:
    all_rows: list[dict[str, str]] = []
    failures: list[str] = []
    for index, spec in enumerate(specs, start=1):
        print(f"[{index}/{len(specs)}] Inventory {spec.cell_line} {spec.arm} {spec.accession}")
        try:
            all_rows.extend(inventory_experiment(spec, output_root))
        except Exception as error:
            failures.append(f"{spec.accession}: {error}")
            print(f"ERROR: {failures[-1]}", file=sys.stderr)
    write_csv(output_root / "manifest.csv", all_rows)
    if failures:
        (output_root / "inventory_failures.txt").write_text("\n".join(failures) + "\n", encoding="utf-8")
        raise RuntimeError(f"Inventory failed for {len(failures)} experiment(s)")
    failure_file = output_root / "inventory_failures.txt"
    if failure_file.exists():
        failure_file.unlink()
    return all_rows


def md5_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def download_one(
    row: dict[str, str], overwrite: bool = False, attempts: int = 6
) -> tuple[str, str]:
    url = row["download_url"]
    destination = Path(row["local_path"])
    expected_md5 = row["md5sum"].lower()
    if not url or not row["local_path"]:
        return row["file_accession"], "missing URL or destination"
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not overwrite:
        if not expected_md5 or md5_file(destination) == expected_md5:
            return row["file_accession"], "already present and valid"
    part = destination.with_suffix(destination.suffix + ".part")
    if overwrite and part.exists():
        part.unlink()
    for attempt in range(1, attempts + 1):
        try:
            start = part.stat().st_size if part.exists() else 0
            headers = {"User-Agent": USER_AGENT}
            if start:
                headers["Range"] = f"bytes={start}-"
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=300) as response:
                response_status = getattr(response, "status", 200)
                mode = "ab" if start and response_status == 206 else "wb"
                with part.open(mode) as handle:
                    shutil.copyfileobj(response, handle, length=8 * 1024 * 1024)
            if expected_md5 and md5_file(part) != expected_md5:
                part.unlink(missing_ok=True)
                raise OSError("MD5 mismatch; retrying from the beginning")
            os.replace(part, destination)
            retry_note = f" after {attempt} attempts" if attempt > 1 else ""
            return row["file_accession"], f"downloaded and valid{retry_note}"
        except urllib.error.HTTPError as error:
            if error.code == 416 and part.exists():
                if expected_md5 and md5_file(part) == expected_md5:
                    os.replace(part, destination)
                    return row["file_accession"], "resumed partial was already complete and valid"
                part.unlink(missing_ok=True)
            if attempt == attempts:
                raise
            delay = min(60, 5 * 2 ** (attempt - 1))
            print(
                f"RETRY {row['file_accession']}: attempt {attempt}/{attempts} failed "
                f"({error}); waiting {delay}s",
                flush=True,
            )
            time.sleep(delay)
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as error:
            if attempt == attempts:
                raise
            delay = min(60, 5 * 2 ** (attempt - 1))
            print(
                f"RETRY {row['file_accession']}: attempt {attempt}/{attempts} failed "
                f"({error}); waiting {delay}s",
                flush=True,
            )
            time.sleep(delay)
    raise RuntimeError("unreachable")


def download_selected(rows: list[dict[str, str]], workers: int, overwrite: bool) -> None:
    selected = [row for row in rows if row["selected"] == "yes"]
    expected_bytes = sum(int(row.get("file_size") or 0) for row in selected)
    print("=" * 88, flush=True)
    print("ENCODE DATA DOWNLOAD", flush=True)
    print(f"Selected files: {len(selected)}", flush=True)
    print(f"Expected data: {expected_bytes / 1024**3:.2f} GiB", flush=True)
    print(f"Parallel workers: {workers}", flush=True)
    print("Existing validated files will be skipped; .part files will resume.", flush=True)
    print("=" * 88, flush=True)
    failures: list[str] = []
    processed = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(download_one, row, overwrite): row for row in selected}
        for future in as_completed(futures):
            row = futures[future]
            processed += 1
            try:
                accession, status = future.result()
                if "mismatch" in status or "missing" in status:
                    failures.append(f"{accession}: {status}")
            except Exception as error:
                status = f"ERROR: {error}"
                failures.append(f"{row['file_accession']}: {error}")
            completed_count, completed_bytes, partial_count, partial_bytes = progress_totals(selected)
            present_bytes = completed_bytes + partial_bytes
            percent = (100 * present_bytes / expected_bytes) if expected_bytes else 0
            expected_file_bytes = int(row.get("file_size") or 0)
            print("-" * 88, flush=True)
            print(f"RESULT {processed}/{len(selected)}: {row['file_accession']}", flush=True)
            print(f"Status: {status}", flush=True)
            print(
                f"Sample: {row['cell_line']} | {row['arm']} | {row['experiment_accession']} | "
                f"biological replicate {row['biological_replicates'] or '?'} | "
                f"technical replicate {row['technical_replicates'] or '?'}",
                flush=True,
            )
            print(
                f"File: {row['output_type']} | {row['file_format']} | "
                f"mate {row['read_mate'] or 'not applicable'} | "
                f"expected {expected_file_bytes / 1024**2:.1f} MiB",
                flush=True,
            )
            if row.get("paired_with"):
                print(f"Paired with: {row['paired_with']}", flush=True)
            print(f"Saved as: {row['local_path']}", flush=True)
            print(
                f"OVERALL: {percent:.2f}% | validated {completed_count}/{len(selected)} files "
                f"({completed_bytes / 1024**3:.2f} GiB) | "
                f"active partials {partial_count} ({partial_bytes / 1024**3:.2f} GiB) | "
                f"failures {len(failures)}",
                flush=True,
            )
    if failures:
        print("=" * 88, file=sys.stderr, flush=True)
        print("FAILED FILES", file=sys.stderr, flush=True)
        for failure in failures:
            print(f"  {failure}", file=sys.stderr, flush=True)
        raise RuntimeError(f"Download failed for {len(failures)} file(s)")
    print("=" * 88, flush=True)
    print("ALL SELECTED FILES DOWNLOADED AND MD5-VALIDATED", flush=True)
    print("=" * 88, flush=True)


def verify_manifest(manifest: Path) -> int:
    with manifest.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    failures = 0
    for row in rows:
        if row.get("selected") != "yes":
            continue
        path = Path(row["local_path"])
        if not path.exists():
            print(f"MISSING {row['file_accession']} {path}")
            failures += 1
        elif row.get("md5sum") and md5_file(path) != row["md5sum"].lower():
            print(f"BAD_MD5 {row['file_accession']} {path}")
            failures += 1
        else:
            print(f"OK {row['file_accession']} {path}")
    return failures


def progress_totals(rows: list[dict[str, str]]) -> tuple[int, int, int, int]:
    completed_count = 0
    completed_bytes = 0
    partial_count = 0
    partial_bytes = 0
    for row in rows:
        destination = Path(row["local_path"])
        part = Path(f"{row['local_path']}.part")
        if destination.is_file():
            completed_count += 1
            completed_bytes += destination.stat().st_size
        elif part.is_file():
            partial_count += 1
            partial_bytes += part.stat().st_size
    return completed_count, completed_bytes, partial_count, partial_bytes


def show_status(manifest: Path) -> int:
    with manifest.open(newline="", encoding="utf-8-sig") as handle:
        selected = [row for row in csv.DictReader(handle) if row.get("selected") == "yes"]
    completed = []
    partial = []
    for row in selected:
        destination = Path(row["local_path"])
        part = Path(f"{row['local_path']}.part")
        if destination.is_file():
            completed.append((row, destination.stat().st_size))
        elif part.is_file():
            partial.append((row, part.stat().st_size))

    expected_bytes = sum(int(row.get("file_size") or 0) for row in selected)
    completed_bytes = sum(size for _, size in completed)
    partial_bytes = sum(size for _, size in partial)
    present_bytes = completed_bytes + partial_bytes
    percent = (100 * present_bytes / expected_bytes) if expected_bytes else 0
    print(f"Progress: {percent:.2f}%")
    print(f"Completed files: {len(completed)} / {len(selected)}")
    print(f"Validated data: {completed_bytes / 1024**3:.2f} GiB")
    print(f"Partial data: {partial_bytes / 1024**3:.2f} GiB")
    print(f"Total expected: {expected_bytes / 1024**3:.2f} GiB")
    print(f"Partial/resumable files: {len(partial)}")
    for row, size in sorted(partial, key=lambda item: item[1], reverse=True):
        print(
            f"  {row['file_accession']} | {row['arm']} {row['experiment_accession']} | "
            f"replicate {row['biological_replicates'] or '?'} | "
            f"{row['output_type']} | {size / 1024**2:.1f} MiB"
        )
    if len(completed) == len(selected) and not partial:
        print("STATUS: DOWNLOAD COMPLETE")
        return 0
    print("STATUS: DOWNLOAD INCOMPLETE OR IN PROGRESS")
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("inventory", "download", "verify", "status"))
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.workers < 1:
        raise ValueError("--workers must be at least 1")
    if args.mode == "status":
        manifest = args.output / "manifest.csv"
        if not manifest.exists():
            raise FileNotFoundError(f"Run inventory first; missing {manifest}")
        return show_status(manifest)
    if args.mode == "verify":
        manifest = args.output / "manifest.csv"
        if not manifest.exists():
            raise FileNotFoundError(f"Run inventory first; missing {manifest}")
        return 1 if verify_manifest(manifest) else 0
    specs = load_specs(args.config)
    rows = build_inventory(specs, args.output)
    if args.mode == "download":
        download_selected(rows, args.workers, args.overwrite)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted; partial downloads were retained for resume", file=sys.stderr)
        raise SystemExit(130)
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
