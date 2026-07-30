"""Command-line interface for the v0.5.0 Splice Provenance Layer."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from neoag.splice.coordinates import file_sha256
from neoag.utils import read_tsv, write_json

from .pipeline import SpliceLayer, build_splice_provenance_layer
from .schemas import OUTPUT_FILENAMES


def _tool_version(value: str) -> tuple[str, str]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("tool version must use TOOL=VERSION")
    tool, version = value.split("=", 1)
    return tool.strip(), version.strip()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="NeoAg v0.5.0 formal Splice Provenance Layer")
    sub = parser.add_subparsers(dest="command", required=True)

    build = sub.add_parser("build", help="Build junction→event→transcript→ORF→peptide provenance tables")
    build.add_argument("--sample-id", required=True)
    build.add_argument("--outdir", required=True)
    build.add_argument("--genome-build", default="GRCh38")
    build.add_argument("--disease-profile", default="default")
    build.add_argument("--junctions")
    build.add_argument("--junction-coordinate-system", default="auto")
    build.add_argument("--junction-source-assay-id", default="")
    build.add_argument("--star-junctions")
    build.add_argument("--star-junction-source-assay-id", default="")
    build.add_argument("--spladder-gff3", action="append", default=[])
    build.add_argument("--spladder-txt", action="append", default=[])
    build.add_argument("--irfinder", action="append", default=[])
    build.add_argument("--irfinder-coordinate-system", default="UNSPECIFIED")
    build.add_argument("--immunopepper-meta", action="append", default=[])
    build.add_argument("--immunopepper-kmers", action="append", default=[])
    build.add_argument("--pvacbind", action="append", default=[])
    build.add_argument("--pvacbind-fasta-map")
    build.add_argument("--normal-junctions", action="append", default=[])
    build.add_argument("--normal-coordinate-system", default="auto")
    build.add_argument("--normal-coverage", action="append", default=[])
    build.add_argument("--high-order-evidence", action="append", default=[])
    build.add_argument("--tool-version", action="append", type=_tool_version, default=[])
    build.add_argument("--strict", action="store_true")

    fasta = sub.add_parser("write-pvacbind-fasta", help="Regenerate pVACbind FASTA and exact index map from formal ORFs")
    fasta.add_argument("--layer-dir", required=True)
    fasta.add_argument("--sample-id", required=True)
    fasta.add_argument("--outdir")

    validate = sub.add_parser("validate", help="Validate referential integrity and manifest hashes")
    validate.add_argument("--layer-dir", required=True)
    validate.add_argument("--strict", action="store_true")
    validate.add_argument("--report")
    return parser


def _load_layer(layer_dir: Path, sample_id: str) -> SpliceLayer:
    layer = SpliceLayer(sample_id=sample_id)
    for table, filename in OUTPUT_FILENAMES.items():
        path = layer_dir / filename
        if table in {"manifest", "pvacbind_fasta", "raw_events", "raw_peptides", "rna_junction_evidence", "qc"}:
            continue
        if path.is_file() and path.stat().st_size:
            layer.tables[table] = read_tsv(path)
    return layer


def _validate(layer_dir: Path) -> dict[str, object]:
    manifest_path = layer_dir / OUTPUT_FILENAMES["manifest"]
    errors: list[str] = []
    warnings: list[str] = []
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else {}
    for key, meta in manifest.get("outputs", {}).items():
        path = Path(meta.get("path", ""))
        if not path.is_absolute():
            path = layer_dir / path.name
        if not path.is_file():
            errors.append(f"missing output: {key}: {path}")
        elif meta.get("sha256") and file_sha256(path) != meta["sha256"]:
            errors.append(f"sha256 mismatch: {key}: {path}")
    qc_path = layer_dir / OUTPUT_FILENAMES["qc"]
    qc = read_tsv(qc_path) if qc_path.is_file() else []
    errors.extend(f"QC failure: {row.get('metric')}={row.get('value')}" for row in qc if row.get("status") == "FAIL")
    if not manifest:
        warnings.append("provenance_manifest.json is missing")
    return {"status": "PASS" if not errors else "FAIL", "errors": errors, "warnings": warnings, "layer_dir": str(layer_dir)}


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "build":
        outputs = build_splice_provenance_layer(
            sample_id=args.sample_id, outdir=args.outdir, genome_build=args.genome_build,
            disease_profile=args.disease_profile, junctions=args.junctions,
            junction_coordinate_system=args.junction_coordinate_system,
            junction_source_assay_id=args.junction_source_assay_id,
            star_junctions=args.star_junctions,
            star_junction_source_assay_id=args.star_junction_source_assay_id,
            spladder_gff3=args.spladder_gff3, spladder_txt=args.spladder_txt,
            irfinder=args.irfinder, irfinder_coordinate_system=args.irfinder_coordinate_system,
            immunopepper_meta=args.immunopepper_meta, immunopepper_kmers=args.immunopepper_kmers,
            pvacbind=args.pvacbind, pvacbind_fasta_map=args.pvacbind_fasta_map,
            normal_junctions=args.normal_junctions, normal_coordinate_system=args.normal_coordinate_system,
            normal_coverage=args.normal_coverage, high_order_evidence=args.high_order_evidence,
            tool_versions=dict(args.tool_version), strict=args.strict,
        )
        print(json.dumps({key: str(path) for key, path in outputs.items()}, indent=2, ensure_ascii=False))
        return 0
    if args.command == "write-pvacbind-fasta":
        layer_dir = Path(args.layer_dir)
        layer = _load_layer(layer_dir, args.sample_id)
        fasta, mapping = layer.write_pvacbind_fasta(args.outdir or layer_dir)
        print(json.dumps({"fasta": str(fasta), "map": str(mapping)}, indent=2))
        return 0
    if args.command == "validate":
        report = _validate(Path(args.layer_dir))
        if args.report:
            write_json(args.report, report)
        print(json.dumps(report, indent=2, ensure_ascii=False))
        return 1 if args.strict and report["status"] != "PASS" else 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
