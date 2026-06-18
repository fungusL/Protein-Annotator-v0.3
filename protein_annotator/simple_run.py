#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from collections import Counter
from pathlib import Path

from protein_annotator.source_annotation import add_source_annotation
from protein_annotator.scoring import enrich_annotation_table


def run_cmd(cmd):
    print("[CMD]", " ".join(map(str, cmd)), flush=True)
    subprocess.run(cmd, check=True)


def export_subset(final_out: Path, out_tsv: Path, predicate) -> int:
    count = 0
    with final_out.open("r", encoding="utf-8", errors="ignore", newline="") as fin, \
         out_tsv.open("w", encoding="utf-8", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        fieldnames = reader.fieldnames or []

        writer = csv.DictWriter(
            fout,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()

        for row in reader:
            if predicate(row):
                writer.writerow(row)
                count += 1

    return count


def append_v03_summary(final_out: Path, summary_out: Path, high_conf_count: int, review_count: int) -> None:
    counters = {
        "PACS_annotation_confidence_level": Counter(),
        "EGCS_graph_status": Counter(),
        "PAI_paralog_resolution_level": Counter(),
        "DACS_domain_architecture_status": Counter(),
        "ACAD_annotation_consistency_status": Counter(),
        "WSCC_calibrated_confidence_level": Counter(),
        "CFNR_recommended_name_level": Counter(),
        "manual_review_priority": Counter(),
    }

    total = 0
    with final_out.open("r", encoding="utf-8", errors="ignore", newline="") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        for row in reader:
            total += 1
            for col, counter in counters.items():
                counter[row.get(col, "") or "NA"] += 1

    with summary_out.open("a", encoding="utf-8") as f:
        f.write("\n")
        f.write("v0.3_metric\tvalue\n")
        f.write(f"total_rows\t{total}\n")
        f.write(f"high_confidence_rows\t{high_conf_count}\n")
        f.write(f"review_needed_rows\t{review_count}\n")

        for col, counter in counters.items():
            f.write("\n")
            f.write(f"{col}\tcount\n")
            for k, v in counter.most_common():
                f.write(f"{k}\t{v}\n")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Run the full protein annotation workflow and generate simplified "
            "annotation, high-confidence, review-needed and summary tables."
        )
    )

    parser.add_argument("--protein-fasta", required=True, help="Input protein FASTA file.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--prefix", required=True, help="Output prefix.")

    parser.add_argument("--diamond-db", required=True, help="DIAMOND database.")
    parser.add_argument("--pfam-hmm", required=True, help="Pfam-A.hmm file.")
    parser.add_argument("--eggnog-data-dir", required=True, help="eggNOG-mapper data directory.")

    parser.add_argument("--threads", type=int, default=8, help="Number of threads. Default: 8.")

    parser.add_argument(
        "--no-source-annotation",
        action="store_true",
        help="Do not add original FASTA header/product annotation to the final table.",
    )

    parser.add_argument(
        "--keep-work",
        action="store_true",
        help="Keep intermediate workflow directory for debugging.",
    )

    args = parser.parse_args()

    protein_fasta = Path(args.protein_fasta).resolve()
    outdir = Path(args.outdir).resolve()
    workdir = outdir / "_work"

    outdir.mkdir(parents=True, exist_ok=True)

    if workdir.exists():
        print(f"[INFO] Removing old work directory: {workdir}", flush=True)
        shutil.rmtree(workdir)

    final_out = outdir / f"{args.prefix}.annotation.tsv"
    high_conf_out = outdir / f"{args.prefix}.high_confidence.tsv"
    review_out = outdir / f"{args.prefix}.review_needed.tsv"
    summary_out = outdir / f"{args.prefix}.summary.tsv"

    cmd = [
        "protein-annotate",
        "run",
        "--protein-fasta",
        str(protein_fasta),
        "--outdir",
        str(workdir),
        "--prefix",
        args.prefix,
        "--diamond-db",
        args.diamond_db,
        "--pfam-hmm",
        args.pfam_hmm,
        "--eggnog-data-dir",
        args.eggnog_data_dir,
        "--threads",
        str(args.threads),
    ]

    run_cmd(cmd)

    merged_tsv = workdir / "04_merged" / f"{args.prefix}.final_merged_annotation.tsv"
    summary_tsv = workdir / "05_high_confidence" / f"{args.prefix}.summary.tsv"

    if not merged_tsv.exists():
        raise FileNotFoundError(f"Final merged annotation not found: {merged_tsv}")

    if args.no_source_annotation:
        shutil.copy2(merged_tsv, final_out)
    else:
        stats = add_source_annotation(
            protein_fasta=protein_fasta,
            annotation_tsv=merged_tsv,
            out_tsv=final_out,
            gene_id_col="gene_id",
        )
        print("[INFO] Source annotation added:", flush=True)
        for k, v in stats.items():
            print(f"  {k}: {v}", flush=True)

    n_scored = enrich_annotation_table(final_out, final_out)
    print(f"[INFO] v0.3 EGCS/PAI/WSCC scoring added to {n_scored} rows", flush=True)

    high_conf_count = export_subset(
        final_out,
        high_conf_out,
        lambda row: (
            row.get("annotation_evidence_level") == "diamond_eggnog_pfam_supported"
            and row.get("WSCC_calibrated_confidence_level") in {"high_confidence", "very_high_confidence"}
            and float(row.get("PAI_paralog_ambiguity_index", "1") or "1") <= 0.60
            and row.get("manual_review_priority") != "high"
        ),
    )

    review_count = export_subset(
        final_out,
        review_out,
        lambda row: row.get("manual_review_priority") in {"high", "medium"},
    )

    print(f"[INFO] High-confidence annotations written: {high_conf_out}", flush=True)
    print(f"[INFO] High-confidence count: {high_conf_count}", flush=True)
    print(f"[INFO] Review-needed annotations written: {review_out}", flush=True)
    print(f"[INFO] Review-needed count: {review_count}", flush=True)

    if summary_tsv.exists():
        shutil.copy2(summary_tsv, summary_out)
    else:
        with summary_out.open("w", encoding="utf-8") as f:
            f.write("metric\tvalue\n")
            f.write("warning\tsummary table was not generated by the internal workflow\n")

    append_v03_summary(final_out, summary_out, high_conf_count, review_count)

    if not args.keep_work:
        print(f"[INFO] Removing intermediate work directory: {workdir}", flush=True)
        shutil.rmtree(workdir)

    print()
    print("===== Simple annotation finished =====")
    print(f"Final annotation table: {final_out}")
    print(f"High-confidence table:  {high_conf_out}")
    print(f"Review-needed table:    {review_out}")
    print(f"Summary table:          {summary_out}")
    if args.keep_work:
        print(f"Intermediate files:     {workdir}")
    else:
        print("Intermediate files:     removed")
    print()


if __name__ == "__main__":
    main()
