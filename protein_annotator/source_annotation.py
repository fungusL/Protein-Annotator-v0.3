#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


UNINFORMATIVE_PATTERNS = [
    r"\bhypothetical protein\b",
    r"\buncharacterized protein\b",
    r"\bpredicted protein\b",
    r"\blow quality protein\b",
    r"\bunnamed protein\b",
    r"\bunknown protein\b",
    r"\bnovel protein\b",
]


def normalize_text(x):
    if x is None:
        return ""
    x = str(x).lower()
    x = re.sub(r"\s+", " ", x).strip()
    return x


def is_uninformative_product(product):
    text = normalize_text(product)
    if not text:
        return True
    return any(re.search(pat, text) for pat in UNINFORMATIVE_PATTERNS)


def parse_fasta_headers(protein_fasta):
    protein_fasta = Path(protein_fasta)
    header_info = {}

    with protein_fasta.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith(">"):
                continue

            header = line[1:].strip()
            if not header:
                continue

            first_token = header.split()[0]
            rest = header[len(first_token):].strip()

            source_protein_id = first_token
            source_product = rest
            source_species = ""

            # UniProt style:
            # >sp|Q13591|SEM5A_HUMAN Semaphorin-5A OS=Homo sapiens ...
            if "|" in first_token and first_token.count("|") >= 2:
                parts = first_token.split("|")
                if len(parts) >= 3:
                    source_protein_id = parts[1]

                m_os = re.search(r"\bOS=([^=]+?)(?:\s+[A-Z]{2}=|$)", rest)
                if m_os:
                    source_species = m_os.group(1).strip()
                    source_product = rest[:m_os.start()].strip()

            # NCBI style:
            # >XP_035657343.1 semaphorin-5A-like, partial [Branchiostoma floridae]
            else:
                m_species = re.search(r"\s*\[([^\]]+)\]\s*$", source_product)
                if m_species:
                    source_species = m_species.group(1).strip()
                    source_product = source_product[:m_species.start()].strip()

            source_product = source_product.strip()

            info = {
                "source_protein_id": source_protein_id,
                "source_product": source_product,
                "source_species": source_species,
                "source_annotation_informative": "no" if is_uninformative_product(source_product) else "yes",
                "source_header": header,
            }

            header_info[first_token] = info
            header_info[source_protein_id] = info

    return header_info


def read_tsv(path):
    path = Path(path)
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames or []
        rows = [dict(r) for r in reader]
    return fieldnames, rows


def write_tsv(path, fieldnames, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def add_source_annotation(protein_fasta, annotation_tsv, out_tsv, gene_id_col="gene_id"):
    source_map = parse_fasta_headers(protein_fasta)
    fieldnames, rows = read_tsv(annotation_tsv)

    if not fieldnames:
        raise ValueError(f"No header found in annotation table: {annotation_tsv}")

    if gene_id_col not in fieldnames:
        raise ValueError(
            f"Column '{gene_id_col}' not found in annotation table. "
            f"Available columns: {', '.join(fieldnames)}"
        )

    source_cols = [
        "source_protein_id",
        "source_product",
        "source_species",
        "source_annotation_informative",
        "source_header",
    ]

    new_fieldnames = []
    for col in fieldnames:
        new_fieldnames.append(col)
        if col == gene_id_col:
            for sc in source_cols:
                if sc not in new_fieldnames:
                    new_fieldnames.append(sc)

    for sc in source_cols:
        if sc not in new_fieldnames:
            new_fieldnames.append(sc)

    total = 0
    matched = 0
    informative = 0
    new_rows = []

    for row in rows:
        total += 1
        gid = row.get(gene_id_col, "")

        info = source_map.get(gid)

        if info is None and gid:
            gid_no_version = gid.split(".")[0]
            info = source_map.get(gid_no_version)

        if info is None:
            info = {
                "source_protein_id": "",
                "source_product": "",
                "source_species": "",
                "source_annotation_informative": "no",
                "source_header": "",
            }
        else:
            matched += 1
            if info.get("source_annotation_informative") == "yes":
                informative += 1

        for k, v in info.items():
            row[k] = v

        new_rows.append(row)

    write_tsv(out_tsv, new_fieldnames, new_rows)

    return {
        "total_rows": total,
        "matched_source_headers": matched,
        "unmatched_source_headers": total - matched,
        "informative_source_products": informative,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Add source protein FASTA header annotation to an existing annotation TSV."
    )
    parser.add_argument("--protein-fasta", required=True)
    parser.add_argument("--annotation-tsv", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--gene-id-col", default="gene_id")

    args = parser.parse_args()

    stats = add_source_annotation(
        protein_fasta=args.protein_fasta,
        annotation_tsv=args.annotation_tsv,
        out_tsv=args.out_tsv,
        gene_id_col=args.gene_id_col,
    )

    print("===== Add source annotation finished =====")
    print(f"Input annotation: {args.annotation_tsv}")
    print(f"Input FASTA:      {args.protein_fasta}")
    print(f"Output TSV:       {args.out_tsv}")
    print()
    for k, v in stats.items():
        print(f"{k}: {v}")


if __name__ == "__main__":
    main()
