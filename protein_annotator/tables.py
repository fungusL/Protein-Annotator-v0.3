from pathlib import Path
import csv
import re


FUNC_COLS = [
    "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
    "KEGG_Reaction", "BRITE", "KEGG_TC", "CAZy",
    "BiGG_Reaction", "COG_category", "annotation_text",
    "function_class_primary", "function_class_all",
]

MARKER_KEYWORDS = {
    "endothelial_vascular": [
        "fli", "erg", "ets", "etv", "kdr", "flk", "flt1", "vegfr",
        "tek", "tie", "pecam", "cdh5", "cadherin 5", "vwf", "egfl7",
        "endothelial", "vascular", "pdgfr",
    ],
    "hematopoietic_progenitor": [
        "tal1", "scl", "gata", "runx", "lmo2", "myb", "meis",
        "hhex", "hematopoietic", "hemogenic",
    ],
    "myeloid_phagocyte": [
        "spi1", "pu.1", "csf1r", "mpeg", "lyz", "lysozyme",
        "ctss", "cathepsin", "c1q", "complement", "mrc", "marco",
        "scavenger", "lectin", "phagocyt",
    ],
    "erythroid": [
        "gata1", "klf", "alas2", "tfrc", "globin", "hemoglobin",
        "hba", "hbb", "erythroid",
    ],
    "tf_grn_core": [
        "transcription factor", "zinc finger", "homeobox", "bhlh",
        "runx", "gata", "ets", "sox", "irf", "cebpa", "cebpb",
        "nuclear receptor", "forkhead", "dna-binding",
    ],
}

SEARCH_COLS = [
    "gene_id", "Preferred_name", "Description",
    "annotation_text", "function_class_primary", "function_class_all",
    "PFAMs", "pfam_domains_hmmscan", "domain_class_primary",
    "domain_class_all", "diamond_best_hit",
]

PREFERRED_COLS = [
    "gene_id",
    "multi_evidence_class",
    "evidence_score",
    "recommended_for_marker",
    "diamond_best_hit",
    "diamond_pident",
    "diamond_align_len",
    "diamond_evalue",
    "diamond_bitscore",
    "diamond_qcovhsp",
    "diamond_scovhsp",
    "seed_ortholog",
    "evalue",
    "score",
    "eggNOG_OGs",
    "max_annot_lvl",
    "COG_category",
    "Description",
    "Preferred_name",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "BRITE",
    "PFAMs",
    "annotation_text",
    "function_class_primary",
    "function_class_all",
    "matched_category_count",
    "n_pfam_domains",
    "pfam_domains_hmmscan",
    "domain_class_primary",
    "domain_class_all",
    "annotation_evidence_level",
    "diamond_supported",
    "diamond_strong_supported",
    "eggnog_supported",
    "eggnog_strong_supported",
    "pfam_supported",
    "functional_supported",
    "is_endothelial_vascular_candidate",
    "is_hematopoietic_progenitor_candidate",
    "is_myeloid_phagocyte_candidate",
    "is_erythroid_candidate",
    "is_tf_grn_core_candidate",
    "is_key_marker_candidate",
]

KEYINFO_COLS = [
    "gene_id",
    "Preferred_name",
    "Description",
    "diamond_best_hit",
    "seed_ortholog",
    "eggNOG_OGs",
    "max_annot_lvl",
    "diamond_pident",
    "diamond_qcovhsp",
    "diamond_scovhsp",
    "diamond_evalue",
    "diamond_bitscore",
    "annotation_text",
    "function_class_primary",
    "function_class_all",
    "COG_category",
    "GOs",
    "KEGG_ko",
    "KEGG_Pathway",
    "BRITE",
    "PFAMs",
    "pfam_domains_hmmscan",
    "n_pfam_domains",
    "domain_class_primary",
    "domain_class_all",
    "multi_evidence_class",
    "evidence_score",
    "annotation_evidence_level",
]


def not_empty(x):
    if x is None:
        return False
    x = str(x).strip()
    return x not in ["", "-", "NA", "NaN", "nan", "None", "none", "no_annotation"]


def to_float(x):
    try:
        if not_empty(x):
            return float(str(x).replace(",", ""))
    except Exception:
        pass
    return None


def bstr(x):
    return "TRUE" if x else "FALSE"


def add_evidence_columns(row):
    diamond_has_hit = not_empty(row.get("diamond_best_hit"))

    pident = to_float(row.get("diamond_pident"))
    align_len = to_float(row.get("diamond_align_len"))
    evalue = to_float(row.get("diamond_evalue"))
    bitscore = to_float(row.get("diamond_bitscore"))
    qcov = to_float(row.get("diamond_qcovhsp"))
    scov = to_float(row.get("diamond_scovhsp"))

    diamond_supported = (
        diamond_has_hit and
        evalue is not None and evalue <= 1e-20 and
        bitscore is not None and bitscore >= 80 and
        align_len is not None and align_len >= 80 and
        qcov is not None and qcov >= 50
    )

    diamond_strong_supported = (
        diamond_supported and
        pident is not None and pident >= 30 and
        scov is not None and scov >= 50
    )

    eggnog_has_og = not_empty(row.get("eggNOG_OGs"))
    eggnog_has_name = not_empty(row.get("Preferred_name"))
    eggnog_has_desc = not_empty(row.get("Description"))
    eggnog_has_seed = not_empty(row.get("seed_ortholog"))

    eggnog_supported = eggnog_has_og or eggnog_has_name or eggnog_has_desc or eggnog_has_seed
    eggnog_strong_supported = (
        eggnog_supported and
        (eggnog_has_name or eggnog_has_desc) and
        (eggnog_has_og or eggnog_has_seed)
    )

    pfam_supported = (
        not_empty(row.get("PFAMs")) or
        not_empty(row.get("pfam_domains_hmmscan")) or
        ((to_float(row.get("n_pfam_domains")) or 0) > 0)
    )

    functional_supported = any(not_empty(row.get(c)) for c in FUNC_COLS)

    evidence_score = 0
    for flag in [
        diamond_supported,
        diamond_strong_supported,
        eggnog_supported,
        eggnog_strong_supported,
        pfam_supported,
        functional_supported,
    ]:
        if flag:
            evidence_score += 1

    level = str(row.get("annotation_evidence_level", "")).lower()
    if re.search(r"high|高|eggnog_and_pfam_supported|eggnog.*pfam|diamond_eggnog_pfam", level):
        evidence_score += 1

    if diamond_supported and eggnog_supported and pfam_supported:
        multi_class = "strict_high_confidence"
    elif eggnog_supported and pfam_supported:
        multi_class = "eggnog_pfam_supported"
    elif diamond_supported and eggnog_supported:
        multi_class = "diamond_eggnog_supported"
    elif diamond_supported and pfam_supported:
        multi_class = "diamond_pfam_supported"
    elif diamond_supported or eggnog_supported or pfam_supported:
        multi_class = "single_or_weak_evidence"
    else:
        multi_class = "no_annotation"

    recommended_for_marker = (
        multi_class == "strict_high_confidence" or
        (multi_class == "eggnog_pfam_supported" and functional_supported)
    )

    search_text = " ".join(str(row.get(c, "")) for c in SEARCH_COLS).lower()
    is_key_marker = False

    for group, kws in MARKER_KEYWORDS.items():
        flag = any(kw.lower() in search_text for kw in kws)
        row[f"is_{group}_candidate"] = bstr(flag)
        is_key_marker = is_key_marker or flag

    row["diamond_supported"] = bstr(diamond_supported)
    row["diamond_strong_supported"] = bstr(diamond_strong_supported)
    row["eggnog_supported"] = bstr(eggnog_supported)
    row["eggnog_strong_supported"] = bstr(eggnog_strong_supported)
    row["pfam_supported"] = bstr(pfam_supported)
    row["functional_supported"] = bstr(functional_supported)

    row["evidence_score"] = str(evidence_score)
    row["multi_evidence_class"] = multi_class
    row["recommended_for_marker"] = bstr(recommended_for_marker)
    row["is_key_marker_candidate"] = bstr(is_key_marker)

    return row


def write_tsv(outfile, rows, cols):
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    with open(outfile, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def filter_annotation_tables(annotation_tsv, outdir, prefix):
    annotation_tsv = Path(annotation_tsv)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rows = []

    with open(annotation_tsv, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(add_evidence_columns(row))

    if not rows:
        raise RuntimeError(f"No rows found in {annotation_tsv}")

    out_cols = [c for c in PREFERRED_COLS if c in rows[0]]

    subsets = {
        "all_with_scores": rows,
        "strict_high_confidence": [r for r in rows if r["multi_evidence_class"] == "strict_high_confidence"],
        "eggnog_pfam_supported": [r for r in rows if r["multi_evidence_class"] == "eggnog_pfam_supported"],
        "diamond_eggnog_supported": [r for r in rows if r["multi_evidence_class"] == "diamond_eggnog_supported"],
        "diamond_pfam_supported": [r for r in rows if r["multi_evidence_class"] == "diamond_pfam_supported"],
        "key_marker_high_conf": [
            r for r in rows
            if r["recommended_for_marker"] == "TRUE" and r["is_key_marker_candidate"] == "TRUE"
        ],
        "key_marker_manual_check": [
            r for r in rows
            if r["is_key_marker_candidate"] == "TRUE"
            and r["recommended_for_marker"] != "TRUE"
            and r["multi_evidence_class"] != "no_annotation"
        ],
        "no_annotation": [r for r in rows if r["multi_evidence_class"] == "no_annotation"],
    }

    result_files = {}

    for name, sub in subsets.items():
        outfile = outdir / f"{prefix}.{name}.tsv"
        write_tsv(outfile, sub, out_cols)
        result_files[name] = str(outfile)
        print(f"{name}: {len(sub)} -> {outfile}")

    keyinfo_cols = [c for c in KEYINFO_COLS if c in rows[0]]
    keyinfo = outdir / f"{prefix}.strict_high_confidence.keyinfo.tsv"
    write_tsv(keyinfo, subsets["strict_high_confidence"], keyinfo_cols)
    result_files["strict_high_confidence_keyinfo"] = str(keyinfo)
    print(f"strict_high_confidence.keyinfo: {len(subsets['strict_high_confidence'])} -> {keyinfo}")

    summary = {}
    for r in rows:
        summary[r["multi_evidence_class"]] = summary.get(r["multi_evidence_class"], 0) + 1

    summary_file = outdir / f"{prefix}.summary.tsv"
    with open(summary_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["multi_evidence_class", "n_genes"])
        for k, v in sorted(summary.items(), key=lambda x: x[1], reverse=True):
            writer.writerow([k, v])

    print("\n===== Summary =====")
    for k, v in sorted(summary.items(), key=lambda x: x[1], reverse=True):
        print(f"{k}\t{v}")

    result_files["summary"] = str(summary_file)
    return result_files
