from pathlib import Path
import csv

from .parsers import parse_fasta_ids


def not_empty(x):
    if x is None:
        return False
    x = str(x).strip()
    return x not in ["", "-", "NA", "NaN", "nan", "None", "none", "no_annotation"]


def read_table_by_gene(path):
    data = {}
    path = Path(path)

    if not path.exists():
        return data

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gid = row.get("gene_id", "")
            if gid:
                data[gid] = row
    return data


FUNCTION_RULES = [
    ("transcription_regulation", ["transcription factor", "zinc finger", "homeobox", "bhlh", "dna-binding", "gata", "runx", "ets", "sox"]),
    ("signal_transduction", ["kinase", "phosphatase", "receptor", "signaling", "mapk", "wnt", "notch", "tgf", "jak", "stat"]),
    ("metabolism", ["metabolism", "metabolic", "dehydrogenase", "oxidase", "transferase", "synthetase", "synthase", "enzyme"]),
    ("transport", ["transporter", "transport", "channel", "solute carrier"]),
    ("immune_phagocyte", ["immune", "complement", "lectin", "cathepsin", "lysozyme", "phagocyt", "scavenger", "macrophage", "c1q"]),
    ("cell_cycle", ["cell cycle", "cyclin", "mitotic", "mitosis", "spindle", "dna replication"]),
    ("cytoskeleton_ecm_adhesion", ["actin", "tubulin", "collagen", "cadherin", "integrin", "extracellular matrix", "adhesion"]),
    ("protein_processing", ["ribosomal", "ribosome", "translation", "proteasome", "ubiquitin", "chaperone"]),
    ("rna_processing", ["rna-binding", "splice", "splicing", "helicase"]),
]

DOMAIN_RULES = [
    ("zinc_finger_or_DNA_binding", ["zinc finger", "gata", "homeobox", "hmg", "ets", "runt", "bhlh", "forkhead"]),
    ("protein_kinase", ["protein kinase", "tyrosine kinase", "pkinase"]),
    ("receptor_or_membrane", ["7tm", "gpcr", "receptor", "transmembrane", "ig-like", "immunoglobulin"]),
    ("enzyme_domain", ["dehydrogenase", "transferase", "hydrolase", "oxidase", "peptidase", "protease", "catalytic"]),
    ("rna_dna_processing", ["helicase", "rna recognition", "rrm", "polymerase", "nuclease"]),
    ("protein_interaction_repeat", ["wd40", "ankyrin", "leucine rich repeat", "tpr"]),
    ("ecm_adhesion", ["cadherin", "integrin", "collagen", "fibronectin", "egf-like"]),
]


def classify(text, rules):
    text = text.lower()
    hits = []
    for name, kws in rules:
        if any(kw.lower() in text for kw in kws):
            hits.append(name)

    if not hits:
        return "other_or_unknown", "other_or_unknown"

    return hits[0], ",".join(hits)


def make_annotation_text(row):
    parts = []
    for c in ["Preferred_name", "Description", "KEGG_Pathway", "BRITE", "PFAMs", "pfam_domains_hmmscan"]:
        if not_empty(row.get(c)):
            parts.append(str(row.get(c)))
    return " | ".join(parts)


def add_integrated_fields(row):
    annotation_text = make_annotation_text(row)
    row["annotation_text"] = annotation_text

    f_primary, f_all = classify(annotation_text, FUNCTION_RULES)
    row["function_class_primary"] = f_primary
    row["function_class_all"] = f_all
    row["matched_category_count"] = str(0 if f_all == "other_or_unknown" else len(f_all.split(",")))

    domain_text = " | ".join([str(row.get("PFAMs", "")), str(row.get("pfam_domains_hmmscan", ""))])
    d_primary, d_all = classify(domain_text, DOMAIN_RULES)
    row["domain_class_primary"] = d_primary
    row["domain_class_all"] = d_all

    has_diamond = not_empty(row.get("diamond_best_hit"))
    has_eggnog = any(not_empty(row.get(c)) for c in ["Preferred_name", "Description", "eggNOG_OGs", "seed_ortholog"])
    has_pfam = any(not_empty(row.get(c)) for c in ["PFAMs", "pfam_domains_hmmscan"])

    if has_diamond and has_eggnog and has_pfam:
        level = "diamond_eggnog_pfam_supported"
    elif has_eggnog and has_pfam:
        level = "eggnog_and_pfam_supported"
    elif has_diamond and has_eggnog:
        level = "diamond_eggnog_supported"
    elif has_diamond and has_pfam:
        level = "diamond_pfam_supported"
    elif has_eggnog:
        level = "eggnog_supported_only"
    elif has_pfam:
        level = "pfam_supported_only"
    elif has_diamond:
        level = "diamond_supported_only"
    else:
        level = "no_annotation"

    row["annotation_evidence_level"] = level
    return row


def merge_annotations(protein_fasta, diamond_file, eggnog_file, pfam_file, outfile):
    gene_ids = parse_fasta_ids(protein_fasta)
    diamond = read_table_by_gene(diamond_file)
    eggnog = read_table_by_gene(eggnog_file)
    pfam = read_table_by_gene(pfam_file)

    cols = [
        "gene_id",
        "diamond_best_hit", "diamond_pident", "diamond_align_len", "diamond_evalue",
        "diamond_bitscore", "diamond_qcovhsp", "diamond_scovhsp",
        "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
        "COG_category", "Description", "Preferred_name",
        "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
        "KEGG_Reaction", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs",
        "annotation_text", "function_class_primary", "function_class_all", "matched_category_count",
        "n_pfam_domains", "pfam_domains_hmmscan",
        "domain_class_primary", "domain_class_all", "annotation_evidence_level",
    ]

    rows = []
    for gid in gene_ids:
        row = {c: "" for c in cols}
        row["gene_id"] = gid

        if gid in diamond:
            row.update(diamond[gid])
        if gid in eggnog:
            row.update(eggnog[gid])
        if gid in pfam:
            row.update(pfam[gid])

        rows.append(add_integrated_fields(row))

    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    with open(outfile, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] Final merged annotation: {len(rows)} genes -> {outfile}")
    return outfile
