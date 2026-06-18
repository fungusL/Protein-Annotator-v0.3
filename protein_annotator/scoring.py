#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


BAD_VALUES = {"", "-", "NA", "na", "nan", "None", "none"}

UNINFORMATIVE_PATTERNS = [
    r"\bhypothetical protein\b",
    r"\buncharacterized protein\b",
    r"\bpredicted protein\b",
    r"\blow quality protein\b",
    r"\bunnamed protein\b",
    r"\bunknown protein\b",
    r"\bnovel protein\b",
]

STOPWORDS = {
    "protein", "predicted", "putative", "hypothetical", "uncharacterized",
    "like", "isoform", "partial", "low", "quality", "chain",
    "precursor", "fragment", "family", "member", "homolog", "domain",
    "containing", "branchiostoma", "floridae", "probable", "possible",
    "related", "type", "subunit",
}


NEW_FIELDS = [
    "PACS_annotation_confidence_score",
    "PACS_annotation_confidence_level",

    "EGCS_evidence_graph_score",
    "EGCS_evidence_graph_density",
    "EGCS_supporting_edges",
    "EGCS_conflicting_edges",
    "EGCS_graph_status",

    "PAI_paralog_ambiguity_index",
    "PAI_paralog_resolution_level",
    "PAI_detected_families",
    "PAI_reason",

    "DACS_domain_architecture_status",

    "CFNR_recommended_name",
    "CFNR_recommended_name_level",
    "CFNR_recommended_name_reason",

    "ACAD_annotation_consistency_status",

    "WSCC_pseudo_label",
    "WSCC_calibration_bin",
    "WSCC_bin_positive_rate",
    "WSCC_calibrated_confidence_score",
    "WSCC_calibrated_confidence_level",

    "manual_review_priority",
]


def clean_value(x) -> str:
    if x is None:
        return ""
    x = str(x).strip()
    if x in BAD_VALUES:
        return ""
    return x


def normalize_text(x) -> str:
    x = clean_value(x).lower()
    x = re.sub(r"\[[^\]]+\]", " ", x)
    x = re.sub(r"[^a-z0-9]+", " ", x)
    x = re.sub(r"\s+", " ", x).strip()
    return x


def compact_text(x) -> str:
    return re.sub(r"[^a-z0-9]+", "", normalize_text(x))


def safe_float(x, default=None):
    x = clean_value(x)
    if not x:
        return default
    try:
        return float(x)
    except Exception:
        return default


def is_uninformative(x) -> bool:
    text = normalize_text(x)
    if not text:
        return True
    return any(re.search(pat, text) for pat in UNINFORMATIVE_PATTERNS)


def tokens(x):
    text = normalize_text(x)
    out = set()
    for t in text.split():
        if len(t) < 3:
            continue
        if t in STOPWORDS:
            continue
        if re.fullmatch(r"\d+", t):
            continue
        out.add(t)
    return out


def row_text(row, cols):
    vals = []
    for c in cols:
        v = clean_value(row.get(c, ""))
        if v:
            vals.append(v)
    return " | ".join(vals)


def source_product(row):
    return clean_value(row.get("source_product", ""))


def annotation_level(row):
    return clean_value(row.get("annotation_evidence_level", ""))


def evidence_text(row):
    return row_text(row, [
        "Preferred_name",
        "Description",
        "diamond_best_hit",
        "annotation_text",
        "PFAMs",
        "pfam_domains_hmmscan",
        "function_class_primary",
        "function_class_all",
        "domain_class_primary",
        "domain_class_all",
    ])


def alias_tokens(text):
    """
    Conservative synonym / abbreviation expansion for biological names.
    This improves matching between full protein names, gene symbols and domain labels.
    """
    n = normalize_text(text)
    c = compact_text(text)
    out = set(tokens(text))
    upper = str(text).upper()

    alias_rules = [
        ("semaphorin", "SEMA"),
        ("sema", "SEMA"),
        ("tgf beta receptor", "TGFBR"),
        ("transforming growth factor beta receptor", "TGFBR"),
        ("tgfbr", "TGFBR"),
        ("tgfr", "TGFBR"),
        ("potassium voltage gated channel", "KCN"),
        ("shaker", "KCN"),
        ("kcn", "KCN"),
        ("gamma aminobutyric acid receptor", "GABA_RECEPTOR"),
        ("gaba receptor", "GABA_RECEPTOR"),
        ("neurotransmitter gated ion channel", "GABA_RECEPTOR"),
        ("kat8 regulatory nsl", "KANSL"),
        ("kansl", "KANSL"),
        ("polycystic kidney disease protein", "PKD"),
        ("pkd", "PKD"),
        ("tripartite motif", "TRIM"),
        ("trim", "TRIM"),
        ("krueppel", "ZINC_FINGER"),
        ("zinc finger", "ZINC_FINGER"),
        ("low density lipoprotein receptor", "LRP"),
        ("ldl receptor", "LRP"),
        ("scavenger receptor", "SCAVENGER_RECEPTOR"),
        ("deleted in malignant brain tumors", "DMBT"),
        ("dmbt", "DMBT"),
        ("plasminogen", "KRINGLE"),
        ("myosin", "MYOSIN"),
        ("collagen", "COLLAGEN"),
        ("cadherin", "CADHERIN"),
        ("mucin", "MUCIN"),
        ("kinase", "KINASE"),
        ("homeobox", "HOMEOBOX"),
        ("wnt", "WNT"),
        ("notch", "NOTCH"),
    ]

    for phrase, alias in alias_rules:
        if phrase.replace(" ", "") in c or phrase in n:
            out.add(alias.lower())

    if re.search(r"\bSEMA\d*[A-Z]?\b", upper):
        out.add("sema")
    if re.search(r"\bTGFBR\d+\b|\bTGFR\d+\b", upper):
        out.add("tgfbr")
    if re.search(r"\bKCN[A-Z0-9]+\b|\bSHK[-_]?\d*\b", upper):
        out.add("kcn")
    if re.search(r"\bKANSL\d+\b", upper):
        out.add("kansl")
    if re.search(r"\bPKD\d*[A-Z0-9]*\b", upper):
        out.add("pkd")
    if re.search(r"\bTRIM\d+\b", upper):
        out.add("trim")
    if re.search(r"\bKLF\d*\b", upper):
        out.add("zinc_finger")

    return out


def jaccard(a, b) -> float:
    a = set(a)
    b = set(b)
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def domain_relation_bonus(a_text, b_text) -> float:
    """
    Detect name-domain compatibility.
    Returns a small bonus when one evidence node describes a family name
    and the other describes its expected domain architecture.
    """
    a = normalize_text(a_text)
    b = normalize_text(b_text)
    joined = a + " | " + b

    rules = [
        (r"semaphorin|sema", r"\bsema\b|thrombospondin|tsp|plexin|psi"),
        (r"tgf|tgfbr|tgfr|transforming growth factor", r"pkinase|pkinase_tyr|ectbetar2|protein kinase"),
        (r"potassium|shaker|voltage gated|kcn", r"ion_trans|ion transport|btb"),
        (r"gamma aminobutyric|gaba|neurotransmitter gated", r"neur_chan|ion channel|ion_trans"),
        (r"tripartite motif|trim", r"ring|b box|b-box|nhl|zf"),
        (r"krueppel|zinc finger", r"zinc finger|c2h2|zf"),
        (r"deleted in malignant brain tumors|dmbt|scavenger receptor", r"srcr|scavenger receptor|lectin|tsp|c-type"),
        (r"polycystic kidney disease|pkd", r"pkd|gps|gain|egf|7tm|fn3"),
        (r"low density lipoprotein receptor|ldl|lrp", r"ldl_recept|egf_ca|egf"),
        (r"plasminogen", r"kringle|sushi"),
        (r"myosin", r"myosin"),
        (r"collagen", r"collagen"),
        (r"cadherin", r"cadherin"),
        (r"mucin", r"mucin"),
        (r"kinase", r"pkinase|protein kinase"),
    ]

    for name_pat, dom_pat in rules:
        if re.search(name_pat, joined) and re.search(dom_pat, joined):
            return 0.25

    return 0.0


def evidence_nodes(row):
    nodes = {}

    src = source_product(row)
    if src:
        nodes["source"] = src

    diamond = row_text(row, ["diamond_best_hit"])
    if diamond:
        nodes["diamond"] = diamond

    eggnog = row_text(row, ["Preferred_name", "Description"])
    if eggnog:
        nodes["eggnog"] = eggnog

    domain = row_text(row, ["PFAMs", "pfam_domains_hmmscan"])
    if domain:
        nodes["domain"] = domain

    function = row_text(row, ["GOs", "KEGG_ko", "KEGG_Pathway", "COG_category"])
    if function:
        nodes["function"] = function

    return nodes


def evidence_graph_consistency(row):
    """
    EGCS: Evidence Graph Consistency Score.

    Each evidence source is treated as a node. Edges are weighted by text,
    synonym and domain-compatibility similarity. The final score reflects
    how coherently all evidence sources support the same annotation.
    """
    nodes = evidence_nodes(row)
    node_names = list(nodes)

    if len(node_names) < 2:
        return {
            "score": 0.0,
            "density": 0.0,
            "supporting_edges": "",
            "conflicting_edges": "",
            "status": "insufficient_evidence",
        }

    edge_scores = []
    supporting = []
    conflicting = []

    for i in range(len(node_names)):
        for j in range(i + 1, len(node_names)):
            a_name = node_names[i]
            b_name = node_names[j]
            a_text = nodes[a_name]
            b_text = nodes[b_name]

            sim1 = jaccard(tokens(a_text), tokens(b_text))
            sim2 = jaccard(alias_tokens(a_text), alias_tokens(b_text))
            bonus = domain_relation_bonus(a_text, b_text)

            sim = max(sim1, sim2) + bonus
            sim = min(1.0, sim)

            # Edges involving source/diamond/eggnog are slightly more informative.
            weight = 1.0
            pair = {a_name, b_name}
            if pair <= {"diamond", "eggnog"}:
                weight = 1.3
            elif "source" in pair:
                weight = 1.2
            elif "domain" in pair:
                weight = 1.1

            edge_scores.append((sim, weight, a_name, b_name))

            if sim >= 0.25:
                supporting.append(f"{a_name}-{b_name}:{sim:.2f}")
            elif sim < 0.05 and a_name in {"source", "diamond", "eggnog"} and b_name in {"source", "diamond", "eggnog"}:
                conflicting.append(f"{a_name}-{b_name}:{sim:.2f}")

    weighted_sum = sum(sim * w for sim, w, _, _ in edge_scores)
    weight_sum = sum(w for _, w, _, _ in edge_scores)
    score = 100 * weighted_sum / weight_sum if weight_sum else 0

    density = sum(1 for sim, _, _, _ in edge_scores if sim >= 0.25) / len(edge_scores)

    if conflicting and score < 30:
        status = "conflicting_evidence"
    elif score >= 70:
        status = "high_graph_agreement"
    elif score >= 45:
        status = "moderate_graph_agreement"
    elif score >= 25:
        status = "weak_graph_agreement"
    else:
        status = "low_graph_agreement"

    return {
        "score": round(score, 1),
        "density": round(density, 3),
        "supporting_edges": ";".join(supporting),
        "conflicting_edges": ";".join(conflicting),
        "status": status,
    }


def classify_domain_architecture(row):
    """
    DACS: Domain Architecture Consistency Status.
    """
    src = source_product(row)
    ev = evidence_text(row)
    all_text = normalize_text(src + " | " + ev)
    pfam_text = normalize_text(row_text(row, ["PFAMs", "pfam_domains_hmmscan"]))

    if not pfam_text:
        return "no_domain"

    rules = [
        (r"semaphorin|sema", r"\bsema\b|thrombospondin|tsp|plexin|psi", "canonical_domain_architecture"),
        (r"tgf|tgfbr|tgfr|transforming growth factor", r"pkinase|pkinase_tyr|ectbetar2|protein kinase", "canonical_domain_architecture"),
        (r"potassium|shaker|voltage gated|kcn", r"ion_trans|ion transport|btb", "canonical_domain_architecture"),
        (r"gamma aminobutyric|gaba|neurotransmitter gated", r"neur_chan|neurotransmitter gated|ion channel|ion_trans", "canonical_domain_architecture"),
        (r"tripartite motif|trim", r"ring|b box|b-box|nhl|zf", "family_level_domain_support"),
        (r"krueppel|zinc finger", r"zinc finger|c2h2|zf", "canonical_domain_architecture"),
        (r"deleted in malignant brain tumors|dmbt|scavenger receptor", r"srcr|scavenger receptor|lectin|tsp|c-type", "family_level_domain_support"),
        (r"polycystic kidney disease|pkd", r"pkd|gps|gain|egf|7tm|fn3", "canonical_domain_architecture"),
        (r"low density lipoprotein receptor|ldl|lrp", r"ldl_recept|egf_ca|egf", "family_level_domain_support"),
        (r"plasminogen", r"kringle|sushi", "family_level_domain_support"),
        (r"myosin", r"myosin", "canonical_domain_architecture"),
        (r"collagen", r"collagen", "canonical_domain_architecture"),
        (r"cadherin", r"cadherin", "canonical_domain_architecture"),
        (r"mucin", r"mucin", "canonical_domain_architecture"),
        (r"kinase", r"pkinase|protein kinase", "domain_support_only"),
    ]

    for name_pat, dom_pat, status in rules:
        if re.search(name_pat, all_text) and re.search(dom_pat, pfam_text):
            return status

    return "domain_support_only"


def extract_family_calls(text):
    """
    Extract coarse family + paralog identifiers from free text.

    Examples:
    TRIM2_BOVIN -> ("TRIM", "2")
    tripartite motif-containing protein 3-like -> ("TRIM", "3")
    TGF-beta receptor type-2-like -> ("TGFBR", "2")
    SEMA5A -> ("SEMA", "5A")
    """
    raw = str(text)
    upper = raw.upper()
    norm = normalize_text(raw)
    calls = set()

    regex_rules = [
        (r"\bTRIM[-_ ]?(\d+[A-Z]?)\b", "TRIM"),
        (r"tripartite motif.*?(\d+[a-z]?)", "TRIM"),
        (r"\bTGFBR[-_ ]?(\d+[A-Z]?)\b", "TGFBR"),
        (r"\bTGFR[-_ ]?(\d+[A-Z]?)\b", "TGFBR"),
        (r"tgf beta receptor type[-_ ]?(\d+[a-z]?)", "TGFBR"),
        (r"\bSEMA[-_ ]?(\d+[A-Z]?)\b", "SEMA"),
        (r"semaphorin[-_ ]?(\d+[a-z]?)", "SEMA"),
        (r"\bKCN[A-Z]*[-_ ]?(\d+[A-Z]?)\b", "KCN"),
        (r"shaker", "KCN"),
        (r"\bKANSL[-_ ]?(\d+[A-Z]?)\b", "KANSL"),
        (r"kat8 regulatory nsl complex subunit[-_ ]?(\d+[a-z]?)", "KANSL"),
        (r"\bPKD[-_ ]?(\d+[A-Z0-9]*)\b", "PKD"),
        (r"polycystic kidney disease protein[-_ ]?(\d+[a-z0-9]*)", "PKD"),
        (r"\bLRP[-_ ]?(\d+[A-Z]?)\b", "LRP"),
        (r"low density lipoprotein receptor related protein[-_ ]?(\d+[a-z]?)", "LRP"),
        (r"\bKLF[-_ ]?(\d+[A-Z]?)\b", "KLF"),
    ]

    for pat, family in regex_rules:
        m = re.search(pat, upper)
        if not m:
            m = re.search(pat, norm)
        if m:
            try:
                num = m.group(1).upper()
            except Exception:
                num = ""
            calls.add((family, num))

    return calls


def paralog_ambiguity(row, egcs):
    """
    PAI: Paralog Ambiguity Index.

    PAI estimates whether the evidence supports an exact gene-level name
    or only a family-level annotation.
    """
    src = source_product(row)
    ev = evidence_text(row)
    src_calls = extract_family_calls(src)
    ev_calls = extract_family_calls(ev)

    detected = sorted(src_calls | ev_calls)

    if not detected:
        # No gene-family-style paralog numbers were detected.
        if egcs["score"] >= 70:
            return 0.15, "clear_or_not_applicable", "", "No obvious paralog ambiguity and high evidence agreement."
        if egcs["score"] >= 45:
            return 0.30, "likely_gene_or_family_level", "", "Moderate evidence agreement without explicit paralog conflict."
        return 0.50, "unresolved", "", "Insufficient agreement to resolve annotation specificity."

    # Same family but different numbered paralogs.
    same_family_conflict = False
    for sf, sn in src_calls:
        for ef, en in ev_calls:
            if sf == ef and sn and en and sn != en:
                same_family_conflict = True

    if same_family_conflict:
        return 0.75, "family_level_only", ";".join(f"{f}:{n}" for f, n in detected), "Different paralog numbers detected within the same family."

    if src_calls and ev_calls and src_calls & ev_calls:
        return 0.10, "clear_gene_level", ";".join(f"{f}:{n}" for f, n in detected), "Source and evidence support the same family/paralog label."

    if src_calls and not ev_calls:
        return 0.60, "source_specific_evidence_family_level", ";".join(f"{f}:{n}" for f, n in detected), "Source product is specific, but supporting evidence does not resolve the paralog."

    if ev_calls and not src_calls:
        return 0.40, "evidence_specific_source_unspecific", ";".join(f"{f}:{n}" for f, n in detected), "Evidence suggests a family label, while source product is less specific."

    return 0.50, "unresolved", ";".join(f"{f}:{n}" for f, n in detected), "Paralog resolution remains uncertain."


def classify_annotation_consistency(row, dacs_status, egcs, pai_index):
    """
    ACAD: Annotation Conflict and Ambiguity Detector.
    """
    src = source_product(row)
    ev = evidence_text(row)
    lev = annotation_level(row)

    src_informative = not is_uninformative(src)
    ev_present = bool(normalize_text(ev)) and lev != "no_annotation"

    if not src_informative and not ev_present:
        return "no_external_support"

    if src_informative and not ev_present:
        return "source_annotation_more_specific_than_evidence"

    if not src_informative and ev_present:
        return "evidence_more_specific_than_source"

    if pai_index >= 0.70:
        return "paralog_ambiguous"

    if egcs["status"] == "conflicting_evidence":
        return "potential_conflict"

    if egcs["score"] >= 70:
        return "consistent"

    if egcs["score"] >= 45:
        if dacs_status == "canonical_domain_architecture":
            return "domain_level_consistent"
        if dacs_status == "family_level_domain_support":
            return "family_level_consistent"
        return "synonym_or_symbol_difference"

    if dacs_status == "canonical_domain_architecture":
        return "domain_level_consistent"

    if dacs_status == "family_level_domain_support":
        return "family_level_consistent"

    if lev == "diamond_eggnog_pfam_supported" and egcs["score"] < 25:
        return "potential_conflict"

    return "unclear"


def compute_raw_pacs(row, egcs, pai_index, dacs_status, acad_status):
    """
    PACS: Protein Annotation Confidence Score.
    Raw interpretable score before weak-supervised calibration.
    """
    score = 0.0
    lev = annotation_level(row)

    evidence_base = {
        "diamond_eggnog_pfam_supported": 50,
        "eggnog_and_pfam_supported": 36,
        "diamond_eggnog_supported": 34,
        "diamond_pfam_supported": 34,
        "eggnog_supported_only": 20,
        "diamond_supported_only": 18,
        "pfam_supported_only": 18,
        "single_or_weak_evidence": 10,
        "no_annotation": 0,
    }
    score += evidence_base.get(lev, 6)

    pident = safe_float(row.get("diamond_pident"))
    qcov = safe_float(row.get("diamond_qcovhsp"))
    scov = safe_float(row.get("diamond_scovhsp"))
    evalue = safe_float(row.get("diamond_evalue"))

    if clean_value(row.get("diamond_best_hit", "")):
        score += 5

    if pident is not None:
        if pident >= 70:
            score += 7
        elif pident >= 50:
            score += 5
        elif pident >= 30:
            score += 3

    covs = [x for x in [qcov, scov] if x is not None]
    if covs:
        avg_cov = sum(covs) / len(covs)
        if avg_cov >= 70:
            score += 6
        elif avg_cov >= 40:
            score += 3

    if evalue is not None:
        if evalue <= 1e-30:
            score += 4
        elif evalue <= 1e-10:
            score += 2
        elif evalue <= 1e-5:
            score += 1

    if clean_value(row.get("Preferred_name", "")):
        score += 4
    if clean_value(row.get("Description", "")):
        score += 4
    if clean_value(row.get("GOs", "")):
        score += 2
    if clean_value(row.get("KEGG_ko", "")) or clean_value(row.get("KEGG_Pathway", "")):
        score += 2
    if clean_value(row.get("COG_category", "")):
        score += 1

    if clean_value(row.get("PFAMs", "")) or clean_value(row.get("pfam_domains_hmmscan", "")):
        score += 5

    n_pfam = safe_float(row.get("n_pfam_domains"), 0)
    if n_pfam:
        score += min(4, n_pfam)

    # EGCS contribution
    score += 0.18 * float(egcs["score"])

    # DACS contribution
    if dacs_status == "canonical_domain_architecture":
        score += 5
    elif dacs_status == "family_level_domain_support":
        score += 3
    elif dacs_status == "domain_support_only":
        score += 1

    # ACAD contribution
    acad_bonus = {
        "consistent": 8,
        "synonym_or_symbol_difference": 7,
        "domain_level_consistent": 5,
        "family_level_consistent": 3,
        "evidence_more_specific_than_source": 3,
        "paralog_ambiguous": -2,
        "unclear": -3,
        "source_annotation_more_specific_than_evidence": -8,
        "potential_conflict": -15,
        "no_external_support": -15,
    }
    score += acad_bonus.get(acad_status, 0)

    # PAI penalty
    score -= 18 * float(pai_index)

    src = source_product(row)
    if re.search(r"\bpartial\b", normalize_text(src)):
        score -= 3
    if re.search(r"low quality", normalize_text(src)):
        score -= 8

    score = max(0, min(100, round(score, 1)))
    return score


def confidence_level(score):
    score = float(score)
    if score >= 85:
        return "very_high_confidence"
    if score >= 70:
        return "high_confidence"
    if score >= 50:
        return "medium_confidence"
    if score >= 30:
        return "low_confidence"
    return "weak_or_no_annotation"


def first_pfam_name(row):
    text = clean_value(row.get("PFAMs", ""))
    if text:
        return text.split(",")[0].split(";")[0].strip()
    dom = clean_value(row.get("pfam_domains_hmmscan", ""))
    if dom:
        return dom.split("|")[0].split(";")[0].strip()
    return ""


def classify_recommended_name(row, score, acad_status):
    """
    CFNR: Conservative Functional Name Recommendation.
    """
    src = source_product(row)
    pref = clean_value(row.get("Preferred_name", ""))
    desc = clean_value(row.get("Description", ""))
    pfam = first_pfam_name(row)
    src_informative = not is_uninformative(src)

    if acad_status == "potential_conflict":
        return {
            "CFNR_recommended_name": src if src_informative else (pref or "uncharacterized protein"),
            "CFNR_recommended_name_level": "conflict_candidate",
            "CFNR_recommended_name_reason": "Source annotation and multi-evidence annotation may be inconsistent; manual review recommended.",
        }

    if src_informative and float(score) >= 50 and acad_status in {
        "consistent",
        "synonym_or_symbol_difference",
        "domain_level_consistent",
    }:
        return {
            "CFNR_recommended_name": src,
            "CFNR_recommended_name_level": "gene_like",
            "CFNR_recommended_name_reason": "Source product is supported by sequence, orthology and/or domain evidence.",
        }

    if src_informative and acad_status in {"family_level_consistent", "paralog_ambiguous"}:
        conservative = re.sub(r"\b\d+[a-z]?\b", "", src, flags=re.I)
        conservative = re.sub(r"\s+", " ", conservative).strip()
        if not conservative:
            conservative = src
        return {
            "CFNR_recommended_name": conservative,
            "CFNR_recommended_name_level": "family_like",
            "CFNR_recommended_name_reason": "Evidence supports the protein family, but exact paralog identity is uncertain.",
        }

    if pref and not is_uninformative(pref):
        return {
            "CFNR_recommended_name": pref + "-like",
            "CFNR_recommended_name_level": "orthology_like",
            "CFNR_recommended_name_reason": "eggNOG preferred name provides the most informative conservative label.",
        }

    if desc and not is_uninformative(desc):
        return {
            "CFNR_recommended_name": desc,
            "CFNR_recommended_name_level": "functional_class",
            "CFNR_recommended_name_reason": "Functional description is available but exact gene-level identity is uncertain.",
        }

    if pfam:
        return {
            "CFNR_recommended_name": pfam + " domain-containing protein",
            "CFNR_recommended_name_level": "domain_containing",
            "CFNR_recommended_name_reason": "Only domain-level evidence is available.",
        }

    return {
        "CFNR_recommended_name": "uncharacterized protein",
        "CFNR_recommended_name_level": "uncharacterized",
        "CFNR_recommended_name_reason": "No reliable functional evidence was detected.",
    }


def pseudo_label_for_calibration(row):
    """
    Weak-supervised calibration label.

    Positive labels are inferred from informative source products that agree
    with multi-evidence annotation. Negative labels are inferred from informative
    source products with no support or possible conflict.

    This is not a gold-standard label; it is a weak label for score calibration.
    """
    src = source_product(row)
    if is_uninformative(src):
        return "NA"

    acad = row.get("ACAD_annotation_consistency_status", "")
    egcs = safe_float(row.get("EGCS_evidence_graph_score"), 0)
    lev = annotation_level(row)

    positive_states = {
        "consistent",
        "synonym_or_symbol_difference",
        "domain_level_consistent",
        "family_level_consistent",
    }

    negative_states = {
        "potential_conflict",
        "source_annotation_more_specific_than_evidence",
        "no_external_support",
    }

    if acad in positive_states or egcs >= 45:
        return "positive"

    if acad in negative_states or lev == "no_annotation":
        return "negative"

    return "NA"


def manual_review_priority(row):
    acad = row.get("ACAD_annotation_consistency_status", "")
    pai = safe_float(row.get("PAI_paralog_ambiguity_index"), 0.5)
    calibrated = safe_float(row.get("WSCC_calibrated_confidence_score"), None)
    raw = safe_float(row.get("PACS_annotation_confidence_score"), 0)

    if acad in {"potential_conflict", "source_annotation_more_specific_than_evidence"}:
        return "high"

    if pai is not None and pai >= 0.70:
        return "high"

    if acad in {"paralog_ambiguous", "family_level_consistent", "unclear"}:
        return "medium"

    score_to_use = calibrated if calibrated is not None else raw
    if score_to_use < 30:
        return "medium"

    return "low"


def score_row_before_calibration(row):
    dacs = classify_domain_architecture(row)
    egcs = evidence_graph_consistency(row)
    pai_index, pai_level, pai_families, pai_reason = paralog_ambiguity(row, egcs)
    acad = classify_annotation_consistency(row, dacs, egcs, pai_index)
    pacs = compute_raw_pacs(row, egcs, pai_index, dacs, acad)
    pacs_level = confidence_level(pacs)
    rec = classify_recommended_name(row, pacs, acad)

    out = dict(row)
    out["PACS_annotation_confidence_score"] = str(pacs)
    out["PACS_annotation_confidence_level"] = pacs_level

    out["EGCS_evidence_graph_score"] = str(egcs["score"])
    out["EGCS_evidence_graph_density"] = str(egcs["density"])
    out["EGCS_supporting_edges"] = egcs["supporting_edges"]
    out["EGCS_conflicting_edges"] = egcs["conflicting_edges"]
    out["EGCS_graph_status"] = egcs["status"]

    out["PAI_paralog_ambiguity_index"] = str(round(pai_index, 3))
    out["PAI_paralog_resolution_level"] = pai_level
    out["PAI_detected_families"] = pai_families
    out["PAI_reason"] = pai_reason

    out["DACS_domain_architecture_status"] = dacs
    out["ACAD_annotation_consistency_status"] = acad
    out.update(rec)

    out["WSCC_pseudo_label"] = "NA"
    out["WSCC_calibration_bin"] = "NA"
    out["WSCC_bin_positive_rate"] = "NA"
    out["WSCC_calibrated_confidence_score"] = "NA"
    out["WSCC_calibrated_confidence_level"] = "NA"

    out["manual_review_priority"] = "NA"
    return out


def apply_weak_supervised_calibration(rows):
    """
    WSCC: Weak-Supervised Confidence Calibration.

    Uses informative source products as weak labels to empirically calibrate
    the raw PACS score. No external ML dependency is required.
    """
    bins = defaultdict(lambda: Counter())

    for row in rows:
        label = pseudo_label_for_calibration(row)
        row["WSCC_pseudo_label"] = label

        raw = safe_float(row.get("PACS_annotation_confidence_score"), 0)
        bin_id = int(min(9, max(0, math.floor(raw / 10))))
        row["WSCC_calibration_bin"] = str(bin_id)

        if label in {"positive", "negative"}:
            bins[bin_id][label] += 1

    global_counts = Counter()
    for c in bins.values():
        global_counts.update(c)

    enough_labels = (global_counts["positive"] + global_counts["negative"]) >= 50

    for row in rows:
        raw = safe_float(row.get("PACS_annotation_confidence_score"), 0)
        bin_id = int(row.get("WSCC_calibration_bin", "0"))

        c = bins[bin_id]
        # Laplace smoothing.
        if enough_labels:
            pos_rate = (c["positive"] + 1) / (c["positive"] + c["negative"] + 2)
            calibrated = 0.65 * raw + 0.35 * (pos_rate * 100)
            row["WSCC_bin_positive_rate"] = f"{pos_rate:.3f}"
            row["WSCC_calibrated_confidence_score"] = f"{round(calibrated, 1)}"
            row["WSCC_calibrated_confidence_level"] = confidence_level(calibrated)
        else:
            row["WSCC_bin_positive_rate"] = "NA"
            row["WSCC_calibrated_confidence_score"] = f"{round(raw, 1)}"
            row["WSCC_calibrated_confidence_level"] = confidence_level(raw)

    for row in rows:
        row["manual_review_priority"] = manual_review_priority(row)

    return rows


def enrich_annotation_table(in_tsv, out_tsv):
    in_tsv = Path(in_tsv)
    out_tsv = Path(out_tsv)
    tmp_out = out_tsv.with_suffix(out_tsv.suffix + ".tmp")

    with in_tsv.open("r", encoding="utf-8", errors="ignore", newline="") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        original_fields = reader.fieldnames or []
        rows = [score_row_before_calibration(row) for row in reader]

    rows = apply_weak_supervised_calibration(rows)

    fieldnames = list(original_fields)
    for f in NEW_FIELDS:
        if f not in fieldnames:
            fieldnames.append(f)

    with tmp_out.open("w", encoding="utf-8", newline="") as fout:
        writer = csv.DictWriter(
            fout,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

    tmp_out.replace(out_tsv)
    return len(rows)
