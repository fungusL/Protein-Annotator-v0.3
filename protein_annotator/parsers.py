from pathlib import Path
import csv


def parse_fasta_ids(fasta):
    ids = []
    with open(fasta, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def to_float(x):
    try:
        return float(x)
    except Exception:
        return None


def parse_diamond(diamond_raw, outfile):
    best = {}

    with open(diamond_raw, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip():
                continue

            arr = line.rstrip("\n").split("\t")
            if len(arr) < 8:
                continue

            qseqid, sseqid, pident, length, evalue, bitscore, qcovhsp, scovhsp = arr[:8]
            score = to_float(bitscore) or 0

            if qseqid not in best or score > (to_float(best[qseqid]["diamond_bitscore"]) or 0):
                best[qseqid] = {
                    "gene_id": qseqid,
                    "diamond_best_hit": sseqid,
                    "diamond_pident": pident,
                    "diamond_align_len": length,
                    "diamond_evalue": evalue,
                    "diamond_bitscore": bitscore,
                    "diamond_qcovhsp": qcovhsp,
                    "diamond_scovhsp": scovhsp,
                }

    cols = [
        "gene_id",
        "diamond_best_hit",
        "diamond_pident",
        "diamond_align_len",
        "diamond_evalue",
        "diamond_bitscore",
        "diamond_qcovhsp",
        "diamond_scovhsp",
    ]

    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    with open(outfile, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        for gid in sorted(best):
            writer.writerow(best[gid])

    print(f"[OK] Parsed DIAMOND: {len(best)} hits -> {outfile}")
    return outfile


def parse_eggnog(emapper_file, outfile):
    header = None
    rows = []

    with open(emapper_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("#query"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                break

        if header is None:
            raise RuntimeError(f"Cannot find #query header in eggNOG file: {emapper_file}")

        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            arr = line.rstrip("\n").split("\t")
            if len(arr) < len(header):
                arr += [""] * (len(header) - len(arr))

            row = dict(zip(header, arr))
            row["gene_id"] = row.get("query", "")
            rows.append(row)

    wanted = [
        "gene_id",
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
        "KEGG_Reaction",
        "BRITE",
        "KEGG_TC",
        "CAZy",
        "BiGG_Reaction",
        "PFAMs",
    ]

    existing = [c for c in wanted if rows and c in rows[0]]

    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    with open(outfile, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=existing, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"[OK] Parsed eggNOG: {len(rows)} genes -> {outfile}")
    return outfile


def parse_hmmscan(domtblout, outfile, ievalue_cutoff=1e-5):
    domains = {}

    with open(domtblout, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            arr = line.rstrip("\n").split(maxsplit=22)
            if len(arr) < 22:
                continue

            target_name = arr[0]
            target_acc = arr[1]
            query_name = arr[3]
            i_evalue = arr[12]
            dom_score = arr[13]
            desc = arr[22] if len(arr) > 22 else ""

            ie = to_float(i_evalue)
            if ie is None or ie > ievalue_cutoff:
                continue

            item = f"{target_name}|{target_acc}|iE={i_evalue}|score={dom_score}|{desc}"
            domains.setdefault(query_name, []).append(item)

    rows = []
    for gid, doms in sorted(domains.items()):
        rows.append({
            "gene_id": gid,
            "n_pfam_domains": str(len(doms)),
            "pfam_domains_hmmscan": ";".join(doms),
        })

    cols = ["gene_id", "n_pfam_domains", "pfam_domains_hmmscan"]

    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    with open(outfile, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] Parsed Pfam/hmmscan: {len(rows)} genes -> {outfile}")
    return outfile
