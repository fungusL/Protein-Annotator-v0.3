from pathlib import Path

from .checks import doctor
from .runners import run_diamond, run_eggnog, run_hmmscan
from .parsers import parse_diamond, parse_eggnog, parse_hmmscan
from .merge import merge_annotations
from .tables import filter_annotation_tables


def annotate_proteome(
    protein_fasta,
    outdir,
    prefix,
    diamond_db,
    pfam_hmm,
    eggnog_data_dir=None,
    threads=8,
    force=False,
    skip_doctor=False,
):
    protein_fasta = Path(protein_fasta)
    outdir = Path(outdir)

    if not protein_fasta.exists():
        raise FileNotFoundError(f"Protein fasta not found: {protein_fasta}")

    if not skip_doctor:
        doctor(
            diamond_db=diamond_db,
            pfam_hmm=pfam_hmm,
            eggnog_data_dir=eggnog_data_dir,
        )

    diamond_dir = outdir / "01_diamond"
    eggnog_dir = outdir / "02_eggnog"
    pfam_dir = outdir / "03_pfam"
    merged_dir = outdir / "04_merged"
    highconf_dir = outdir / "05_high_confidence"

    for d in [diamond_dir, eggnog_dir, pfam_dir, merged_dir, highconf_dir]:
        d.mkdir(parents=True, exist_ok=True)

    print("\n===== Step 1. DIAMOND best-hit =====")
    diamond_raw = diamond_dir / f"{prefix}.diamond.raw.tsv"
    diamond_best = diamond_dir / f"{prefix}.diamond.besthit.tsv"

    run_diamond(protein_fasta, diamond_db, diamond_raw, threads=threads, force=force)
    parse_diamond(diamond_raw, diamond_best)

    print("\n===== Step 2. eggNOG-mapper =====")
    eggnog_raw = run_eggnog(
        protein_fasta=protein_fasta,
        outdir=eggnog_dir,
        prefix=prefix,
        eggnog_data_dir=eggnog_data_dir,
        threads=threads,
        force=force,
    )

    eggnog_parsed = eggnog_dir / f"{prefix}.eggnog.parsed.tsv"
    parse_eggnog(eggnog_raw, eggnog_parsed)

    print("\n===== Step 3. hmmscan / Pfam =====")
    hmmscan_raw = pfam_dir / f"{prefix}.hmmscan.domtblout"
    pfam_parsed = pfam_dir / f"{prefix}.pfam.parsed.tsv"

    run_hmmscan(protein_fasta, pfam_hmm, hmmscan_raw, threads=threads, force=force)
    parse_hmmscan(hmmscan_raw, pfam_parsed)

    print("\n===== Step 4. Merge annotations =====")
    merged_file = merged_dir / f"{prefix}.final_merged_annotation.tsv"

    merge_annotations(
        protein_fasta=protein_fasta,
        diamond_file=diamond_best,
        eggnog_file=eggnog_parsed,
        pfam_file=pfam_parsed,
        outfile=merged_file,
    )

    print("\n===== Step 5. High-confidence filtering =====")
    result_files = filter_annotation_tables(
        annotation_tsv=merged_file,
        outdir=highconf_dir,
        prefix=prefix,
    )

    print("\n===== Done =====")
    print(f"Final merged annotation: {merged_file}")
    print(f"High-confidence tables: {highconf_dir}")

    return {
        "merged_annotation": str(merged_file),
        "high_confidence_dir": str(highconf_dir),
        **result_files,
    }
