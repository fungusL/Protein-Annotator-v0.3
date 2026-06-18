from pathlib import Path
import subprocess


def run_command(cmd, log_file=None):
    print("[CMD]", " ".join(map(str, cmd)))
    if log_file:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        with open(log_file, "w", encoding="utf-8") as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
    else:
        subprocess.run(cmd, check=True)


def skip_if_exists(outfile, force=False):
    outfile = Path(outfile)
    if outfile.exists() and outfile.stat().st_size > 0 and not force:
        print(f"[SKIP] Existing file found: {outfile}")
        return True
    return False


def run_diamond(protein_fasta, diamond_db, outfile, threads=8, force=False):
    outfile = Path(outfile)
    if skip_if_exists(outfile, force):
        return outfile

    outfile.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "diamond", "blastp",
        "--query", str(protein_fasta),
        "--db", str(diamond_db),
        "--out", str(outfile),
        "--outfmt", "6",
        "qseqid", "sseqid", "pident", "length", "evalue",
        "bitscore", "qcovhsp", "scovhsp",
        "--max-target-seqs", "1",
        "--evalue", "1e-5",
        "--threads", str(threads),
    ]
    run_command(cmd)
    return outfile


def run_eggnog(protein_fasta, outdir, prefix, eggnog_data_dir=None, threads=8, force=False):
    """
    Run eggNOG-mapper.

    Default mode:
    - use DIAMOND
    - disable DIAMOND iterative search
    - use fast sensitivity mode
    """
    from pathlib import Path
    import shutil

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ann_file = outdir / f"{prefix}.emapper.annotations"
    hits_file = outdir / f"{prefix}.emapper.hits"
    temp_dir = outdir / f"{prefix}.emapper_tmp"
    temp_dir.mkdir(parents=True, exist_ok=True)

    if ann_file.exists() and ann_file.stat().st_size > 0 and not force:
        print(f"[SKIP] eggNOG-mapper output exists: {ann_file}", flush=True)
        return ann_file

    if force:
        for f in [ann_file, hits_file]:
            if f.exists():
                print(f"[REMOVE] Existing eggNOG output: {f}", flush=True)
                f.unlink()

    if shutil.which("emapper.py") is None:
        raise RuntimeError("emapper.py not found. Run: protein-annotate install-tools")

    cmd = [
        "emapper.py",
        "-i", str(protein_fasta),
        "--itype", "proteins",
        "-m", "diamond",
        "-o", prefix,
        "--output_dir", str(outdir),
        "--cpu", str(threads),
        "--override",
        "--temp_dir", str(temp_dir),

        # speed settings
        "--dmnd_iterate", "no",
        "--sensmode", "fast",
    ]

    if eggnog_data_dir:
        cmd += ["--data_dir", str(eggnog_data_dir)]

    run_command(cmd)

    if not ann_file.exists():
        raise RuntimeError(f"eggNOG-mapper finished, but annotation file not found: {ann_file}")

    print(f"[OK] eggNOG-mapper annotation ready: {ann_file}", flush=True)
    return ann_file

def run_hmmscan(protein_fasta, pfam_hmm, outfile, threads=8, force=False):
    outfile = Path(outfile)
    if skip_if_exists(outfile, force):
        return outfile

    outfile.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "hmmscan",
        "--cpu", str(threads),
        "--domtblout", str(outfile),
        str(pfam_hmm),
        str(protein_fasta),
    ]

    run_command(cmd, log_file=str(outfile) + ".log")
    return outfile
