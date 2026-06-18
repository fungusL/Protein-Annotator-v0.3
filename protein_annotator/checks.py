from pathlib import Path
import shutil
import subprocess


TOOLS = {
    "diamond": "diamond",
    "emapper.py": "emapper.py",
    "hmmscan": "hmmscan",
    "hmmpress": "hmmpress",
}

CONDA_PACKAGES = {
    "diamond": "diamond",
    "emapper.py": "eggnog-mapper",
    "hmmscan": "hmmer",
    "hmmpress": "hmmer",
}


def which(cmd):
    return shutil.which(cmd)


def _capture(cmd):
    try:
        return subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True).strip()
    except Exception as e:
        return str(e)


def doctor(diamond_db=None, pfam_hmm=None, eggnog_data_dir=None):
    print("===== Protein Annotator: doctor =====")
    ok = True

    print("\n[1] Tools")
    for name, cmd in TOOLS.items():
        path = which(cmd)
        if path:
            print(f"[OK] {name}: {path}")
        else:
            print(f"[MISSING] {name}")
            ok = False

    print("\n[2] Versions")
    if which("diamond"):
        print("[diamond]", _capture(["diamond", "version"]))
    if which("emapper.py"):
        print("[emapper.py]", _capture(["emapper.py", "--version"]))
    if which("hmmscan"):
        print("[hmmscan]", _capture(["hmmscan", "-h"]).splitlines()[0])

    print("\n[3] Databases")

    if diamond_db:
        p = Path(diamond_db)
        if p.exists():
            print(f"[OK] DIAMOND db: {p}")
        else:
            print(f"[MISSING] DIAMOND db: {p}")
            ok = False
    else:
        print("[SKIP] DIAMOND db not provided")

    if pfam_hmm:
        p = Path(pfam_hmm)
        if p.exists():
            print(f"[OK] Pfam HMM: {p}")
            idx = [Path(str(p) + s) for s in [".h3f", ".h3i", ".h3m", ".h3p"]]
            miss = [x for x in idx if not x.exists()]
            if miss:
                print("[WARNING] Pfam hmmpress index missing")
                print(f"Run: protein-annotate prepare-pfam --pfam-hmm {p}")
            else:
                print("[OK] Pfam hmmpress index found")
        else:
            print(f"[MISSING] Pfam HMM: {p}")
            ok = False
    else:
        print("[SKIP] Pfam HMM not provided")

    if eggnog_data_dir:
        d = Path(eggnog_data_dir)
        if d.exists():
            has_db = any(x.suffix in [".db", ".dmnd"] for x in d.glob("*"))
            if has_db:
                print(f"[OK] eggNOG data dir: {d}")
            else:
                print(f"[WARNING] eggNOG dir exists but no .db/.dmnd found: {d}")
        else:
            print(f"[MISSING] eggNOG data dir: {d}")
            ok = False
    else:
        print("[SKIP] eggNOG data dir not provided")

    print("\n===== Summary =====")
    print("[OK] Environment looks usable." if ok else "[WARNING] Some tools/files are missing.")
    return ok


def install_missing_tools(manager="mamba", dry_run=False):
    missing = [name for name, cmd in TOOLS.items() if which(cmd) is None]

    if not missing:
        print("[OK] All required tools are already installed.")
        return

    packages = sorted(set(CONDA_PACKAGES[x] for x in missing))
    cmd = [
        manager, "install", "-y",
        "-c", "conda-forge",
        "-c", "bioconda",
    ] + packages

    print("[INFO] Missing tools:", ", ".join(missing))
    print("[INFO] Packages:", " ".join(packages))
    print("[CMD]", " ".join(cmd))

    if dry_run:
        print("[DRY-RUN] Not executed.")
        return

    subprocess.run(cmd, check=True)


def prepare_pfam(pfam_hmm):
    pfam = Path(pfam_hmm)

    if not pfam.exists():
        raise FileNotFoundError(f"Pfam HMM not found: {pfam}")

    index_files = [Path(str(pfam) + s) for s in [".h3f", ".h3i", ".h3m", ".h3p"]]

    if all(x.exists() for x in index_files):
        print("[OK] Pfam hmmpress index already exists.")
        return

    if which("hmmpress") is None:
        raise RuntimeError("hmmpress not found. Run: protein-annotate install-tools")

    cmd = ["hmmpress", str(pfam)]
    print("[CMD]", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print("[OK] Pfam hmmpress index generated.")
