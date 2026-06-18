from pathlib import Path
import gzip
import shutil
import subprocess
import tarfile
import urllib.request


UNIPROT_SPROT_URL = (
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
    "knowledgebase/complete/uniprot_sprot.fasta.gz"
)

PFAM_A_HMM_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
)

# eggNOG-mapper v2.1.x expects emapperdb-5.0.2.
# Old host eggnogdb.embl.de may fail; use eggnog5.embl.de directly.
EGGNOG_BASE_URLS = [
    "http://eggnog5.embl.de/download/emapperdb-5.0.2",
]


def run_command(cmd):
    print("[CMD]", " ".join(map(str, cmd)), flush=True)
    subprocess.run(cmd, check=True)


def is_nonempty_file(path):
    path = Path(path)
    return path.exists() and path.is_file() and path.stat().st_size > 0


def clean_zero_size_files(files):
    for f in files:
        f = Path(f)
        if f.exists() and f.is_file() and f.stat().st_size == 0:
            print(f"[WARNING] Remove zero-byte file: {f}", flush=True)
            f.unlink()


def download_file(urls, outfile, force=False, max_attempts=50):
    """
    Robust resumable downloader.

    Important behavior:
    - Download to outfile.part first.
    - Keep outfile.part after network interruption.
    - Next run resumes from outfile.part.
    - Rename to outfile only after wget exits successfully and file is non-empty.
    """
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    if isinstance(urls, str):
        urls = [urls]

    if is_nonempty_file(outfile) and not force:
        print(f"[SKIP] Existing non-empty file found: {outfile}", flush=True)
        return outfile

    if outfile.exists() and force:
        print(f"[REMOVE] Existing target file because --force: {outfile}", flush=True)
        outfile.unlink()

    part = Path(str(outfile) + ".part")

    # 0-byte partial files are useless; non-empty partial files should be kept for resume.
    if part.exists() and part.stat().st_size == 0:
        print(f"[REMOVE] Zero-byte partial download: {part}", flush=True)
        part.unlink()

    last_error = None

    for url in urls:
        for attempt in range(1, max_attempts + 1):
            print(f"[DOWNLOAD] attempt {attempt}/{max_attempts}: {url}", flush=True)
            print(f"[TO] {outfile}", flush=True)

            if part.exists():
                print(f"[RESUME] Existing partial file: {part} ({part.stat().st_size} bytes)", flush=True)

            try:
                if shutil.which("wget"):
                    cmd = [
                        "wget",
                        "--continue",
                        "--tries=3",
                        "--read-timeout=90",
                        "--timeout=60",
                        "--waitretry=10",
                        "--retry-connrefused",
                        "-O", str(part),
                        url,
                    ]
                    run_command(cmd)
                else:
                    # Python fallback cannot reliably resume; wget is strongly preferred.
                    with urllib.request.urlopen(url, timeout=90) as response, open(part, "ab") as out:
                        shutil.copyfileobj(response, out)

                if not is_nonempty_file(part):
                    raise RuntimeError(f"Downloaded file is empty: {part}")

                part.rename(outfile)

                if not is_nonempty_file(outfile):
                    raise RuntimeError(f"Downloaded file is empty after rename: {outfile}")

                print(f"[OK] Downloaded: {outfile}", flush=True)
                return outfile

            except Exception as e:
                last_error = e
                print(f"[WARNING] Download interrupted or failed: {url}", flush=True)
                print(f"[WARNING] {e}", flush=True)

                if part.exists() and part.stat().st_size > 0:
                    print(f"[KEEP] Partial file preserved for resume: {part} ({part.stat().st_size} bytes)", flush=True)
                else:
                    print(f"[WARNING] No usable partial file yet: {part}", flush=True)

                print("[INFO] Will retry automatically.", flush=True)

    raise RuntimeError(
        f"All download attempts failed for {outfile}. "
        f"Partial file is kept if present: {part}. "
        f"Last error: {last_error}"
    )

def gunzip_file(gz_file, outfile, force=False):
    gz_file = Path(gz_file)
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    if is_nonempty_file(outfile) and not force:
        print(f"[SKIP] Existing decompressed file found: {outfile}", flush=True)
        return outfile

    if not is_nonempty_file(gz_file):
        raise RuntimeError(f"Compressed file missing or empty: {gz_file}")

    tmp = Path(str(outfile) + ".part")
    if tmp.exists():
        tmp.unlink()

    print(f"[GUNZIP] {gz_file} -> {outfile}", flush=True)

    try:
        with gzip.open(gz_file, "rb") as f_in, open(tmp, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        if not is_nonempty_file(tmp):
            raise RuntimeError(f"Decompressed file is empty: {tmp}")

        tmp.rename(outfile)

    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise

    print(f"[OK] Decompressed: {outfile}", flush=True)
    return outfile


def safe_extract_tar(tar, path):
    path = Path(path).resolve()

    for member in tar.getmembers():
        member_path = (path / member.name).resolve()
        try:
            member_path.relative_to(path)
        except ValueError:
            raise RuntimeError(f"Unsafe path in tar archive: {member.name}")

    tar.extractall(path=path)


def extract_tar_gz(tar_gz, outdir, force=False):
    tar_gz = Path(tar_gz)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    marker = outdir / ".eggnog_taxa_extracted.ok"

    if marker.exists() and not force:
        print(f"[SKIP] eggNOG taxa archive already extracted: {outdir}", flush=True)
        return outdir

    if not is_nonempty_file(tar_gz):
        print(f"[WARNING] Taxa archive missing or empty, skip extraction: {tar_gz}", flush=True)
        return outdir

    print(f"[TAR] {tar_gz} -> {outdir}", flush=True)

    try:
        with tarfile.open(tar_gz, "r:gz") as tar:
            safe_extract_tar(tar, outdir)
    except tarfile.ReadError as e:
        print(f"[WARNING] Failed to extract taxa archive: {tar_gz}", flush=True)
        print(f"[WARNING] {e}", flush=True)
        print("[WARNING] Continue if eggnog.db and eggnog_proteins.dmnd are available.", flush=True)
        return outdir

    marker.write_text("ok\n")
    print(f"[OK] Extracted: {tar_gz}", flush=True)
    return outdir


def check_pfam_index(pfam_hmm):
    pfam_hmm = Path(pfam_hmm)
    index_files = [
        Path(str(pfam_hmm) + ".h3f"),
        Path(str(pfam_hmm) + ".h3i"),
        Path(str(pfam_hmm) + ".h3m"),
        Path(str(pfam_hmm) + ".h3p"),
    ]
    return all(is_nonempty_file(x) for x in index_files)


def prepare_diamond_swissprot(db_dir, force=False):
    db_dir = Path(db_dir)
    diamond_dir = db_dir / "diamond"
    diamond_dir.mkdir(parents=True, exist_ok=True)

    fasta_gz = diamond_dir / "uniprot_sprot.fasta.gz"
    fasta = diamond_dir / "uniprot_sprot.fasta"
    dmnd = diamond_dir / "uniprot_sprot.dmnd"

    if is_nonempty_file(dmnd) and not force:
        print(f"[OK] DIAMOND database already exists: {dmnd}", flush=True)
        return dmnd

    if shutil.which("diamond") is None:
        raise RuntimeError("diamond not found. Run: protein-annotate install-tools")

    download_file(UNIPROT_SPROT_URL, fasta_gz, force=force)
    gunzip_file(fasta_gz, fasta, force=force)

    cmd = [
        "diamond",
        "makedb",
        "--in", str(fasta),
        "--db", str(diamond_dir / "uniprot_sprot"),
    ]

    run_command(cmd)

    if not is_nonempty_file(dmnd):
        raise RuntimeError(f"DIAMOND makedb finished, but database not found or empty: {dmnd}")

    print(f"[OK] DIAMOND db ready: {dmnd}", flush=True)
    return dmnd


def prepare_pfam_db(db_dir, force=False):
    db_dir = Path(db_dir)
    pfam_dir = db_dir / "pfam"
    pfam_dir.mkdir(parents=True, exist_ok=True)

    pfam_gz = pfam_dir / "Pfam-A.hmm.gz"
    pfam_hmm = pfam_dir / "Pfam-A.hmm"

    if is_nonempty_file(pfam_hmm) and check_pfam_index(pfam_hmm) and not force:
        print(f"[OK] Pfam database and hmmpress index already exist: {pfam_hmm}", flush=True)
        return pfam_hmm

    if shutil.which("hmmpress") is None:
        raise RuntimeError("hmmpress not found. Run: protein-annotate install-tools")

    download_file(PFAM_A_HMM_URL, pfam_gz, force=force)
    gunzip_file(pfam_gz, pfam_hmm, force=force)

    if not check_pfam_index(pfam_hmm) or force:
        run_command(["hmmpress", str(pfam_hmm)])

    if not check_pfam_index(pfam_hmm):
        raise RuntimeError(f"hmmpress finished, but Pfam index files are missing: {pfam_hmm}")

    print(f"[OK] Pfam HMM ready: {pfam_hmm}", flush=True)
    return pfam_hmm


def eggnog_file_urls(filename):
    return [f"{base}/{filename}" for base in EGGNOG_BASE_URLS]


def clean_zero_size_eggnog_files(eggnog_dir):
    eggnog_dir = Path(eggnog_dir)
    clean_zero_size_files([
        eggnog_dir / "eggnog.db",
        eggnog_dir / "eggnog.db.gz",
        eggnog_dir / "eggnog_proteins.dmnd",
        eggnog_dir / "eggnog_proteins.dmnd.gz",
        eggnog_dir / "eggnog.taxa.tar.gz",
        eggnog_dir / "eggnog.db.gz.part",
        eggnog_dir / "eggnog_proteins.dmnd.gz.part",
        eggnog_dir / "eggnog.taxa.tar.gz.part",
    ])


def decompress_eggnog_files(eggnog_dir, force=False):
    eggnog_dir = Path(eggnog_dir)
    eggnog_dir.mkdir(parents=True, exist_ok=True)

    eggnog_db_gz = eggnog_dir / "eggnog.db.gz"
    eggnog_db = eggnog_dir / "eggnog.db"

    proteins_gz = eggnog_dir / "eggnog_proteins.dmnd.gz"
    proteins_dmnd = eggnog_dir / "eggnog_proteins.dmnd"

    taxa_tar = eggnog_dir / "eggnog.taxa.tar.gz"

    if is_nonempty_file(eggnog_db_gz):
        gunzip_file(eggnog_db_gz, eggnog_db, force=force)

    if is_nonempty_file(proteins_gz):
        gunzip_file(proteins_gz, proteins_dmnd, force=force)

    if is_nonempty_file(taxa_tar):
        extract_tar_gz(taxa_tar, eggnog_dir, force=force)

    return eggnog_dir


def eggnog_db_exists(eggnog_dir):
    eggnog_dir = Path(eggnog_dir)
    eggnog_db = eggnog_dir / "eggnog.db"
    proteins_dmnd = eggnog_dir / "eggnog_proteins.dmnd"

    return is_nonempty_file(eggnog_db) and is_nonempty_file(proteins_dmnd)


def prepare_eggnog_db(db_dir, force=False):
    """
    Prepare eggNOG-mapper v2 database without calling interactive download_eggnog_data.py.
    This avoids old hard-coded eggnogdb.embl.de URLs.
    """
    db_dir = Path(db_dir)
    eggnog_dir = db_dir / "eggnog"
    eggnog_dir.mkdir(parents=True, exist_ok=True)

    eggnog_db_gz = eggnog_dir / "eggnog.db.gz"
    eggnog_db = eggnog_dir / "eggnog.db"

    proteins_gz = eggnog_dir / "eggnog_proteins.dmnd.gz"
    proteins_dmnd = eggnog_dir / "eggnog_proteins.dmnd"

    taxa_tar = eggnog_dir / "eggnog.taxa.tar.gz"

    clean_zero_size_eggnog_files(eggnog_dir)

    if eggnog_db_exists(eggnog_dir) and not force:
        print(f"[OK] eggNOG data dir already exists and is decompressed: {eggnog_dir}", flush=True)
        return eggnog_dir

    print("[INFO] Download eggNOG v5.0.2 database directly from eggnog5.embl.de", flush=True)

    if force:
        for f in [eggnog_db, proteins_dmnd, eggnog_db_gz, proteins_gz, taxa_tar]:
            if f.exists():
                print(f"[REMOVE] Existing eggNOG file because --force: {f}", flush=True)
                f.unlink()

    if not is_nonempty_file(eggnog_db):
        download_file(eggnog_file_urls("eggnog.db.gz"), eggnog_db_gz, force=force)
        gunzip_file(eggnog_db_gz, eggnog_db, force=True)

    if not is_nonempty_file(proteins_dmnd):
        download_file(eggnog_file_urls("eggnog_proteins.dmnd.gz"), proteins_gz, force=force)
        gunzip_file(proteins_gz, proteins_dmnd, force=True)

    # taxa is useful but not allowed to break the whole pipeline if unavailable.
    try:
        if not is_nonempty_file(taxa_tar) or force:
            download_file(eggnog_file_urls("eggnog.taxa.tar.gz"), taxa_tar, force=force)
        extract_tar_gz(taxa_tar, eggnog_dir, force=force)
    except Exception as e:
        print(f"[WARNING] eggNOG taxa download/extraction failed: {e}", flush=True)
        print("[WARNING] Continue because main files eggnog.db and eggnog_proteins.dmnd are present.", flush=True)

    clean_zero_size_eggnog_files(eggnog_dir)

    if not eggnog_db_exists(eggnog_dir):
        raise RuntimeError(
            "eggNOG database is incomplete. Required non-empty files are missing: "
            f"{eggnog_db} and {proteins_dmnd}. "
            "Please check network access to eggnog5.embl.de."
        )

    print(f"[OK] eggNOG data dir ready: {eggnog_dir}", flush=True)
    return eggnog_dir


def prepare_databases(
    db_dir,
    prepare_diamond=True,
    prepare_pfam=True,
    prepare_eggnog=True,
    force=False,
):
    db_dir = Path(db_dir)
    db_dir.mkdir(parents=True, exist_ok=True)

    print("===== Protein Annotator: prepare-db =====", flush=True)
    print(f"[DB DIR] {db_dir}", flush=True)

    result = {}

    if prepare_diamond:
        print("\n===== Prepare DIAMOND Swiss-Prot database =====", flush=True)
        result["diamond_db"] = str(prepare_diamond_swissprot(db_dir, force=force))
    else:
        print("[SKIP] DIAMOND database preparation", flush=True)

    if prepare_pfam:
        print("\n===== Prepare Pfam database =====", flush=True)
        result["pfam_hmm"] = str(prepare_pfam_db(db_dir, force=force))
    else:
        print("[SKIP] Pfam database preparation", flush=True)

    if prepare_eggnog:
        print("\n===== Prepare eggNOG-mapper database =====", flush=True)
        result["eggnog_data_dir"] = str(prepare_eggnog_db(db_dir, force=force))
    else:
        print("[SKIP] eggNOG database preparation", flush=True)

    manifest = db_dir / "protein_annotator_db_paths.sh"

    # Preserve existing paths if some components are skipped.
    old_lines = {}
    if manifest.exists():
        for line in manifest.read_text().splitlines():
            if line.startswith("export PROTEIN_ANNOTATOR_DIAMOND_DB="):
                old_lines["diamond_db"] = line
            elif line.startswith("export PROTEIN_ANNOTATOR_PFAM_HMM="):
                old_lines["pfam_hmm"] = line
            elif line.startswith("export PROTEIN_ANNOTATOR_EGGNOG_DATA_DIR="):
                old_lines["eggnog_data_dir"] = line

    with open(manifest, "w", encoding="utf-8") as out:
        out.write("# Auto-generated by protein-annotate prepare-db\n")

        if "diamond_db" in result:
            out.write(f'export PROTEIN_ANNOTATOR_DIAMOND_DB="{result["diamond_db"]}"\n')
        elif "diamond_db" in old_lines:
            out.write(old_lines["diamond_db"] + "\n")

        if "pfam_hmm" in result:
            out.write(f'export PROTEIN_ANNOTATOR_PFAM_HMM="{result["pfam_hmm"]}"\n')
        elif "pfam_hmm" in old_lines:
            out.write(old_lines["pfam_hmm"] + "\n")

        if "eggnog_data_dir" in result:
            out.write(f'export PROTEIN_ANNOTATOR_EGGNOG_DATA_DIR="{result["eggnog_data_dir"]}"\n')
        elif "eggnog_data_dir" in old_lines:
            out.write(old_lines["eggnog_data_dir"] + "\n")

    print("\n===== Database paths =====", flush=True)
    for k, v in result.items():
        print(f"{k}: {v}", flush=True)

    print(f"\n[OK] Path manifest written to: {manifest}", flush=True)
    print(f"To load paths later, run:\n  source {manifest}", flush=True)

    return result
