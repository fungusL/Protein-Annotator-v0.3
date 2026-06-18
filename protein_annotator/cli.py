import argparse

from .checks import doctor, install_missing_tools, prepare_pfam
from .pipeline import annotate_proteome
from .tables import filter_annotation_tables
from .db import prepare_databases


def build_parser():
    parser = argparse.ArgumentParser(
        prog="protein-annotate",
        description="Protein functional annotation workflow.",
    )

    sub = parser.add_subparsers(dest="command", required=True)

    p_doctor = sub.add_parser("doctor", help="Check tools and database paths.")
    p_doctor.add_argument("--diamond-db", default=None)
    p_doctor.add_argument("--pfam-hmm", default=None)
    p_doctor.add_argument("--eggnog-data-dir", default=None)

    p_install = sub.add_parser("install-tools", help="Install missing tools with mamba/conda.")
    p_install.add_argument("--manager", default="mamba", choices=["mamba", "conda"])
    p_install.add_argument("--dry-run", action="store_true")

    p_pfam = sub.add_parser("prepare-pfam", help="Run hmmpress for an existing Pfam-A.hmm.")
    p_pfam.add_argument("--pfam-hmm", required=True)

    p_db = sub.add_parser("prepare-db", help="Download/build DIAMOND, Pfam and eggNOG databases if missing.")
    p_db.add_argument("--db-dir", required=True)
    p_db.add_argument("--skip-diamond", action="store_true")
    p_db.add_argument("--skip-pfam", action="store_true")
    p_db.add_argument("--skip-eggnog", action="store_true")
    p_db.add_argument("--force", action="store_true")

    p_filter = sub.add_parser("filter", help="Filter an existing final merged annotation table.")
    p_filter.add_argument("--annotation-tsv", required=True)
    p_filter.add_argument("--outdir", required=True)
    p_filter.add_argument("--prefix", required=True)

    p_run = sub.add_parser("run", help="Run full annotation pipeline from protein fasta.")
    p_run.add_argument("--protein-fasta", required=True)
    p_run.add_argument("--outdir", required=True)
    p_run.add_argument("--prefix", required=True)

    p_run.add_argument("--diamond-db", default=None)
    p_run.add_argument("--pfam-hmm", default=None)
    p_run.add_argument("--eggnog-data-dir", default=None)

    p_run.add_argument("--auto-prepare-db", action="store_true")
    p_run.add_argument("--db-dir", default=None)

    p_run.add_argument("--threads", type=int, default=8)
    p_run.add_argument("--force", action="store_true")
    p_run.add_argument("--skip-doctor", action="store_true")

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "doctor":
        doctor(args.diamond_db, args.pfam_hmm, args.eggnog_data_dir)

    elif args.command == "install-tools":
        install_missing_tools(manager=args.manager, dry_run=args.dry_run)

    elif args.command == "prepare-pfam":
        prepare_pfam(args.pfam_hmm)

    elif args.command == "prepare-db":
        prepare_databases(
            db_dir=args.db_dir,
            prepare_diamond=not args.skip_diamond,
            prepare_pfam=not args.skip_pfam,
            prepare_eggnog=not args.skip_eggnog,
            force=args.force,
        )

    elif args.command == "filter":
        filter_annotation_tables(args.annotation_tsv, args.outdir, args.prefix)

    elif args.command == "run":
        diamond_db = args.diamond_db
        pfam_hmm = args.pfam_hmm
        eggnog_data_dir = args.eggnog_data_dir

        if args.auto_prepare_db:
            if not args.db_dir:
                parser.error("--auto-prepare-db requires --db-dir")

            db_paths = prepare_databases(
                db_dir=args.db_dir,
                prepare_diamond=(diamond_db is None),
                prepare_pfam=(pfam_hmm is None),
                prepare_eggnog=(eggnog_data_dir is None),
                force=False,
            )

            diamond_db = diamond_db or db_paths.get("diamond_db")
            pfam_hmm = pfam_hmm or db_paths.get("pfam_hmm")
            eggnog_data_dir = eggnog_data_dir or db_paths.get("eggnog_data_dir")

        if not diamond_db:
            parser.error("Missing --diamond-db. Provide it or use --auto-prepare-db --db-dir.")
        if not pfam_hmm:
            parser.error("Missing --pfam-hmm. Provide it or use --auto-prepare-db --db-dir.")

        annotate_proteome(
            protein_fasta=args.protein_fasta,
            outdir=args.outdir,
            prefix=args.prefix,
            diamond_db=diamond_db,
            pfam_hmm=pfam_hmm,
            eggnog_data_dir=eggnog_data_dir,
            threads=args.threads,
            force=args.force,
            skip_doctor=args.skip_doctor,
        )


if __name__ == "__main__":
    main()
