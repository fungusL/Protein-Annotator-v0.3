from .checks import doctor, install_missing_tools, prepare_pfam
from .pipeline import annotate_proteome
from .tables import filter_annotation_tables
from .db import prepare_databases

__version__ = "0.3.0"

__all__ = [
    "doctor",
    "install_missing_tools",
    "prepare_pfam",
    "annotate_proteome",
    "filter_annotation_tables",
    "prepare_databases",
]
