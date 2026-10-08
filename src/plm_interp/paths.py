"""Project-relative paths, resolved from this file so scripts work from any cwd."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

MODELS_DIR = ROOT / "models"
DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
INTERIM_DIR = DATA_DIR / "interim"
PROCESSED_DIR = DATA_DIR / "processed"
OUTPUTS_DIR = ROOT / "outputs"
RUNS_DIR = OUTPUTS_DIR / "runs"
FIGURES_DIR = OUTPUTS_DIR / "figures"
EXPLANATIONS_DIR = OUTPUTS_DIR / "explanations"
CONFIGS_DIR = ROOT / "configs"

PROTEINGYM_DIR = RAW_DIR / "proteingym"
DMS_DIR = PROTEINGYM_DIR / "dms" / "DMS_ProteinGym_substitutions"
CLINICAL_DIR = PROTEINGYM_DIR / "clinical"
REFERENCE_DIR = PROTEINGYM_DIR / "reference"
BENCHMARK_DIR = PROTEINGYM_DIR / "benchmark"
STRUCTURES_DIR = RAW_DIR / "structures" / "ProteinGym_AF2_structures"
MSA_ZIP = RAW_DIR / "msas" / "DMS_msa_files.zip"
MSA_DIR = INTERIM_DIR / "msas"
ANNOTATIONS_DIR = PROCESSED_DIR / "annotations"
SCORES_DIR = PROCESSED_DIR / "scores"
