from pathlib import Path
import os
CODE_ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = Path(os.environ.get("GENES_INPUT_ROOT", CODE_ROOT / "data"))
OUTPUT_ROOT = Path(os.environ.get("GENES_OUTPUT_ROOT", CODE_ROOT / "results"))
from pathlib import Path
import importlib.util
import numpy as np
import pandas as pd
from anndata import AnnData
PROJECT_DIR = INPUT_ROOT / 'bulk_training'
EXPRESSION_FILE = PROJECT_DIR / '00_preprocess' / '06_expanded_training_ComBat_HGNC_expression.csv'
METADATA_FILE = PROJECT_DIR / '00_preprocess' / '07_expanded_training_sample_info.tsv'
DRUGREFLECTOR_UTILS = INPUT_ROOT / 'external_assets/utils.py'
OUTPUT_DIR = PROJECT_DIR / '09_DrugReflector' / '01_input'
SAMPLE_COLUMN = 'Sample'
GROUP_COLUMN = 'Group'
FROM_GROUP = 'OA'
TO_GROUP = 'Control'

def load_drugreflector_utils(path: Path):
    spec = importlib.util.spec_from_file_location('drugreflector_utils', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f'Expression matrix: {EXPRESSION_FILE}')
    print(f'Sample metadata: {METADATA_FILE}')
    print(f'Output directory: {OUTPUT_DIR}')
    expression = pd.read_csv(EXPRESSION_FILE, index_col=0)
    metadata = pd.read_csv(METADATA_FILE, sep='\t')
    sample_order = metadata[SAMPLE_COLUMN].astype(str).tolist()
    expression.columns = expression.columns.astype(str)
    matched_samples = [sample for sample in sample_order if sample in expression.columns]
    expression = expression.loc[:, matched_samples]
    expression.index = expression.index.astype(str).str.strip().str.upper()
    expression = expression.loc[expression.index != '']
    expression = expression.groupby(level=0, sort=False).mean()
    metadata = metadata.copy()
    metadata[SAMPLE_COLUMN] = metadata[SAMPLE_COLUMN].astype(str)
    metadata = metadata.set_index(SAMPLE_COLUMN).loc[matched_samples].copy()
    adata = AnnData(X=expression.T.to_numpy(dtype=np.float32), obs=metadata, var=pd.DataFrame(index=expression.index))
    drugreflector_utils = load_drugreflector_utils(DRUGREFLECTOR_UTILS)
    vscores = drugreflector_utils.compute_vscores_adata(adata=adata, group_col=GROUP_COLUMN, group1_value=FROM_GROUP, group2_value=TO_GROUP)
    transition_name = f'{GROUP_COLUMN}:{FROM_GROUP}->{TO_GROUP}'
    vscore_table = pd.DataFrame({'Gene': vscores.index, 'v_score': vscores.values})
    vscore_csv = OUTPUT_DIR / '01_OA_to_Control_vscore_all_genes.csv'
    vscore_table.to_csv(vscore_csv, index=False)
    vscore_adata = AnnData(X=vscores.to_numpy(dtype=np.float32).reshape(1, -1), obs=pd.DataFrame({'transition': [transition_name], 'from_group': [FROM_GROUP], 'to_group': [TO_GROUP]}, index=[transition_name]), var=pd.DataFrame(index=vscores.index))
    vscore_h5ad = OUTPUT_DIR / '02_OA_to_Control_DrugReflector_input.h5ad'
    vscore_adata.write_h5ad(vscore_h5ad)
    group_counts = metadata[GROUP_COLUMN].value_counts().to_dict()
    summary = pd.DataFrame({'metric': ['input_expression_file', 'input_metadata_file', 'transition', 'matched_samples', 'OA_samples', 'Control_samples', 'genes_in_vscore_input', 'vscore_formula'], 'value': [str(EXPRESSION_FILE), str(METADATA_FILE), transition_name, len(matched_samples), group_counts.get(FROM_GROUP, 0), group_counts.get(TO_GROUP, 0), len(vscores), '(mean_Control - mean_OA) / sqrt(var_OA + var_Control)']})
    summary_file = OUTPUT_DIR / '03_OA_to_Control_vscore_summary.csv'
    summary.to_csv(summary_file, index=False)
    print(f'Prepared {len(vscores)} genes from {len(matched_samples)} samples.')
    print(f'OA={group_counts.get(FROM_GROUP, 0)}; Control={group_counts.get(TO_GROUP, 0)}')
    print(f'CSV input: {vscore_csv}')
    print(f'H5AD input: {vscore_h5ad}')
if __name__ == '__main__':
    main()
