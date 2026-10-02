from pathlib import Path
import os
CODE_ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = Path(os.environ.get("GENES_INPUT_ROOT", CODE_ROOT / "data"))
OUTPUT_ROOT = Path(os.environ.get("GENES_OUTPUT_ROOT", CODE_ROOT / "results"))
from pathlib import Path
import pandas as pd
PROJECT_DIR = INPUT_ROOT / 'bulk_training'
PREDICTION_DIR = PROJECT_DIR / '09_DrugReflector' / '02_prediction'
TOP50_FILE = PREDICTION_DIR / '02_DrugReflector_OA_to_Control_top50.csv'
LINCS_INFO_FILE = PROJECT_DIR / '09_DrugReflector' / '00_reference' / 'GSE92742_Broad_LINCS_pert_info.txt.gz'

def main():
    top50 = pd.read_csv(TOP50_FILE)
    lincs = pd.read_csv(LINCS_INFO_FILE, sep='\t', compression='gzip')
    lincs = lincs.drop_duplicates(subset='pert_id', keep='first')
    annotated = top50.merge(lincs, left_on='compound_id', right_on='pert_id', how='left')
    annotated['has_readable_compound_name'] = annotated['pert_iname'].notna() & ~annotated['pert_iname'].astype(str).str.startswith('BRD-')
    annotated.to_csv(PREDICTION_DIR / '05_DrugReflector_top50_LINCS_annotated.csv', index=False)
    named_unique = annotated.loc[(annotated['pert_type'] == 'trt_cp') & annotated['has_readable_compound_name']].sort_values('rank', kind='stable').drop_duplicates(subset='pert_iname', keep='first').loc[:, ['compound_id', 'rank', 'rank_model_zero_based', 'pert_iname', 'logit', 'probability', 'pubchem_cid', 'canonical_smiles', 'inchi_key', 'pert_type']]
    named_unique.to_csv(PREDICTION_DIR / '06_DrugReflector_top50_named_unique.csv', index=False)
    summary = pd.DataFrame({'metric': ['top50_total', 'LINCS_annotated', 'small_molecule_entries', 'named_unique_small_molecules'], 'value': [len(annotated), int(annotated['pert_iname'].notna().sum()), int((annotated['pert_type'] == 'trt_cp').sum()), len(named_unique)]})
    summary.to_csv(PREDICTION_DIR / '07_DrugReflector_top50_annotation_summary.csv', index=False)
    print(named_unique.to_string(index=False))
if __name__ == '__main__':
    main()
