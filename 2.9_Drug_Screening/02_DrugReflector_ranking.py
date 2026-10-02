from pathlib import Path
import os
CODE_ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = Path(os.environ.get("GENES_INPUT_ROOT", CODE_ROOT / "data"))
OUTPUT_ROOT = Path(os.environ.get("GENES_OUTPUT_ROOT", CODE_ROOT / "results"))
from pathlib import Path
import pandas as pd
PROJECT_DIR = INPUT_ROOT / 'bulk_training'
PREDICTION_FILE = PROJECT_DIR / '09_DrugReflector' / '02_prediction' / '01_DrugReflector_OA_to_Control_all_compounds.csv'
OUTPUT_DIR = PROJECT_DIR / '09_DrugReflector' / '02_prediction'

def main():
    raw = pd.read_csv(PREDICTION_FILE, header=[0, 1], index_col=0)
    transition = raw.columns.get_level_values(1)[0]
    result = pd.DataFrame({'compound_id': raw.index, 'rank_model_zero_based': pd.to_numeric(raw['rank', transition]), 'logit': pd.to_numeric(raw['logit', transition]), 'probability': pd.to_numeric(raw['prob', transition])}).sort_values('rank_model_zero_based', kind='stable')
    result.insert(1, 'rank', result['rank_model_zero_based'] + 1)
    result.head(50).to_csv(OUTPUT_DIR / '02_DrugReflector_OA_to_Control_top50.csv', index=False)
    result.head(10).to_csv(OUTPUT_DIR / '03_DrugReflector_OA_to_Control_top10.csv', index=False)
    summary = pd.DataFrame({'metric': ['transition', 'n_ranked_compounds', 'model_rank_indexing', 'top_compound_id'], 'value': [transition, len(result), 'raw model rank starts at 0; exported rank starts at 1', result.iloc[0]['compound_id']]})
    summary.to_csv(OUTPUT_DIR / '04_DrugReflector_prediction_summary.csv', index=False)
    print(result.head(10).to_string(index=False))
if __name__ == '__main__':
    main()
