from pathlib import Path
import os
import subprocess
import sys
code_root = Path(__file__).resolve().parents[1]
input_root = Path(os.environ.get("GENES_INPUT_ROOT", code_root / "data"))
package_root = Path(os.environ.get("DRUGREFLECTOR_ROOT", input_root / "external_assets"))
checkpoints = Path(os.environ.get("DRUGREFLECTOR_CHECKPOINTS", input_root / "external_assets" / "checkpoints"))
project = input_root / "bulk_training" / "09_DrugReflector"
prediction = project / "02_prediction" / "01_DrugReflector_OA_to_Control_all_compounds.csv"
prediction.parent.mkdir(parents=True, exist_ok=True)
subprocess.run([sys.executable, str(package_root / "drugreflector" / "predict.py"), str(project / "01_input" / "02_OA_to_Control_DrugReflector_input.h5ad"),
    "--model1", str(checkpoints / "model_fold_0.pt"), "--model2", str(checkpoints / "model_fold_1.pt"), "--model3", str(checkpoints / "model_fold_2.pt"),
    "--top-n", "0", "-o", str(prediction)], check=True)
