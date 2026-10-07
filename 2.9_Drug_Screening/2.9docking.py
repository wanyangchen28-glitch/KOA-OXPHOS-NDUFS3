from pathlib import Path
import subprocess

BIN = Path("/home/yourname/opt/gene-drug-automd/env/bin")
RECEPTOR = Path("/home/yourname/projects/9TI4_NDUFS3_interface_screen_20260905/receptor/9TI4_NDUFS3_NDUFA5_NDUFS2.pdb")
ROOT = Path("/home/yourname/projects/9TI4_NDUFS3_interface_screen_20260905/final_docking")
CENTER = (179.274, 234.805, 289.337)
SIZE = (23.228, 23.228, 23.228)
SEEDS = (9701, 9702, 9703)
COMPOUNDS = {
    "Triciribine_CID16759157": Path("/home/yourname/projects/9TI4_CID16759157/output/02_ligand"),
    "ST4066738_CID1043828": Path("/home/yourname/projects/9TI4_ST4066738_CID1043828/output/02_ligand"),
    "Calmidazolium_CID2531": Path("/home/yourname/projects/9TI4_Calmidazolium_CID2531/output/02_ligand"),
    "I606051_CID23891056": Path("/home/yourname/projects/9TI4_I606051_CID23891056/output/02_ligand"),
    "VU0418947_2_CID46742354": Path("/home/yourname/projects/9TI4_VU0418947_2_CID46742354/output/02_ligand"),
}


def run(command, cwd, log):
    result = subprocess.run([str(x) for x in command], cwd=cwd, text=True, capture_output=True)
    log.write_text(result.stdout + result.stderr, encoding="utf-8")
    result.check_returncode()


ROOT.mkdir(parents=True, exist_ok=True)
receptor = ROOT / "receptor.pdbqt"
run(
    [BIN / "mk_prepare_receptor.py", "--read_pdb", RECEPTOR, "-a", "-p", receptor],
    ROOT,
    ROOT / "receptor_preparation.log",
)

for compound, source in COMPOUNDS.items():
    for state_sdf in sorted(source.glob("state_*.sdf")):
        state = ROOT / compound / state_sdf.stem
        state.mkdir(parents=True, exist_ok=True)
        ligand = state / "ligand.pdbqt"
        run(
            [BIN / "mk_prepare_ligand.py", "-i", state_sdf, "-o", ligand],
            state,
            state / "ligand_preparation.log",
        )
        for seed in SEEDS:
            out = state / f"seed_{seed}"
            out.mkdir(exist_ok=True)
            poses = out / "poses.pdbqt"
            run(
                [
                    BIN / "vina",
                    "--receptor", receptor,
                    "--ligand", ligand,
                    "--center_x", CENTER[0],
                    "--center_y", CENTER[1],
                    "--center_z", CENTER[2],
                    "--size_x", SIZE[0],
                    "--size_y", SIZE[1],
                    "--size_z", SIZE[2],
                    "--exhaustiveness", 32,
                    "--num_modes", 10,
                    "--energy_range", 3,
                    "--seed", seed,
                    "--out", poses,
                ],
                out,
                out / "vina.log",
            )
            run(
                [BIN / "mk_export.py", poses, "-s", out / "poses.sdf"],
                out,
                out / "export.log",
            )
