from pathlib import Path
import subprocess

ROOT = Path("/home/yourname/projects/9TI4_NDUFS3_NDUFA5_VU0418947_2_MD_5ns_gate_20260906")
GMX = Path("/home/yourname/opt/gene-drug-automd/env/bin/gmx")
BUILD = ROOT / "build"
MDP = ROOT / "mdp"
RUN = ROOT / "run"

PARAMETERS = {
    "em.mdp": """integrator = steep
emtol = 1000.0
emstep = 0.01
nsteps = 50000
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
pbc = xyz
constraints = none
constraint-algorithm = lincs
nstlist = 20
""",
    "nvt.mdp": """define = -DPOSRES1000
integrator = md
dt = 0.002
nsteps = 50000
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
pme-order = 4
fourierspacing = 0.12
tcoupl = V-rescale
tc-grps = System
tau-t = 0.1
ref-t = 300
pcoupl = no
constraints = h-bonds
constraint-algorithm = lincs
lincs-iter = 2
lincs-order = 4
continuation = no
gen-vel = yes
gen-temp = 300
gen-seed = 6501
nstxout-compressed = 1000
nstenergy = 1000
nstlog = 1000
nstlist = 20
pbc = xyz
""",
    "npt.mdp": """define = -DPOSRES500
integrator = md
dt = 0.002
nsteps = 250000
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
pme-order = 4
fourierspacing = 0.12
tcoupl = V-rescale
tc-grps = System
tau-t = 0.1
ref-t = 300
pcoupl = Parrinello-Rahman
pcoupltype = isotropic
tau-p = 2.0
ref-p = 1.0
compressibility = 4.5e-5
constraints = h-bonds
constraint-algorithm = lincs
lincs-iter = 2
lincs-order = 4
continuation = yes
gen-vel = no
nstxout-compressed = 1000
nstenergy = 1000
nstlog = 1000
nstlist = 20
pbc = xyz
refcoord-scaling = all
""",
    "md.mdp": """integrator = md
dt = 0.002
nsteps = 2500000
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
pme-order = 4
fourierspacing = 0.12
tcoupl = V-rescale
tc-grps = System
tau-t = 0.1
ref-t = 300
pcoupl = Parrinello-Rahman
pcoupltype = isotropic
tau-p = 2.0
ref-p = 1.0
compressibility = 4.5e-5
constraints = h-bonds
constraint-algorithm = lincs
lincs-iter = 2
lincs-order = 4
continuation = yes
gen-vel = no
nstxout-compressed = 5000
nstenergy = 5000
nstlog = 5000
nstlist = 20
pbc = xyz
comm-mode = Linear
nstcomm = 100
""",
}


def run(command):
    subprocess.run([str(x) for x in command], cwd=RUN, check=True)


MDP.mkdir(parents=True, exist_ok=True)
RUN.mkdir(parents=True, exist_ok=True)
for name, text in PARAMETERS.items():
    (MDP / name).write_text(text, encoding="utf-8")

run([GMX, "grompp", "-f", MDP / "em.mdp", "-c", BUILD / "system.gro", "-p", BUILD / "system.top", "-o", "em.tpr"])
run([GMX, "mdrun", "-deffnm", "em", "-ntmpi", 1, "-ntomp", 8])
run([GMX, "grompp", "-f", MDP / "nvt.mdp", "-c", "em.gro", "-r", BUILD / "system.gro", "-p", BUILD / "system.top", "-o", "nvt.tpr"])
run([GMX, "mdrun", "-deffnm", "nvt", "-ntmpi", 1, "-nb", "gpu", "-pme", "gpu", "-bonded", "gpu", "-update", "gpu", "-pin", "on"])
run([GMX, "grompp", "-f", MDP / "npt.mdp", "-c", "nvt.gro", "-t", "nvt.cpt", "-r", BUILD / "system.gro", "-p", BUILD / "system.top", "-o", "npt.tpr"])
run([GMX, "mdrun", "-deffnm", "npt", "-ntmpi", 1, "-nb", "gpu", "-pme", "gpu", "-bonded", "gpu", "-update", "gpu", "-pin", "on"])
run([GMX, "grompp", "-f", MDP / "md.mdp", "-c", "npt.gro", "-t", "npt.cpt", "-p", BUILD / "system.top", "-o", "md.tpr"])
run([GMX, "mdrun", "-deffnm", "md", "-s", "md.tpr", "-ntmpi", 1, "-nb", "gpu", "-pme", "gpu", "-bonded", "gpu", "-update", "gpu", "-pin", "on"])
run([GMX, "convert-tpr", "-s", "md.tpr", "-extend", 95000, "-o", "md_100ns.tpr"])
run([GMX, "mdrun", "-deffnm", "md", "-s", "md_100ns.tpr", "-cpi", "md.cpt", "-append", "-ntmpi", 1, "-nb", "gpu", "-pme", "gpu", "-bonded", "gpu", "-update", "gpu", "-pin", "on"])
