# Importa il modulo ChimeraX
from chimerax.core.commands import run
import os

with open("temp.txt", "r") as f:
    lines = [line.strip() for line in f]
    ligand_path=lines[0]
    crystal_path=lines[1]
    crystal_code=lines[2]
    ligand_common_name=lines[3]
    ligand_temp_path=lines[4]
    ligand_code=lines[5]

# Carica il file di struttura (può essere un PDB, MOL2, ecc.)
run(session, f"open {ligand_temp_path}")

# Aggiungi gli idrogeni
run(session, "addh")

# Salva la struttura con gli idrogeni aggiunti in un nuovo file PDB
run(session, f"save {ligand_temp_path}")

os.rename(f"{crystal_path}/{crystal_code}.pdb", f"{crystal_path}/{ligand_code}_{crystal_code}.pdb")
os.rename(ligand_temp_path, f"{ligand_path}/{ligand_code}_{ligand_common_name}.pdb")

renamed_crystal_path_pdb_file=f"{crystal_path}/{ligand_code}_{crystal_code}.pdb"
renamed_ligand_path_pdb_file=f"{ligand_path}/{ligand_code}_{ligand_common_name}.pdb"
renamed_ligand_pdb_file=f"{ligand_code}_{ligand_common_name}.pdb"

with open("temp.txt", "w") as f:
    f.write(f"{ligand_path}\n{crystal_path}\n{crystal_code}\n{ligand_common_name}\n{ligand_temp_path}\n{ligand_code}\n{renamed_crystal_path_pdb_file}\n{renamed_ligand_path_pdb_file}\n{renamed_ligand_pdb_file}")
    
run(session, "quit")


