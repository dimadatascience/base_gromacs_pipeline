import os

with open("temp.txt", "r") as f:
    lines = [line.strip() for line in f]
    ligand_path=lines[0]
    crystal_path=lines[1]
    crystal_code=lines[2]
    ligand_common_name=lines[3]
    ligand_temp_path=lines[4]
    ligand_code=lines[5]
    renamed_crystal_path_pdb_file=lines[6]
    renamed_ligand_path_pdb_file=lines[7]
    renamed_ligand_pdb_file=lines[8]

os.chdir(ligand_path) # change the current working directory

os.system(f"acpype -i {renamed_ligand_pdb_file} -b {ligand_code}")