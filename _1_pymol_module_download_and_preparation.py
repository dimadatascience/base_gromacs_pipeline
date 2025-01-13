import pymol
from pymol import cmd
import urllib.request
from collections import Counter

def download_crystal_pdb_in_crystal_path(pdb_code, crystal_path):
    # URL del PDB server
    pdb_url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
    
    # Scarica il file PDB
    urllib.request.urlretrieve(pdb_url, f"{crystal_path}/{pdb_code}.pdb")
    print(f"File {pdb_code}.pdb scaricato con successo!")

def pymol_module(crystal_pdb, save_path):
    # Inizializzazione di PyMOL

    pymol.finish_launching(['pymol', '-qc'])

    cmd.load(crystal_pdb, "cry")
    cmd.remove("cry and !organic")

    # Modifica il ligando: cambio catena e numero residuo
    cmd.alter(f"cry", f"chain = 'A'")
    cmd.alter(f"cry", f"resi = '1'")

    # Applica le modifiche
    cmd.sort()

    ligand_name=f"lig_temp.pdb"
    # Seleziona e salva solo il ligando
    cmd.select("protein", "!organic")
    cmd.remove("protein")
    ligand_name_temp_path=f"{save_path}/{ligand_name}"
    cmd.save(ligand_name_temp_path)

    list_lig=[]
    with open(ligand_name_temp_path, 'r') as f:
        # Leggi tutte le righe del file
        for line in f:
            columns = line.split()
            if len(columns) >= 4:
                lig_code = columns[3]
                list_lig.append(lig_code)
    conteggi = Counter(list_lig)
    ligand_code = conteggi.most_common(1)[0][0]

    cmd.remove(f"!resn {ligand_code}")
    
    cmd.save(ligand_name_temp_path)

    cmd.quit()  

    return ligand_name_temp_path, ligand_code


def main():
    with open("temp.txt", "r") as f:
        lines = [line.strip() for line in f]
        ligand_path=lines[0]
        crystal_path=lines[1]
        crystal_code=lines[2]
        ligand_common_name=lines[3]

    download_crystal_pdb_in_crystal_path(crystal_code, crystal_path)
    results=pymol_module(f"{crystal_path}/{crystal_code}.pdb", ligand_path)
    ligand_temp_path=results[0]
    ligand_code=results[1]

    with open("temp.txt", "w") as f:
        f.write(f"{ligand_path}\n{crystal_path}\n{crystal_code}\n{ligand_common_name}\n{ligand_temp_path}\n{ligand_code}")

if __name__ == "__main__":
    main()

