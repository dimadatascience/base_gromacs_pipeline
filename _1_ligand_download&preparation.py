import pymol
from pymol import cmd
import os
import argparse
import urllib.request
from collections import Counter

def download_crystal_pdb_in_crystal_path(pdb_code, crystal_path):
    # URL del PDB server
    pdb_url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
    
    # Scarica il file PDB
    urllib.request.urlretrieve(pdb_url, f"{crystal_path}/{pdb_code}.pdb")
    print(f"File {pdb_code}.pdb scaricato con successo!")

def get_ligand_code(ligand_temp_name_path):
    list_lig=[]
    with open(ligand_temp_name_path, 'r') as f:
        # Leggi tutte le righe del file
        for line in f:
            columns = line.split()
            if len(columns) >= 4:
                lig_code = columns[3]
                list_lig.append(lig_code)
    conteggi = Counter(list_lig)
    moda = conteggi.most_common(1)[0][0]

    return moda



def pymol_module(crystal_pdb, save_path):
    # Inizializzazione di PyMOL

    pymol.finish_launching(['pymol', '-qc'])

    cmd.load(crystal_pdb, "cry")
    
    cmd.remove("cry and !organic")

    # Modifica il ligando: cambio catena e numero residuo
    cmd.alter(f"cry", f"chain = 'A'")
    cmd.alter(f"cry", f"resi = '1'")
    cmd.h_add()

    # Applica le modifiche
    cmd.sort()

    ligand_name=f"lig_temp.pdb"
    # Seleziona e salva solo il ligando
    cmd.select("protein", "!organic")
    cmd.remove("protein")
    ligand_name_temp_path=f"{save_path}/{ligand_name}"
    cmd.save(ligand_name_temp_path)

    cmd.quit()  

    return ligand_name_temp_path


def main():
    # Creazione dell'oggetto parser per gestire gli argomenti della riga di comando
    parser = argparse.ArgumentParser(description="Copia un file PDB e un ligando in una nuova cartella.")
    
    # Aggiunta dell'argomento -p per il nome del file (senza estensione)
    parser.add_argument('-r', '--reference', required=True, type=str, help="Codice del file del cristallo pdb contenente proteina + ligando (es. '6mx8.pdb')")
    
    # Aggiunta dell'argomento -l per il ligando o per il nome associato al ligando (es. 'LIGAND' o 'NAME')
    parser.add_argument('-l', '--ligand', type=str, help="nome del ligando in minuscolo (es. 'brigatinib')")
    
    # Parsing degli argomenti
    args = parser.parse_args()

    crystal_code=args.reference
    ligand_common_name=args.ligand
    crystal_path="crystal_reference"
    ligand_path="ligands"

    download_crystal_pdb_in_crystal_path(crystal_code, crystal_path)
    ligand_temp_path=pymol_module(f"{crystal_path}/{crystal_code}.pdb", ligand_path)
    ligand_code=get_ligand_code(ligand_temp_path)

    os.rename(f"{crystal_path}/{crystal_code}.pdb", f"{crystal_path}/{ligand_code}_{crystal_code}.pdb")
    os.rename(ligand_temp_path, f"{ligand_path}/{ligand_code}_{ligand_common_name}.pdb")


if __name__ == "__main__":
    main()

