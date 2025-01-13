import pymol
from pymol import cmd
import string
import os
import glob
import sys

def get_next_chain(chain_list):
    # Verifica se la lista è ordinata alfabeticamente
    if chain_list == sorted(chain_list):
        # Trova l'ultima lettera
        ultima_lettera = chain_list[-1]
        
        # Calcola la successiva lettera
        prossimo_elemento = chr(ord(ultima_lettera) + 1)
        
        # Verifica che il prossimo elemento sia una lettera dell'alfabeto
        if prossimo_elemento not in string.ascii_uppercase:
            print("Error: La lista ha raggiunto la fine dell'alfabeto.")
    else:
        print("Error: La lista non è ordinata alfabeticamente.")
    return prossimo_elemento

def get_crystal_pdb(dir_crystals, ligand_code):
    # Itera su tutti i file nella directory
    for filename in os.listdir(dir_crystals):
        # Verifica se il file finisce con .pdb
        if filename.endswith('.pdb'):
            # Verifica se il file inizia con il prefisso "ABC"
            if filename.startswith(ligand_code):
                return filename
    return None  # Restituisce None se non viene trovato alcun file
  
def pymol_module(monomer_list, last_resi, path_pdb_file_protein, path_pdb_file_crystal, path_pdb_file_ligand, dir_new):
    # Inizializzazione di PyMOL
    pymol.finish_launching(['pymol', '-qc'])

    # Carica la proteina multimerica
    cmd.load(path_pdb_file_protein, "prot")

    # Variabili di partenza (ad esempio, lista di monomeri e ligandi)

    # Per ogni monomero della proteina
    for index, n in enumerate(monomer_list):
        # Carica il ligando cristallografico
        # print(index, n, monomer_list, path_pdb_file_crystal)
        cmd.load(path_pdb_file_crystal, f"cry{index}")
        cmd.load(path_pdb_file_ligand, f"lig{index}")

        # Allinea il ligando con il monomero corrispondente
        cmd.align(f"cry{index}", f"prot and chain {n}")
        cmd.align(f"lig{index}", f"cry{index}")

        # Seleziona il ligando e rimuovi tutto ciò che non è ligando
        cmd.remove(f"cry{index}")

    # Rinumerazione e cambiamento di catena del ligando
    resinum = int(last_resi)

    for index, n in enumerate(monomer_list):
        copy = monomer_list
        next_chain = get_next_chain(copy)
        copy.append(next_chain)
        copy.pop(0)

    for index, n in enumerate(copy):
        # Rinumerazione e cambio catena per ciascun ligando
        resinum += 1  # Aumenta il numero del residuo
        
        # Modifica il ligando: cambio catena e numero residuo
        cmd.alter(f"lig{index}", f"chain = '{n}'")
        cmd.alter(f"lig{index}", f"resi = '{resinum}'")

    # Applica le modifiche
    cmd.sort()

    # Salva il complesso finale
    cmd.save(f"{dir_new}/pdb_complex.pdb")

    ligandfile_name=f"pdb_ligand.pdb"
    # Seleziona e salva solo il ligando
    cmd.select("protein", "!organic")
    cmd.remove("protein")
    cmd.save(f"{dir_new}/{ligandfile_name}")

    cmd.delete("all")

    cmd.load(f"{dir_new}/pdb_complex.pdb")

    proteinfile_name=f"pdb_protein.pdb"
    # Seleziona e salva solo la parte proteica
    cmd.select("ligand", "!polymer")
    cmd.remove("ligand")
    cmd.save(f"{dir_new}/{proteinfile_name}")

    cmd.quit()  

    return ligandfile_name, proteinfile_name
    



