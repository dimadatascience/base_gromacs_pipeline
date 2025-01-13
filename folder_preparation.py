from __python_scripts._functions_folder_preparation import *
import os
import shutil
import argparse
import glob

def main():
    # Creazione dell'oggetto parser per gestire gli argomenti della riga di comando
    parser = argparse.ArgumentParser(description="Copia un file PDB e un ligando in una nuova cartella.")
    
    # Aggiunta dell'argomento -p per il nome del file (senza estensione)
    parser.add_argument('-p', '--protein', required=True, type=str, help="Nome del file pdb (es. 'protein.pdb')")
    
    # Aggiunta dell'argomento -l per il ligando o per il nome associato al ligando (es. 'LIGAND' o 'NAME')
    parser.add_argument('-l', '--ligand', required=True, type=str, help="nome del ligando in minuscolo (es. 'brigatinib') o codice pdb (es. '6GY')")
    
    # Parsing degli argomenti
    args = parser.parse_args()

    test=get_ligand_info(args.ligand)
    ligand_code=test[0]
    ligand_name=test[1]
    pdb_file_ligand=test[2]
    pdb_file_protein=args.protein

    dir_acpype_ligand = f"/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/000_reference/ligands/{ligand_code}.acpype"
    dir_mdp=f"/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/000_reference"
    dir_crystals=f"/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/000_reference/crystal_reference_ccp4_method"

    path_pdb_file_ligand=f"/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/000_reference/ligands/{pdb_file_ligand}"
    path_pdb_file_protein=f"/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/000_reference/af_proteins/{pdb_file_protein}"
    path_gro_water=f"/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/000_reference/tip3p.gro"
    path_pdb_file_crystal=f"{dir_crystals}/{str(get_crystal_pdb(dir_crystals, ligand_code))}"

    dir_new = f"/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/{(args.protein).split('.pdb')[0]}_{ligand_name}" 
    dir_new_folder_log=f"/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/000_reference/dir_new_folder"

    with open(f"{dir_new_folder_log}/{(args.protein).split('.pdb')[0]}_{ligand_name}", 'w') as log:
        log.writelines(dir_new)

#============================================FOLDER_PREPARATION=================================================

    result_folder_preparation=folder_preparation(dir_new, dir_acpype_ligand, ligand_code, ligand_name, dir_mdp)
    gmx_file_ligand=f"{os.path.basename(result_folder_preparation[1])}"
    posre_file_ligand=f"{os.path.basename(result_folder_preparation[2])}"
    
#============================================MONOMER_COUNT=================================================

    log_to_file(f"{dir_new}/file.log", f"\nValutazione del numero di monomeri in corso:")
    num_monomers, monomer_list, last_resi, last_atom_in_chain_A = get_monomers_from_pdb(path_pdb_file_protein)
    log_to_file(f"{dir_new}/file.log", f"Numero di monomeri individuati: {num_monomers}")
    log_to_file(f"{dir_new}/file.log", f"Lista dei monomeri individuati: {', '.join(monomer_list)}")
    log_to_file(f"{dir_new}/file.log", f"Ultimo residuo della catena proteica: {last_resi}")
    
#============================================LIGAND_GENERATION=================================================

    log_to_file(f"{dir_new}/file.log", f"\nAcquisizione directory del cristallo di riferimento in corso:")
    log_to_file(f"{dir_new}/file.log", f"Il cristallo di riferimento per il ligando {ligand_name} specificato è {os.path.basename(path_pdb_file_crystal)}")
    log_to_file(f"{dir_new}/file.log", f"Generazione complesso proteina-ligando in corso:")
    result=pymol_module(monomer_list, last_resi, path_pdb_file_protein, path_pdb_file_crystal, path_pdb_file_ligand, dir_new)
    pymol_pdb_file_ligand=result[0]
    pymol_pdb_file_protein=result[1]
    log_to_file(f"{dir_new}/file.log", f"Generazione complesso proteina-ligando completata")

#============================================GROMACS_PREPARATION=================================================

    gromacs_system_preparation(dir_new, pymol_pdb_file_protein, pymol_pdb_file_ligand, gmx_file_ligand, posre_file_ligand, num_monomers, path_gro_water, last_atom_in_chain_A, ligand_code)
    os.system(f"cd {dir_new}")

if __name__ == "__main__":
    main()



