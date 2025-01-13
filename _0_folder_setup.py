import os
import shutil
import glob

def log_to_file(file_path, log_message):
    with open(file_path, 'a') as log:
        log.write(f"{log_message}\n")

def get_ligand_info(ligand_arg):
    allligand_dir = f"/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/000_reference/ligands"
    files = os.listdir(allligand_dir)
    for file in files:
        if file.endswith('.pdb'):
            # Dividi il nome del file usando l'underscore
            parts = file.split('_')
            # print(parts[0])
            if len(parts) >= 2 and parts[0].startswith(ligand_arg):
                ligand_code=parts[0]
                ligand_name=parts[1][:-len(".pdb")]
                pdb_file_ligand=file
            if len(parts) >= 2 and parts[1].startswith(ligand_arg):
                ligand_code=parts[0]
                ligand_name=parts[1][:-len(".pdb")]
                pdb_file_ligand=file
    # print(ligand_code, ligand_name, pdb_file_ligand)
    return ligand_code, ligand_name, pdb_file_ligand

def folder_preparation(dir_new, dir_acpype_ligand, ligand_code, ligand_name, dir_mdp):        
    # Crea la directory se non esiste
    if not os.path.exists(dir_new):
        os.makedirs(dir_new)
        log_to_file(f"{dir_new}/file.log", f"La directory {dir_new} è stata creata con successo")
    else:
        with open(f"{dir_new}/file.log", 'w') as f:
            f.write('Log initialized successfully\n')

    log_to_file(f"{dir_new}/file.log", f"\nInizio copiatura dei file del ligando {ligand_code}_{ligand_name}")       
    # Verifica se la cartella del ligando esiste
    if not os.path.isdir(dir_acpype_ligand):
        print(f"Errore: la cartella del ligando {dir_acpype_ligand} non esiste.")
        return
    
    # Percorsi dei file da copiare
    ligand_gmx_file = os.path.join(dir_acpype_ligand, f"{ligand_code}_GMX.itp")
    # print(ligand_gmx_file)
    posre_ligand_file = os.path.join(dir_acpype_ligand, f"posre_{ligand_code}.itp")
    
    # Verifica se i file esistono
    if not os.path.isfile(ligand_gmx_file):
        print(f"Errore: il file {ligand_gmx_file} non esiste.")
        return
    if not os.path.isfile(posre_ligand_file):
        print(f"Errore: il file {posre_ligand_file} non esiste.")
        return
    
    # Copia dei file del ligando nella nuova cartella
    shutil.copy(ligand_gmx_file, dir_new)
    shutil.copy(posre_ligand_file, dir_new)
    
    log_to_file(f"{dir_new}/file.log", f"I file del ligando sono stati copiati in {dir_new}/")

    # Copia di tutti i file che terminano con .mdp da mdp_dir a new_dir
    log_to_file(f"{dir_new}/file.log", "\nInizio copiatura dei file .mdp")
    mdp_files = glob.glob(os.path.join(dir_mdp, "*.mdp"))

    # Copia ogni file trovato nella directory di destinazione
    for file in mdp_files:
        shutil.copy(file, dir_new) 
        log_to_file(f"{dir_new}/file.log", f"Il file {os.path.basename(file)} e' stato copiato in {dir_new}/")
    
    return dir_new, ligand_gmx_file, posre_ligand_file