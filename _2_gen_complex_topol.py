import os
import re
import shutil

def copy_ligand_files(file_input, directory_destinazione):
    # Costruisci il percorso completo del file di destinazione
    destinazione = f"{directory_destinazione}"
    
    # Copia il file alla nuova destinazione
    shutil.copy(file_input, destinazione)
    print(f"File copiato con successo a: {destinazione}")

def extract_ligand_name(file_path):
    # Definisci il pattern regex per trovare il nome della molecola
    # Supponiamo che il nome della molecola sia una stringa alfanumerica con underscore,
    # seguita da uno o più spazi e poi un numero intero (es. 6gy_Briga_H      3)
    pattern = r'\s*([a-zA-Z0-9_]+)\s+\d+\s*$'

    # Apri il file in modalità lettura
    with open(file_path, 'r') as f:
        # Leggi tutte le righe del file
        lines = f.readlines()

    # Flag per sapere quando siamo nella sezione giusta
    dentro_moleculetype = False
    
    # Itera su tutte le righe del file
    for line in lines:
        # Quando troviamo la sezione [ moleculetype ], attiviamo il flag
        if '[ moleculetype ]' in line:
            dentro_moleculetype = True
            continue  # Salta alla riga successiva, che potrebbe contenere il nome della molecola
        
        # Se siamo nella sezione e troviamo un nome molecola con la regex, estraiamo il nome
        if dentro_moleculetype:
            match = re.match(pattern, line.strip())  # Usa strip per rimuovere spazi iniziali e finali
            if match:
                return match.group(1)  # Estrai solo il nome della molecola trovato

    # Se non è stato trovato il nome della molecola
    return None

def inserisci_blocco_testo(file_input, file_output, condition_line, blocco_testo):
    # Leggi il contenuto del file di input
    with open(file_input, 'r') as f:
        lines = f.readlines()
    
    # Trova la posizione per inserire il nuovo blocco di testo
    nuova_lines = []
    for line in lines:
        nuova_lines.append(line)
        
        # Cerca la linea di "Include forcefield parameters" e inserisci il nuovo blocco subito dopo
        if line.strip() == condition_line:
            nuova_lines.append(f"\n{blocco_testo}\n")
    
    # Scrivi il contenuto modificato nel file di output
    with open(file_output, 'w') as f:
        f.writelines(nuova_lines)

def inserisci_blocco_testo_dopo_ultima_occorrenza(file_input, file_output, pattern, blocco_testo):
    with open(file_input, 'r') as f:
        lines = f.readlines()
        testo = ''.join(lines)
        # Trova tutte le occorrenze della regex nel testo
        occorrenze = re.findall(pattern, testo)
        
        # Se ci sono occorrenze, trova l'ultima e inserisci la stringa dopo di essa

        if occorrenze:
            # Trova la posizione dell'ultima occorrenza
            ultima_occorrenza = occorrenze[-1]
            # Trova l'indice dell'ultima occorrenza nel testo
            index = testo.rfind(ultima_occorrenza)
            
            # Aggiungi la stringa dopo l'ultima occorrenza
            nuovo_testo = f"{testo[:index + len(ultima_occorrenza)]}\n{blocco_testo}{testo[index + len(ultima_occorrenza):]}"
    
    with open(file_output, 'w') as f:
        f.writelines(nuovo_testo)       

    
########################

# ligand_id="6GY"
# ligand_itp=f"{ligand_id}_GMX.itp"
# ligand_top=f"{ligand_id}_GMX.top"
# ligand_posre=f"posre_{ligand_id}.itp"

# copy_ligand_files(f"/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/000_referencePDB/ligands/6GY.acpype/{ligand_itp}", "/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/kindomTseq_dimerwt_af0_plus2Brigamol/")
# copy_ligand_files(f"/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/000_referencePDB/ligands/6GY.acpype/{ligand_top}", "/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/kindomTseq_dimerwt_af0_plus2Brigamol/")
# copy_ligand_files(f"/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/000_referencePDB/ligands/6GY.acpype/{ligand_posre}", "/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/kindomTseq_dimerwt_af0_plus2Brigamol/")

# mol_name = extract_ligand_name(f"/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/000_referencePDB/ligands/6GY.acpype/{ligand_itp}")
# if mol_name:
#     print(f'Molecule name: {mol_name}')
# else:
#     print('Molecule name not found.')

# # Specifica il percorso del file di input e di output
# file_input = '/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/kindomTseq_dimerwt_af0_plus2Brigamol/topol_protein.top'
# file_mod = '/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/kindomTseq_dimerwt_af0_plus2Brigamol/topol_protein_mod.top'
# pattern = r'#include "topol_Protein_chain_\w+.itp"'

# # Chiamata alla funzione
# inserisci_blocco_testo(file_input, file_mod, f"#include \"amber03.ff/forcefield.itp\"", f"; Include {ligand_itp} topology\n#include \"{ligand_itp}\"")
# print("lig_itp added successfully")
# inserisci_blocco_testo_dopo_ultima_occorrenza(file_mod, file_mod, pattern, f"\n; Ligand position restraints\n#ifdef POSRES_LIG\n#include \"{ligand_posre}\"\n#endif")
# print("posre_lig_itp added successfully")
# inserisci_blocco_testo_dopo_ultima_occorrenza(file_mod, file_mod, r'Protein_chain_\w+\s+\d+', f"{mol_name}           2")

# # Save in test file
# working_dir = '/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/kindomTseq_dimerwt_af0_plus2Brigamol'  # Cambia con il percorso della tua working directory
# file_top = os.path.join(working_dir, 'topol_test.top')  # Cambia con il percorso del tuo file topol.top

