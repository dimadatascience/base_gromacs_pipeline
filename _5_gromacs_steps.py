import os
import re
import shutil
from __python_scripts._0_folder_setup import *

def extract_ligand_name(molecule_itp_file_path):
    # Definisci il pattern regex per trovare il nome della molecola
    # Supponiamo che il nome della molecola sia una stringa alfanumerica con underscore,
    # seguita da uno o più spazi e poi un numero intero (es. 6gy_Briga_H      3)
    pattern = r'\s*([a-zA-Z0-9_]+)\s+\d+\s*$'

    # Apri il file in modalità lettura
    with open(molecule_itp_file_path, 'r') as f:
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

def prolig2complex(ligand_path, protein_path, complex_path):
    # Leggi il file del ligando
    with open(ligand_path, "r") as ligand_file:
        ligand_list = []
        found_number = False
        for line in ligand_file:
            # Controlla se la riga contiene un numero che indica il numero di atomi
            if not found_number and re.match(r'^\d+$', line.strip()):
                found_number = True
                print("Number of atoms in ligand file:", line.strip())
                lig_atomnum = int(line.strip())
            
            # Aggiungi alla lista solo se la riga non corrisponde al formato delle coordinate
            if found_number and not re.match(r'^\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+$', line.strip()):
                ligand_list.append(f"\t{line.strip()}\n")
        ligand_list.pop(0)

    # Leggi il file della proteina e trova la posizione della riga di fine (Z)
    with open(protein_path, 'r') as protein_file:
        lines = protein_file.readlines()
        found_number = False
        index_z = None
        index_a = None
        for i, line in enumerate(lines):
            if not found_number and re.match(r'^\d+$', line.strip()):
                found_number = True
                print("Number of atoms in protein file:", line.strip())
                prot_atomnum = int(line.strip())
                index_a = i

            # Identifica la riga finale della proteina (corrisponde al formato delle coordinate)
            if re.match(r'^\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+.*$', line.strip()):
                index_z = i  # Trova l'indice della riga finale
                print(f"Found end of protein section at line {i}: {line.strip()}")
                break

    if index_z is None:
        print(f"Errore: non è stata trovata una riga finale valida nel file {protein_path}.")
        return

    # Aggiorna la prima riga del file complex.gro con il numero totale di atomi
    total_atomnum = lig_atomnum + prot_atomnum
    lines[index_a] = f"{total_atomnum}\n"  # La prima riga del file 'complex_test.gro' contiene il numero totale di atomi
    print("New number of atoms in complex file:", total_atomnum)

    # Inserisci gli atomi del ligando prima della riga trovata (index_z)
    new_lines = lines[:index_z]  # Righe prima della riga Z
    new_lines.extend(ligand_list)  # Aggiungi gli atomi del ligando con tabulazione
    new_lines.extend(lines[index_z:])  # Aggiungi le righe che seguono Z

    # Scrivi il nuovo contenuto nel file complesso
    with open(complex_path, 'w') as complex_file:
        complex_file.writelines(new_lines)
    
    print(f"File successfully modified and saved as '{complex_path}'.")

def gromacs_system_preparation(dir_new, pymol_pdb_file_protein, pymol_pdb_file_ligand, gmx_file_ligand, posre_file_ligand, num_monomers, path_gro_water):
        
    os.chdir(dir_new) # change the current working directory
    log_to_file(f"{dir_new}/file.log", f"\nNew working directory set: {dir_new}")
    log_to_file(f"{dir_new}/file.log", f"Starting simulation preparation:")

    if os.path.exists("start.gro"):
        log_to_file(f"{dir_new}/file.log", f"file start.gro already exists, remove it or check if it is ok")
    else: 
        # echo "1" = use Amber03 forcefield
        os.system(f"echo \"1\" | gmx pdb2gmx -f {pymol_pdb_file_protein} -o start.gro -water tip3p -ignh")

    if os.path.exists("ligand.gro"):
        log_to_file(f"{dir_new}/file.log", f"file ligand.gro already exists, remove it or check if it is ok")
    else:
        os.system(f"gmx editconf -f {pymol_pdb_file_ligand} -o ligand.gro")

    log_to_file(f"{dir_new}/file.log", f"\nGenerating complex.gro file from ligand.gro and start.gro")
    prolig2complex("ligand.gro", "start.gro", "complex.gro")
    log_to_file(f"{dir_new}/file.log", f"Complex.gro file generated successfully")

    if os.path.exists("newbox.gro"):
        log_to_file(f"{dir_new}/file.log", f"file newbox.gro already exists, remove it or check if it is ok")
    else:
        os.system(f"gmx editconf -f complex.gro -o newbox.gro -c -d 1.0 -bt dodecahedron")

    ligand_itp_name=extract_ligand_name(gmx_file_ligand)
    log_to_file(f"{dir_new}/file.log", f"\nThe ligand name extracted from {gmx_file_ligand} is {ligand_itp_name}")

    log_to_file(f"{dir_new}/file.log", f"\nGenerating new topology to include ligand")
    inserisci_blocco_testo("topol.top", "topol_lig.top", f"#include \"amber03.ff/forcefield.itp\"", f"; Include {gmx_file_ligand} topology\n#include \"{gmx_file_ligand}\"")
    log_to_file(f"{dir_new}/file.log", "lig_itp added successfully")
    pattern = r'#include "topol_Protein_chain_\w+.itp"'
    inserisci_blocco_testo_dopo_ultima_occorrenza("topol_lig.top", "topol_lig.top", pattern, f"\n; Ligand position restraints\n#ifdef POSRES_LIG\n#include \"{posre_file_ligand}\"\n#endif")
    log_to_file(f"{dir_new}/file.log", "posre_lig_itp added successfully")
    inserisci_blocco_testo_dopo_ultima_occorrenza("topol_lig.top", "topol_lig.top", r'Protein_chain_\w+\s+\d+', f"{ligand_itp_name}           {num_monomers}")
    log_to_file(f"{dir_new}/file.log", "molecule name and number added successfully")

    log_to_file(f"{dir_new}/file.log", f"\nAdding 200 water molecules for ion substitution")
    if os.path.exists("200water_molecules.gro"):
        log_to_file(f"{dir_new}/file.log", f"file 200water_molecules.gro already exists, remove it or check if it is ok")
    else: 
        os.system(f"gmx insert-molecules -f newbox.gro -ci {path_gro_water} -nmol 200 -o 200water_molecules.gro")
    log_to_file(f"{dir_new}/file.log", f"200 water molecules added successfully")

    log_to_file(f"{dir_new}/file.log", f"\nGenerating new topology to include 200 water molecules")   
    inserisci_blocco_testo_dopo_ultima_occorrenza("topol_lig.top", "topol_lig.top", f"{ligand_itp_name}           {num_monomers}", f"SOL           200")
    log_to_file(f"{dir_new}/file.log", "SOL molecules added successfully")

    if os.path.exists("ions.tpr"):
        log_to_file(f"{dir_new}/file.log", f"file ions.tpr already exists, remove it or check if it is ok")
    else: 
        os.system("gmx grompp -f ions.mdp -c 200water_molecules.gro -p topol_lig.top -o ions.tpr")
        
    if os.path.exists("solv_ions.gro"):
        log_to_file(f"{dir_new}/file.log", f"file solv_ions.gro already exists, remove it or check if it is ok")
    else: 
        # echo "13" = substitute solvent
        os.system("echo \"SOL\" | gmx genion -s ions.tpr -o solv_ions.gro -p topol_lig.top -pname NA -nname CL -neutral")

    if os.path.exists("minim1.tpr"):
        log_to_file(f"{dir_new}/file.log", f"file minim1.tpr already exists, remove it or check if it is ok")
    else: 
        os.system("gmx grompp -f minim.mdp -c solv_ions.gro -p topol_lig.top -o minim1.tpr")

    if os.path.exists("minim1.gro"):
        log_to_file(f"{dir_new}/file.log", f"file minim1.gro already exists, remove it or check if it is ok")
    else: 
        os.system("gmx mdrun -v -deffnm minim1")

    if os.path.exists("solv.gro"):
        log_to_file(f"{dir_new}/file.log", f"file solv.gro already exists, remove it or check if it is ok")
    else: 
        os.system("gmx solvate -cp minim1.gro -cs spc216.gro -o solv.gro -p topol_lig.top")

    if os.path.exists("minim2.tpr"):
        log_to_file(f"{dir_new}/file.log", f"file minim2.tpr already exists, remove it or check if it is ok")
    else: 
        os.system("gmx grompp -f minim.mdp -c solv.gro -p topol_lig.top -o minim2.tpr")

    if os.path.exists("minim2.gro"):
        log_to_file(f"{dir_new}/file.log", f"file minim2.gro already exists, remove it or check if it is ok")
    else: 
        os.system("gmx mdrun -v -deffnm minim2")
        





    # if os.path.exists("nvt.tpr"):
    #     log_to_file(f"{dir_new}/file.log", f"file nvt.tpr already exists, remove it or check if it is ok")
    # else: 
    #     os.system("gmx grompp -f nvt.mdp -c minim2.gro -r minim2.gro -p topol_lig.top -o nvt.tpr")
        
    # if os.path.exists("nvt.gro"):
    #     log_to_file(f"{dir_new}/file.log", f"file nvt.gro already exists, remove it or check if it is ok")
    # else: 
    #     os.system("gmx mdrun -v -deffnm nvt")
        
    # if os.path.exists("npt.tpr"):
    #     log_to_file(f"{dir_new}/file.log", f"file npt.tpr already exists, remove it or check if it is ok")
    # else: 
    #     os.system("gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol_lig.top -o npt.tpr")
        
    # if os.path.exists("npt.gro"):
    #     log_to_file(f"{dir_new}/file.log", f"file npt.gro already exists, remove it or check if it is ok")
    # else: 
    #     os.system("gmx mdrun -v -deffnm npt")
        
    # if os.path.exists("md_out.tpr"):
    #     log_to_file(f"{dir_new}/file.log", f"file md_out.tpr already exists, remove it or check if it is ok")
    # else: 
    #     os.system("gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol_lig.top -o md_out.tpr")


        
    # if os.path.exists("md_out.xtc"):
    #     log_to_file(f"{new_dir}/file.log", f"file md_out.xtc already exists, remove it or check if it is ok")
    # else: 
    #     os.system("gmx mdrun -v -deffnm md_out -nb gpu")
        


        
    # if os.path.exists("index_chainA.ndx"):
    #     log_to_file(f"{new_dir}/file.log", f"file index_chainA.ndx already exists, remove it or check if it is ok")
    # else: 
    #     # os.system(makendx_command_list[i])
    #     os.system("sh chainA.sh ")

    # if os.path.exists("md_out_chainA_centered.xtc"):
    #     log_to_file(f"{new_dir}/file.log", f"file md_out_chainA_centered.xtc already exists, remove it or check if it is ok")
    # else: 
    #     os.system("echo \"17\" \"1\" | gmx trjconv -s md_out.tpr  -f md_out.xtc  -o md_out_chainA_centered.xtc -pbc mol -center -n index_chainA.ndx")
        
    # if os.path.exists("md_out_centered.xtc"):
    #     log_to_file(f"{new_dir}/file.log", f"file md_out_centered.xtc already exists, remove it or check if it is ok")
    # else: 
    #     os.system("gmx trjconv -f md_out_chainA_centered.xtc  -o md_out_centered.xtc  -pbc nojump")
