import os
import shutil
import glob
import re
import string
import pymol
from pymol import cmd


#=============================================FOLDER_PREPARATION=================================================


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

def get_crystal_pdb(dir_crystals, ligand_code):
    # Itera su tutti i file nella directory
    for filename in os.listdir(dir_crystals):
        # Verifica se il file finisce con .pdb
        if filename.endswith('.pdb'):
            # Verifica se il file inizia con il prefisso "ABC"
            if filename.startswith(ligand_code):
                return filename
    return None  # Restituisce None se non viene trovato alcun file



#================================================================================================================
#=============================================FOLDER_PREPARATION=================================================
#================================================================================================================



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


#================================================PYMOL_MODULE====================================================


def get_monomers_from_pdb(pdb_file):
    monomers = []  # Lista per memorizzare le lettere dei monomeri
    # resinumbers = []
    current_monomer = None  # Variabile per tenere traccia dell'attuale monomero
    current_resinumber = 0
    current_atom_num = 0
    continue_counting_atoms = True

    with open(pdb_file, 'r') as f:
        for line in f:
            # Consideriamo solo le righe che contengono atomi
            if line.startswith('ATOM'):
                # Separiamo la riga in base agli spazi
                # print(line)
                columns = line.split()
                # La lettera del monomero è il 5° parametro (indice 4 nella lista)
                atom_num_first_chain = int(columns[1])
                monomer = columns[4]
                resinumber = int(columns [5])
                # Se il monomero è cambiato, aggiungiamolo alla lista
                if monomer != current_monomer:
                    monomers.append(monomer)
                    current_monomer = monomer
                if resinumber >= current_resinumber:
                    # resinumbers.append(resinumber)
                    current_resinumber = resinumber
                if atom_num_first_chain >= current_atom_num and continue_counting_atoms:
                    current_atom_num = atom_num_first_chain
            if line.startswith('TER'):
                continue_counting_atoms = False

    # Restituiamo il numero di monomeri e la lista delle lettere
    return len(monomers), monomers, current_resinumber, current_atom_num

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



#================================================================================================================
#================================================PYMOL_MODULE====================================================
#================================================================================================================



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


#============================================GROMACS_PREPARATION=================================================


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



#================================================================================================================
#============================================GROMACS_PREPARATION=================================================
#================================================================================================================



def gromacs_system_preparation(dir_new, pymol_pdb_file_protein, pymol_pdb_file_ligand, gmx_file_ligand, posre_file_ligand, num_monomers, path_gro_water, last_atom_in_chain_A, ligand_code):
        
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

    complex_atoms_group=f"Protein_{ligand_code}"
    chainA_atoms_group=f"a_1-{last_atom_in_chain_A}"

    with open("chainA.sh", "w") as chainfile:
        chainfile.writelines(
    f"""#!/bin/sh

    # This is in the file heredoc.sh

    gmx make_ndx -f minim2.gro -o index.ndx << EOF
    "Protein" | "{ligand_code}"
    {chainA_atoms_group}
    q
    EOF""")

    os.system(f"bash chainA.sh")

    nsteps = "-1"

    #######         ######
    #########       ######
    ###########     ######
    #############   ######
    ######################
    ######   #############
    ######     ###########
    ######       #########
    ######         #######

    with open("nvt.mdp", "w") as npt:
        npt.writelines(f"""define		= -DPOSRES	; position restrain the protein
    ; Run parameters
    integrator	= md		; leap-frog integrator
    nsteps		= 50000		; 2 * 50000 = 100 ps
    dt		    = 0.002		; 2 fs
    ; Output control
    nstxout		= 500		; save coordinates every 1.0 ps
    nstvout		= 500		; save velocities every 1.0 ps
    nstenergy	= 500		; save energies every 1.0 ps
    nstlog		= 500		; update log file every 1.0 ps
    ; Bond parameters
    continuation	        = no		; first dynamics run
    constraint_algorithm    = lincs	    ; holonomic constraints 
    constraints	            = h-bonds	; all bonds (even heavy atom-H bonds) constrained
    lincs_iter	            = 1		    ; accuracy of LINCS
    lincs_order	            = 4		    ; also related to accuracy
    ; Neighborsearching
    cutoff-scheme   = Verlet
    ns_type		    = grid		; search neighboring grid cells
    nstlist		    = 10		; 20 fs, largely irrelevant with Verlet
    rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)
    rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)
    ; Electrostatics
    coulombtype	    = PME	; Particle Mesh Ewald for long-range electrostatics
    pme_order	    = 4		; cubic interpolation
    fourierspacing	= 0.16	; grid spacing for FFT
    ; Temperature coupling is on
    tcoupl		= V-rescale	            ; modified Berendsen thermostat
    tc-grps		= {complex_atoms_group} Water_and_ions	; two coupling groups - more accurate
    tau_t		= 0.1	  0.1           ; time constant, in ps
    ref_t		= 300 	  300           ; reference temperature, one for each group, in K
    ; Pressure coupling is off
    pcoupl		= no 		; no pressure coupling in NVT
    ; Periodic boundary conditions
    pbc		= xyz		    ; 3-D PBC
    ; Dispersion correction
    DispCorr	= EnerPres	; account for cut-off vdW scheme
    ; Velocity generation
    gen_vel		= yes		; assign velocities from Maxwell distribution
    gen_temp	= 300		; temperature for Maxwell distribution
    gen_seed	= -1		; generate a random seed
    """)
        
    if os.path.exists("nvt.tpr"):
        log_to_file(f"{dir_new}/file.log", f"file nvt.tpr already exists, remove it or check if it is ok")
    else: 
        os.system("gmx grompp -f nvt.mdp -c minim2.gro -r minim2.gro -p topol_lig.top -o nvt.tpr -n index.ndx")
        
    if os.path.exists("nvt.gro"):
        log_to_file(f"{dir_new}/file.log", f"file nvt.gro already exists, remove it or check if it is ok")
    else: 
        os.system("gmx mdrun -v -deffnm nvt")

    #######         ######
    #########       ######
    ###########     ######
    #############   ######
    ######################
    ######   #############
    ######     ###########
    ######       #########
    ######         #######

    with open("npt.mdp", "w") as npt:
        npt.writelines(f"""define		= -DPOSRES	; position restrain the protein
    ; Run parameters
    integrator	= md		; leap-frog integrator
    nsteps		= 50000		; 2 * 50000 = 100 ps
    dt		    = 0.002		; 2 fs
    ; Output control
    nstxout		= 500		; save coordinates every 1.0 ps
    nstvout		= 500		; save velocities every 1.0 ps
    nstenergy	= 500		; save energies every 1.0 ps
    nstlog		= 500		; update log file every 1.0 ps
    ; Bond parameters
    continuation	        = yes		; Restarting after NVT 
    constraint_algorithm    = lincs	    ; holonomic constraints 
    constraints	            = h-bonds	; all bonds (even heavy atom-H bonds) constrained
    lincs_iter	            = 1		    ; accuracy of LINCS
    lincs_order	            = 4		    ; also related to accuracy
    ; Neighborsearching
    cutoff-scheme   = Verlet
    ns_type		    = grid		; search neighboring grid cells
    nstlist		    = 10	    ; 20 fs, largely irrelevant with Verlet scheme
    rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)
    rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)
    ; Electrostatics
    coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics
    pme_order	    = 4		    ; cubic interpolation
    fourierspacing	= 0.16		; grid spacing for FFT
    ; Temperature coupling is on
    tcoupl		= V-rescale	            ; modified Berendsen thermostat
    tc-grps		= {complex_atoms_group} Water_and_ions	; two coupling groups - more accurate
    tau_t		= 0.1	  0.1	        ; time constant, in ps
    ref_t		= 300 	  300	        ; reference temperature, one for each group, in K
    ; Pressure coupling is on
    pcoupl		        = Parrinello-Rahman	    ; Pressure coupling on in NPT
    pcoupltype	        = isotropic	            ; uniform scaling of box vectors
    tau_p		        = 2.0		            ; time constant, in ps
    ref_p		        = 1.0		            ; reference pressure, in bar
    compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1
    refcoord_scaling    = com
    ; Periodic boundary conditions
    pbc		= xyz		; 3-D PBC
    ; Dispersion correction
    DispCorr	= EnerPres	; account for cut-off vdW scheme
    ; Velocity generation
    gen_vel		= no		; Velocity generation is off 
    """)

    if os.path.exists("npt.tpr"):
        log_to_file(f"{dir_new}/file.log", f"file npt.tpr already exists, remove it or check if it is ok")
    else: 
        os.system("gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol_lig.top -o npt.tpr -n index.ndx")
        
    if os.path.exists("npt.gro"):
        log_to_file(f"{dir_new}/file.log", f"file npt.gro already exists, remove it or check if it is ok")
    else: 
        os.system("gmx mdrun -v -deffnm npt")

    ########         #######    ##############
    ##########    ##########    ##################
    ########################    ######       ########
    ###### ########## ######    ######         ######
    ######   ######   ######    ######         ######
    ######     ##     ######    ######         ######
    ######            ######    ######       ########
    ######            ######    ###################
    ######            ######    ################

    with open("md.mdp", "w") as npt:
        npt.writelines(f"""title                   = GROMACS_run 
    ; Run parameters
    integrator              = md        ; leap-frog integrator
    nsteps                  = {nsteps}    ; 2 * 500000 = 1000 ps (1 ns)
    dt                      = 0.002     ; 2 fs
    ; Output control
    nstxout                 = 0         ; suppress bulky .trr file by specifying 
    nstvout                 = 0         ; 0 for output frequency of nstxout,
    nstfout                 = 0         ; nstvout, and nstfout
    nstenergy               = 5000      ; save energies every 10.0 ps
    nstlog                  = 20000      ; update log file every 10.0 ps
    nstxout-compressed      = 5000      ; 1000 = save compressed coordinates every 1.0 ps
    compressed-x-grps       = {complex_atoms_group}    ; What to save
    ; Bond parameters
    continuation            = yes       ; Restarting after NPT 
    constraint_algorithm    = lincs     ; holonomic constraints 
    constraints             = h-bonds   ; bonds involving H are constrained
    lincs_iter              = 1         ; accuracy of LINCS
    lincs_order             = 4         ; also related to accuracy
    ; Neighborsearching
    cutoff-scheme           = Verlet    ; Buffered neighbor searching
    ns_type                 = grid      ; search neighboring grid cells
    nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
    rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
    rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
    ; Electrostatics
    coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
    pme_order               = 4         ; cubic interpolation
    fourierspacing          = 0.16      ; grid spacing for FFT
    ; Temperature coupling is on
    tcoupl                  = V-rescale             ; modified Berendsen thermostat
    tc-grps                 = {complex_atoms_group} Water_and_ions   ; two coupling groups - more accurate
    tau_t                   = 0.1     0.1           ; time constant, in ps
    ref_t                   = 300     300           ; reference temperature, one for each group, in K
    ; Pressure coupling is on
    pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
    pcoupltype              = isotropic             ; uniform scaling of box vectors
    tau_p                   = 2.0                   ; time constant, in ps
    ref_p                   = 1.0                   ; reference pressure, in bar
    compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
    ; Periodic boundary conditions
    pbc                     = xyz       ; 3-D PBC
    ; Dispersion correction
    DispCorr                = EnerPres  ; account for cut-off vdW scheme
    ; Velocity generation
    gen_vel                 = no        ; Velocity generation is off """)
        
    if os.path.exists("md.tpr"):
        log_to_file(f"{dir_new}/file.log", f"file md.tpr already exists, remove it or check if it is ok")
    else: 
        os.system("gmx grompp -f md.mdp -c npt.gro -r npt.gro -p topol_lig.top -o md.tpr -n index.ndx")


