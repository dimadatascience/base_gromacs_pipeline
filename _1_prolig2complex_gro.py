import re

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
            if re.match(r'^\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+$', line.strip()):
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

# # Esegui il programma
# lig2prot("ligand.gro", "protein.gro", "complex_test.gro")      

