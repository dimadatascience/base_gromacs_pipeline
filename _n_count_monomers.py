import string

def get_monomers_from_pdb(pdb_file):
    monomers = []  # Lista per memorizzare le lettere dei monomeri
    # resinumbers = []
    current_monomer = None  # Variabile per tenere traccia dell'attuale monomero
    current_resinumber = 0
    
    with open(pdb_file, 'r') as f:
        for line in f:
            # Consideriamo solo le righe che contengono atomi
            if line.startswith('ATOM'):
                # Separiamo la riga in base agli spazi
                # print(line)
                columns = line.split()
                # La lettera del monomero è il 5° parametro (indice 4 nella lista)
                monomer = columns[4]
                resinumber = int(columns [5])
                # Se il monomero è cambiato, aggiungiamolo alla lista
                if monomer != current_monomer:
                    monomers.append(monomer)
                    current_monomer = monomer
                if resinumber >= current_resinumber:
                    # resinumbers.append(resinumber)
                    current_resinumber = resinumber
    
    # Restituiamo il numero di monomeri e la lista delle lettere
    return len(monomers), monomers, current_resinumber

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

# lista=["A", "B", "C", "D"]
# get_next_chain(lista)


