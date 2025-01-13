import argparse
import subprocess
import sys

def execute_script_in_conda_env(env_name, comando, *args):
    """
    Esegue un comando in un ambiente Conda specificato, mostrando l'output in tempo reale.
    
    :param env_name: Nome dell'ambiente Conda in cui eseguire il comando.
    :param comando: Il comando da eseguire (es. 'python', 'chimerax', 'acpype', ecc.).
    :param args: Argomenti del comando (ad esempio, il percorso dello script).
    """
    try:
        # Prepara il comando completo
        cmd = ["conda", "run", "-n", env_name, comando] + list(args)
        
        # Esegui il comando con Popen per leggere l'output in tempo reale
        process = subprocess.Popen(
            cmd,
            stdout=sys.stdout,  # Invia l'output direttamente al terminale
            stderr=sys.stderr,  # Invia gli errori direttamente al terminale
            text=True
        )
        
        # Aspetta che il processo termini
        process.wait()  # Questo blocca finché il processo non è finito
        
        return process.returncode
        
    except subprocess.CalledProcessError as e:
        print(f"Errore nell'esecuzione del comando: {e}")
        return e.returncode
    except FileNotFoundError:
        print("Errore: il comando 'conda' non è stato trovato.")
        return 1


def main():
    # Creazione dell'oggetto parser per gestire gli argomenti della riga di comando
    parser = argparse.ArgumentParser(description="Copia un file PDB e un ligando in una nuova cartella.")
    
    # Aggiunta dell'argomento -p per il nome del file (senza estensione)
    parser.add_argument('-r', '--reference', type=str, help="Codice del file del cristallo pdb contenente proteina + ligando (es. '6mx8')")
    
    # Aggiunta dell'argomento -l per il ligando o per il nome associato al ligando (es. 'LIGAND' o 'NAME')
    parser.add_argument('-l', '--ligand', type=str, help="nome del ligando in minuscolo da assegnare alla molecola estratta come ligando (es. 'brigatinib')")
    
    # Parsing degli argomenti
    args = parser.parse_args()

    ligand_path="ligands"
    crystal_path="crystal_reference"

    crystal_code=args.reference
    ligand_common_name=args.ligand

    with open("temp.txt", "w") as f:
        f.write(f"{ligand_path}\n{crystal_path}\n{crystal_code}\n{ligand_common_name}")

    print("executing pymol module")
    execute_script_in_conda_env("pygro", "python", "_1_pymol_module_download_and_preparation.py")
    print("pymol module executed")
    print("executing chimerax module")
    execute_script_in_conda_env("chimerax", "chimerax", "--nogui",  "_2_chimeraX_module_addH.py")
    print("chimerax module executed")
    print("executing acpype module")
    execute_script_in_conda_env("acpype", "python", "_3_acpype_module_parametrize_ligand.py")
    print("acpype module executed")
    
if __name__ == "__main__":
    main()