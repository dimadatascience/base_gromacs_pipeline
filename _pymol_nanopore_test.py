import os
import pymol
from pymol import cmd

location=f"/data/tests/ieo7224_gianluca/GROMACS_workingdirectory"
os.chdir(location) # change the current working directory
workD_sh = os.getcwd() # set the main directory again

lista_alldirs=os.listdir(location)
lista_WD=[]
lista_pdb=[]

for file in lista_alldirs:
    if file.endswith(".pdb"):
        lista_pdb.append(file)
print("PDB files are: ")
index=0
index_list=[]
for i in lista_pdb:
    indexstr=str(index)
    element=str(i)
    ind_elem=f"{indexstr} {element}"
    print (ind_elem)
    index_list.append(lista_pdb.index(i))
    index += 1

stop=False
user_input=input("Choose the index of the crystal pdb file to reference: ")
while stop==False:
    if int(user_input) in index_list:
        print(f"{lista_pdb[int(user_input)]} selected successfully as crystal")
        pymol_crystal=lista_pdb[int(user_input)]
        stop=True
    else:
        user_input=input(f"The index {user_input} is not present in the pdb list, try again\n")

for dirs in lista_alldirs:
    if dirs.startswith('W_'):
        lista_WD.append(dirs)
print("Working directories are: ")
index=0
index_list=[]
for i in lista_WD:
    indexstr=str(index)
    element=str(i)
    ind_elem=f"{indexstr} {element}"
    print (ind_elem)
    index_list.append(lista_WD.index(i))
    index += 1
       
stop=False
user_input=input("Choose the index of the alk protein pdb folder to evaluate: ")
while stop==False:
    if int(user_input) in index_list:
        print(f"{lista_WD[int(user_input)]} selected successfully as protein")
        pymol_protein=f"{lista_WD[int(user_input)]}/start.pdb"
        stop=True
    else:
        user_input=input(f"The index {user_input} is not present in the pdb list, try again\n")

pymol.finish_launching(['pymol', '-qc'])
cmd.load(f"{pymol_protein}", "alk")
cmd.load(f"{pymol_crystal}", "crystal")

cmd.align("crystal", "alk")
cmd.select("crystal and polymer")
cmd.remove("sele")

cmd.select("crystal and organic")
cmd.select("polymer within 4.5 of sele")
list = cmd.index(selection="sele")
with open(f"{lista_WD[int(user_input)]}/index_bindpock.ndx", 'w') as bp:
    bp.write("[ Protein ]\n")
    for i in list:
        bp.write(f"{str(i[1])}\n")

cmd.select("resi 1 and name ca")
list = cmd.index(selection="sele")
with open(f"{lista_WD[int(user_input)]}/index_3Nterm.ndx", 'w') as bp:
    bp.write("[ Protein ]\n")
    for i in list:
        bp.write(f"{str(i[1])}\n")

chainA=cmd.count_atoms("alk & chain A")
with open(f"{lista_WD[int(user_input)]}/chainA.sh", 'w') as bp:
    bp.write(f"#!/bin/sh\n\n# This is in the file heredoc.sh\n\ngmx make_ndx -f md_out.gro -o index_chainA.ndx << EOF\na 1-{chainA}\nq\nEOF")