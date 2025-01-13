import argparse 
import numpy as np
import statistics as stt
import re
import os
import pymol
from pymol import cmd
from itertools import combinations

def has_numbers(inputString):
    return any(char.isdigit() for char in inputString)

def rSubset(arr, r):
    # return list of all subsets of length r
    # to deal with duplicate subsets use
    # set(list(combinations(arr, r)))
    return list(combinations(arr, r))

#---------------------------------WORKING_DIRECTORIES_GENERATION/CHOICE---------------------------------

lista_pdb_input=[]

location_GROMACS=(f"{os.getcwd()}/")
lista_allfiles=next(os.walk('.'))[2]
index=0
index_list=[]
print("PDB files are: ")
for i in lista_allfiles:
    if i.endswith(".pdb"):
        lista_pdb_input.append(i)
        indexstr=str(index)
        element=str(i)
        ind_elem=f"{indexstr} {element}"
        print (ind_elem)
        index_list.append(lista_pdb_input.index(i))
        index += 1

workD_final=[]
stop=False
while stop==False:
    workD = input('Choose index of pdb file(s) as input(s) for simulation; ls to list current added pdbs; q to stop adding pds to gromacs simulation pipeline\n')
    if workD=="q":
        stop=True
    elif workD=="ls":
        print(f"The current added folders are: {workD_final}")
    else:
        print(f"{lista_pdb_input[int(workD)]} added to folder list")
        pdb_dir=(lista_pdb_input[int(workD)]).replace(".pdb", "")
        os.system(f"cp -r copy_folder {pdb_dir}")
        os.system(f"mv {lista_pdb_input[int(workD)]} start.pdb")
        os.system(f"mv start.pdb {pdb_dir}/")
        workD_final.append(pdb_dir)


workD_final_indices=[]
for i in workD_final:
    workD_final_indices.append(workD_final.index(i))
# print(workD_final_indices)

#---------------------------------CRYSTAL_CHOICE---------------------------------

lista_pdb=[]

location_temp=(f"{os.getcwd()}/000_referencePDB/")
os.chdir(location_temp) # change the current working directory
workD_sh = os.getcwd() # set the main directory again
lista_allfiles=next(os.walk('.'))[2]
index=0
index_list=[]
print("PDB files are: ")
for i in lista_allfiles:
    if i.endswith(".pdb"):
        lista_pdb.append(i)
        indexstr=str(index)
        element=str(i)
        ind_elem=f"{indexstr} {element}"
        print (ind_elem)
        index_list.append(lista_pdb.index(i))
        index += 1

stop=False
user_input=input("Choose the index of the crystal pdb file to reference for protein alignment and index generation: ")
while stop==False:
    if int(user_input) in index_list:
        print(f"{lista_pdb[int(user_input)]} selected successfully as crystal")
        pymol_crystal=lista_pdb[int(user_input)]
        stop=True
    else:
        user_input=input(f"The index {user_input} is not present in the pdb list, try again\n")

box_command_list=[]
    
for j in workD_final_indices:
    path=f"{location_GROMACS}{workD_final[j]}"
    path_u=f"{location_GROMACS}{workD_final[j]}/"
    os.chdir(path_u) # change the current working directory
    workD_sh = os.getcwd() # set the main directory again
    
    #---------------------------------INDEX_FILES---------------------------------
    
    print(f"{workD_final[j]} selected successfully as protein")
        
    pymol_protein="start.pdb"

    pymol.finish_launching(['pymol', '-qc'])
    cmd.load(f"{pymol_protein}", "alk")
    cmd.load(f"{location_GROMACS}/000_referencePDB/{pymol_crystal}", "crystal")

    cmd.align("crystal", "alk")
    cmd.select("crystal and polymer")
    cmd.remove("sele")

    cmd.select("crystal and organic")
    cmd.select("polymer within 4.5 of sele")
    list = cmd.index(selection="sele")
    with open(f"index_bindpock.ndx", 'w') as bp:
        bp.write("[ Protein ]\n")
        for i in list:
            bp.write(f"{str(i[1])}\n")

    cmd.select("resi 1 and name ca")
    list = cmd.index(selection="sele")
    with open(f"index_3Nterm.ndx", 'w') as bp:
        bp.write("[ Protein ]\n")
        for i in list:
            bp.write(f"{str(i[1])}\n")

    chainA=cmd.count_atoms("alk & chain A")
    with open(f"chainA.sh", 'w') as bp:
        bp.write(f"#!/bin/sh\n\n# This is in the file heredoc.sh\n\ngmx make_ndx -f md_out.gro -o index_chainA.ndx << EOF\na 1-{chainA}\nq\nEOF")

    #---------------------------------BOX_SIZE_CHOICE---------------------------------
    
    box_size=str(input(f"Choose distance between protein {workD_final[j]} and box limit (nm) : "))
    box_command=f"gmx editconf -f start.gro -o newbox.gro -c -d {box_size}.0 -bt dodecahedron"
    box_command_list.append(box_command)

#=============================================================================================================================================================
#====================================================================AUTOMATED_PART===========================================================================
#=============================================================================================================================================================

# for i in workD_final_indices:
#     path=f"{location_GROMACS}{workD_final[i]}"
#     path_u=f"{location_GROMACS}{workD_final[i]}/"
#     os.chdir(path_u) # change the current working directory
#     workD_sh = os.getcwd() # set the main directory again
      
#     #---------------------------------SIMULATION---------------------------------
        
#     if os.path.exists("start.gro"):
#         print(f"file start.gro already exists, remove it or check if it is ok")
#     else: 
#         # echo "1" = use Amber03 forcefield
#         os.system("echo \"1\" | gmx pdb2gmx -f start.pdb -o start.gro -water tip3p")
        
#     if os.path.exists("newbox.gro"):
#         print(f"file newbox.gro already exists, remove it or check if it is ok")
#     else:
#         os.system(box_command_list[i])
        
#     if os.path.exists("solv.gro"):
#         print(f"file solv.gro already exists, remove it or check if it is ok")
#     else: 
#         os.system("gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top")
        
#     if os.path.exists("ions.tpr"):
#         print(f"file ions.tpr already exists, remove it or check if it is ok")
#     else: 
#         os.system("gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr")
        
#     if os.path.exists("solv_ions.gro"):
#         print(f"file solv_ions.gro already exists, remove it or check if it is ok")
#     else: 
#         # echo "13" = substitute solvent
#         os.system("echo \"13\" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral")
        
#     if os.path.exists("em.tpr"):
#         print(f"file em.tpr already exists, remove it or check if it is ok")
#     else: 
#         os.system("gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr")

#     if os.path.exists("em.gro"):
#         print(f"file em.gro already exists, remove it or check if it is ok")
#     else: 
#         os.system("gmx mdrun -v -deffnm em")
        
#     if os.path.exists("nvt.tpr"):
#         print(f"file nvt.tpr already exists, remove it or check if it is ok")
#     else: 
#         os.system("gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr")
        
#     if os.path.exists("nvt.gro"):
#         print(f"file nvt.gro already exists, remove it or check if it is ok")
#     else: 
#         os.system("gmx mdrun -v -deffnm nvt")
        
#     if os.path.exists("npt.tpr"):
#         print(f"file npt.tpr already exists, remove it or check if it is ok")
#     else: 
#         os.system("gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr")
        
#     if os.path.exists("npt.gro"):
#         print(f"file npt.gro already exists, remove it or check if it is ok")
#     else: 
#         os.system("gmx mdrun -v -deffnm npt")
        
#     if os.path.exists("md_out.tpr"):
#         print(f"file md_out.tpr already exists, remove it or check if it is ok")
#     else: 
#         os.system("gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_out.tpr")
        
#     if os.path.exists("md_out.xtc"):
#         print(f"file md_out.xtc already exists, remove it or check if it is ok")
#     else: 
#         os.system("gmx mdrun -v -deffnm md_out -nb gpu")
        
#     if os.path.exists("index_chainA.ndx"):
#         print(f"file index_chainA.ndx already exists, remove it or check if it is ok")
#     else: 
#         # os.system(makendx_command_list[i])
#         os.system("sh chainA.sh ")

#     if os.path.exists("md_out_chainA_centered.xtc"):
#         print(f"file md_out_chainA_centered.xtc already exists, remove it or check if it is ok")
#     else: 
#         os.system("echo \"17\" \"1\" | gmx trjconv -s md_out.tpr  -f md_out.xtc  -o md_out_chainA_centered.xtc -pbc mol -center -n index_chainA.ndx")
        
#     if os.path.exists("md_out_centered.xtc"):
#         print(f"file md_out_centered.xtc already exists, remove it or check if it is ok")
#     else: 
#         os.system("gmx trjconv -f md_out_chainA_centered.xtc  -o md_out_centered.xtc  -pbc nojump")

#     #---------------------------------RMSD_BINDPOCK---------------------------------

#     if os.path.exists("rmsd_bindpock_nofirst10000.xvg"):
#         print(f"file rmsd_bindpock_nofirst10000.xvg already exists, remove it or check if it is ok")
#     else: 
#         os.system("gmx rms -s md_out.tpr -f md_out_centered.xtc -n index_bindpock.ndx -o rmsd_bindpock_nofirst10000.xvg -b 10000")

#     xvg_path_bp=f"{path_u}rmsd_bindpock_nofirst10000.xvg"

#     with open(xvg_path_bp, "r+") as file:
#         xvg = []
#         for line in file.readlines():
#             if "#" not in line:
#                 if '000' in line:
#                     xvg.append(re.sub(r'^.*?000 ', '', line))
            
#         # print(xvg)
#         nums = [float(i.strip()) for i in xvg]
#     # print(sum(nums))
#     # print(len(nums))
#     tot = sum(nums)
#     avg = tot/len(nums)
#     stdev = stt.stdev(nums)

#     print(f"average binding pocket RMSD = {avg}")
#     print(f"average binding pocket RMSD standard deviation = {stdev}")
#     avgRMSD=f"{path_u}avgRMSD_and_stdv.dat"
#     with open(avgRMSD, 'w') as avgRMSD:
#         avgRMSD.write(f"average binding pocket RMSD = {avg} nm\n")
#         avgRMSD.write(f"average binding pocket RMSD standard deviation = {stdev} nm")
        
#     #---------------------------------NTERM_DISTANCE---------------------------------

#     path_index3N=f"{path_u}index_3Nterm.ndx"

#     with open(path_index3N, 'r') as file:
#         index = []
#         for line in file.readlines():
#             if has_numbers(line)==True:
#                 index.append(re.sub(r'[ ,   ]', '', line))
                
#     subset_index=(rSubset(index, 2))      
#     max_value=[]

#     os.system("rm distance_between_3Nterm.dat")
#     for i in subset_index: 
#         atom1=str(int(i[0]))
#         atom2=str(int(i[1]))
#         command_extractRMSD_2atoms=" echo \"atomnr pos1 pos2\" |  gmx distance -f md_out_centered.xtc -s md_out.tpr -n index_3Nterm.ndx -oall pos1-pos2.xvg"
#         command_extractRMSD_2atoms = command_extractRMSD_2atoms.replace("pos1", atom1)
#         command_extractRMSD_2atoms = command_extractRMSD_2atoms.replace("pos2", atom2)
#         filename="pos1-pos2.xvg"
#         filename=filename.replace("pos1", atom1)
#         filename=filename.replace("pos2", atom2)
#         if os.path.exists(filename):
#             print(f"file {filename} already exists")
#         else:
#             os.system(command_extractRMSD_2atoms)
#         with open(filename, 'r') as rmsd_xvg:
#             xvg = []
#             for line in rmsd_xvg.readlines():
#                 if '000' in line:
#                     xvg.append(re.sub(r'^.*?000 ', '', line))
                
#             # print(xvg)
#             nums = [float(i.strip()) for i in xvg]
#             max_value.append(max(nums))
            
#         print(f"the maximum distance between the atom {atom1} and the atom {atom2} is: {max(nums)} nm")
#         with open("distance_between_3Nterm.dat", 'a') as _3Nterm:
#             _3Nterm.write(f"the maximum distance between the atom {atom1} and the atom {atom2} is: {max(nums)} nm\n")

#     print(f"the maximum distance between any Nterm is: {max(max_value)} nm")
#     with open("distance_between_3Nterm.dat", 'a') as _3Nterm:
#         _3Nterm.write(f"the maximum distance between any Nterm is: {max(max_value)} nm")
            
#     #---------------------------------FOLDX_CALCULATIONS---------------------------------

#     t = [5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000]
#     occurrences=["interaction between A and B", "interaction between A and C", "interaction between B and C"]
#     final_avg=[]
#     command_path2 = "/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/foldx_Linux_VM/foldx --command=AnalyseComplex --pdb-dir={path} --pdb=md_out_${}ns.pdb > {path}foldx-${}.out"
#     command_path2 = command_path2.replace("{path}", path_u)

#     for i in t:
#         i2=str(i)
        
#         command1= "echo \"1\" |  gmx trjconv -f md_out_centered.xtc -s md_out.tpr -o md_out_${}ns.pdb -dump ${}"
#         command1 = command1.replace("${}", i2)
#         PDB="md_out_${}ns.pdb"
#         PDB=PDB.replace("${}", i2)
#         if os.path.exists(PDB):
#             print(f"file {PDB} already exists, remove it or check if it is ok")
#         else: 
#             os.system(command1)
        
#         command_path2bis = command_path2.replace("${}", i2)
#         foldxout="foldx-${}.out"
#         foldxout=foldxout.replace("${}", i2)
#         if os.path.exists(foldxout):
#             print(f"file {foldxout} already exists, remove it or check if it is ok")
#         else: 
#             os.system(command_path2bis)
        
#         foldx_path=f"{path_u}{foldxout}"
#         with open(foldx_path, 'r') as txt:
#             lines = txt.readlines()
#             match=[]
#             match_strip=[]
#             for j in occurrences:
#                 condition=False
#                 string="Total          ="
#                 for row in lines:
#                     if j in row:
#                         condition = True
#                     if string in row and condition:
#                         match.append(row.strip()) 
#                         condition=False
#             for k in match:
#                 result=float(re.sub(r'[^0-9,.,-]', '', k))
#                 match_strip.append(result)
#             final_avg.append(sum(match_strip))

#     os.system("rm Indiv*")
#     os.system("rm Interaction*")
#     os.system("rm Interface*")
#     os.system("rm Summary*")

#     tot = sum(final_avg)
#     avg = tot/len(final_avg)
#     avg_str=str(avg)
#     os.system("touch foldx.dat")
#     foldx_dat_path=f"{path_u}foldx.dat"
#     with open(foldx_dat_path, 'w') as dat:
#         dat.write(f"The average complex interaction energy is: {avg_str}")
#         dat.close()
#     print(f"The average complex interaction energy is: {avg_str}")