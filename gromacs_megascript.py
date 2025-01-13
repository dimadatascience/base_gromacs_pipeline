import argparse 
import numpy as np
import statistics as stt
import re
import os
from itertools import combinations

def has_numbers(inputString):
    return any(char.isdigit() for char in inputString)

def rSubset(arr, r):
    # return list of all subsets of length r
    # to deal with duplicate subsets use
    # set(list(combinations(arr, r)))
    return list(combinations(arr, r))

location=(f"{os.getcwd()}/")
print("Working Directories are: ")
lista_wd=next(os.walk('.'))[1]
# print(lista_wd)

index=0
for i in lista_wd:
    indexstr=str(index)
    element=str(i)
    ind_elem=f"{indexstr} {element}"
    print (ind_elem)
    index += 1
    
# workD = input('Choose index of working folder: ')
# workD_final = lista_wd[int(workD)]
# path=f"{location}{workD_final}"
# path_u=f"{location}{workD_final}/"

workD_final=[]
stop=False
while stop==False:
    workD = input('Choose index of working folder(s); ls to list current added folders; q to stop adding folders to gromacs simulation pipeline\n')
    if workD=="q":
        stop=True
    elif workD=="ls":
        print(f"The current added folders are: {workD_final}")
    else:
        print(f"{lista_wd[int(workD)]} added to folder list")
        workD_final.append(lista_wd[int(workD)])
        
workD_final_indices=[]
for i in workD_final:
    workD_final_indices.append(workD_final.index(i))

box_command_list=[]
# makendx_command_list=[]

for i in workD_final_indices:
    path=f"{location}{workD_final[i]}"
    path_u=f"{location}{workD_final[i]}/"
    os.chdir(path_u) # change the current working directory
    workD_sh = os.getcwd() # set the main directory again
    
    #---------------------------------INDEX_FILES---------------------------------
    
    path_indexbp=f"{path}/index_bindpock.ndx"
    print(f"NOW CHECKING {path_u}index_bindpock.ndx FILE:")
    print("Open Pymol -> Drag and Drop structure pdb file and crystal 6mx8 pdb file")
    print("Align alk… , 6mx8 -> click 6mx8 A -> generate selection -> polymer")
    print("Remove atoms -> Select ligand as sele")
    print("select pocket, polymer within 4.5 of sele -> iterate pocket, print(index)")
    print("Copy [ Protein ] and output of last command in new index file")
    print("Did you create the index_bindpock.ndx file manually as stated above? This is a summary of what is present in the current index_bindpock.ndx file:\n")
    with open(f"{path_indexbp}", 'r') as indexBP:
        lines = indexBP.readlines()
        index_bp_file=0
        while index_bp_file<5:
            x = lines[index_bp_file].replace("\n", "")
            print(str(x))
            index_bp_file+=1

    user_input=input("[yes] to continue, ctrl+z to stop the script\n")
    while not user_input==("yes"):
        print ("You didn't say yes, try again or interrupt the script with crtl+z")
        user_input=input("[yes] to continue, ctrl+z to stop the script\n")

    path_index3N=f"{path}/index_3Nterm.ndx"
    with open(f"{path_index3N}", 'r') as index3N:
        print(f"NOW CHECKING {path_u}index_3Nterm.ndx FILE:")
        print("Get atoms indexes: select nterm, resi 1 and name ca -> iterate nterm, print(index)")
        print("Do this for every Nterm (normally 3)")
        print("Copy [ Protein ] and output of last command in new index file")
        print("Did you create the index_3Nterm.ndx file manually as stated above? This is a summary of what is present in the current index_bindpock.ndx file:\n")
        with open(f"{path_index3N}", 'r') as index3N:
            lines = index3N.readlines()
            index_3N_file=0
            while index_3N_file<4:
                x = lines[index_3N_file].replace("\n", "")
                print(str(x))
                index_3N_file+=1
                
    user_input=input("[yes] to continue, ctrl+z to stop the script\n")
    while not user_input==("yes"):
        print ("You didn't say yes, try again or interrupt the script with crtl+z")
        user_input=input("[yes] to continue, ctrl+z to stop the script\n")

    box_command="gmx editconf -f start.gro -o newbox.gro -c -d {d}.0 -bt dodecahedron"
    box_size=input(f"Choose distance between protein {workD_final[i]} and box limit (nm) : ")
    box_command=box_command.replace("{d}", str(box_size))
    box_command_list.append(box_command)
    
    # makendx_command="echo \"a 1-{chain}\" && echo \"q\" | gmx make_ndx -f md_out.gro -o index_chainA.ndx"
    # # echo 'a 1-2579'  | gmx make_ndx -f md_out.gro -o index_chainA.ndx
    # user_makendx_for_centering=input(f"Choose final atom of chain A for protein {workD_final[i]}:\n")
    # makendx_command=makendx_command.replace("{chain}", user_makendx_for_centering)
    # makendx_command_list.append(makendx_command)  

# #---------------------------------INDEX_FILES2---------------------------------

# path_indexbp=f"{path}/index_bindpock.ndx"
# print("NOW CHECKING 'index_bindpock.ndx' FILE:")
# print("Open Pymol -> Drag and Drop structure pdb file and crystal 6mx8 pdb file")
# print("Align alk… , 6mx8 -> click 6mx8 A -> generate selection -> polymer")
# print("Remove atoms -> Select ligand as sele")
# print("select pocket, polymer within 4.5 of sele -> iterate pocket, print(index)")
# print("Copy [ Protein ] and output of last command in new index file")
# print("Did you create the index_bindpock.ndx file manually as stated above? This is a summary of what is present in the current index_bindpock.ndx file:")
# with open(f"{path_indexbp}", 'r') as indexBP:
#     lines = indexBP.readlines()
#     index_bp_file=0
#     while index_bp_file<5:
#         x = lines[index_bp_file].replace("\n", "")
#         print(str(x))
#         index_bp_file+=1

# user_input=input("[yes] to continue, ctrl+z to stop the script")
# while not user_input==("yes"):
#     print ("You didn't say yes, try again or interrupt the script with crtl+z")
#     user_input=input("[yes] to continue, ctrl+z to stop the script")

# path_index3N=f"{path}/index_3Nterm.ndx"
# with open(f"{path_index3N}", 'r') as index3N:
#     print("NOW CHECKING 'index_3Nterm.ndx' FILE:")
#     print("Get atoms indexes: select nterm, resi 1 and name ca -> iterate nterm, print(index)")
#     print("Do this for every Nterm (normally 3)")
#     print("Copy [ Protein ] and output of last command in new index file")
#     print("Did you create the index_3Nterm.ndx file manually as stated above? This is a summary of what is present in the current index_bindpock.ndx file:")
#     with open(f"{path_index3N}", 'r') as index3N:
#         lines = index3N.readlines()
#         index_3N_file=0
#         while index_3N_file<3:
#             x = lines[index_3N_file].replace("\n", "")
#             print(str(x))
#             index_3N_file+=1
            
# user_input=input("[yes] to continue, ctrl+z to stop the script")
# while not user_input==("yes"):
#     print ("You didn't say yes, try again or interrupt the script with crtl+z")
#     user_input=input("[yes] to continue, ctrl+z to stop the script")

for i in workD_final_indices:
    path=f"{location}{workD_final[i]}"
    path_u=f"{location}{workD_final[i]}/"
    os.chdir(path_u) # change the current working directory
    workD_sh = os.getcwd() # set the main directory again
      
    #---------------------------------SIMULATION---------------------------------

    if os.path.exists("start.pdb"):
        print(f"file start.pdb already exists, remove it or check if it is ok")
    else: 
        os.system("mv *.pdb start.pdb")
        
    if os.path.exists("start.gro"):
        print(f"file start.gro already exists, remove it or check if it is ok")
    else: 
        # echo "1" = use Amber03 forcefield
        os.system("echo \"1\" | gmx pdb2gmx -f start.pdb -o start.gro -water tip3p")
        
    if os.path.exists("newbox.gro"):
        print(f"file newbox.gro already exists, remove it or check if it is ok")
    else:
        os.system(box_command_list[i])
        
    if os.path.exists("solv.gro"):
        print(f"file solv.gro already exists, remove it or check if it is ok")
    else: 
        os.system("gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top")
        
    if os.path.exists("ions.tpr"):
        print(f"file ions.tpr already exists, remove it or check if it is ok")
    else: 
        os.system("gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr")
        
    if os.path.exists("solv_ions.gro"):
        print(f"file solv_ions.gro already exists, remove it or check if it is ok")
    else: 
        # echo "13" = substitute solvent
        os.system("echo \"13\" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral")
        
    if os.path.exists("em.tpr"):
        print(f"file em.tpr already exists, remove it or check if it is ok")
    else: 
        os.system("gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr")

    if os.path.exists("em.gro"):
        print(f"file em.gro already exists, remove it or check if it is ok")
    else: 
        os.system("gmx mdrun -v -deffnm em")
        
    if os.path.exists("nvt.tpr"):
        print(f"file nvt.tpr already exists, remove it or check if it is ok")
    else: 
        os.system("gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr")
        
    if os.path.exists("nvt.gro"):
        print(f"file nvt.gro already exists, remove it or check if it is ok")
    else: 
        os.system("gmx mdrun -v -deffnm nvt")
        
    if os.path.exists("npt.tpr"):
        print(f"file npt.tpr already exists, remove it or check if it is ok")
    else: 
        os.system("gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr")
        
    if os.path.exists("npt.gro"):
        print(f"file npt.gro already exists, remove it or check if it is ok")
    else: 
        os.system("gmx mdrun -v -deffnm npt")
        
    if os.path.exists("md_out.tpr"):
        print(f"file md_out.tpr already exists, remove it or check if it is ok")
    else: 
        os.system("gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_out.tpr")
        
    if os.path.exists("md_out.xtc"):
        print(f"file md_out.xtc already exists, remove it or check if it is ok")
    else: 
        os.system("gmx mdrun -v -deffnm md_out -nb gpu")
        
    if os.path.exists("index_chainA.ndx"):
        print(f"file index_chainA.ndx already exists, remove it or check if it is ok")
    else: 
        # os.system(makendx_command_list[i])
        os.system("sh chainA.sh ")

    if os.path.exists("md_out_chainA_centered.xtc"):
        print(f"file md_out_chainA_centered.xtc already exists, remove it or check if it is ok")
    else: 
        os.system("echo \"17\" \"1\" | gmx trjconv -s md_out.tpr  -f md_out.xtc  -o md_out_chainA_centered.xtc -pbc mol -center -n index_chainA.ndx")
        
    if os.path.exists("md_out_centered.xtc"):
        print(f"file md_out_centered.xtc already exists, remove it or check if it is ok")
    else: 
        os.system("gmx trjconv -f md_out_chainA_centered.xtc  -o md_out_centered.xtc  -pbc nojump")

    #---------------------------------RMSD_BINDPOCK---------------------------------

    if os.path.exists("rmsd_bindpock_nofirst10000.xvg"):
        print(f"file rmsd_bindpock_nofirst10000.xvg already exists, remove it or check if it is ok")
    else: 
        os.system("gmx rms -s md_out.tpr -f md_out_centered.xtc -n index_bindpock.ndx -o rmsd_bindpock_nofirst10000.xvg -b 10000")

    xvg_path_bp=f"{path_u}rmsd_bindpock_nofirst10000.xvg"

    with open(xvg_path_bp, "r+") as file:
        xvg = []
        for line in file.readlines():
            if "#" not in line:
                if '000' in line:
                    xvg.append(re.sub(r'^.*?000 ', '', line))
            
        # print(xvg)
        nums = [float(i.strip()) for i in xvg]
    # print(sum(nums))
    # print(len(nums))
    tot = sum(nums)
    avg = tot/len(nums)
    stdev = stt.stdev(nums)

    print(f"average binding pocket RMSD = {avg}")
    print(f"average binding pocket RMSD standard deviation = {stdev}")
    avgRMSD=f"{path_u}avgRMSD_and_stdv.dat"
    with open(avgRMSD, 'w') as avgRMSD:
        avgRMSD.write(f"average binding pocket RMSD = {avg} nm\n")
        avgRMSD.write(f"average binding pocket RMSD standard deviation = {stdev} nm")
        
    #---------------------------------NTERM_DISTANCE---------------------------------

    with open(path_index3N, 'r') as file:
        index = []
        for line in file.readlines():
            if has_numbers(line)==True:
                index.append(re.sub(r'[ ,   ]', '', line))
                
    subset_index=(rSubset(index, 2))      
    max_value=[]

    os.system("rm distance_between_3Nterm.dat")
    for i in subset_index: 
        atom1=str(int(i[0]))
        atom2=str(int(i[1]))
        command_extractRMSD_2atoms=" echo \"atomnr pos1 pos2\" |  gmx distance -f md_out_centered.xtc -s md_out.tpr -n index_3Nterm.ndx -oall pos1-pos2.xvg"
        command_extractRMSD_2atoms = command_extractRMSD_2atoms.replace("pos1", atom1)
        command_extractRMSD_2atoms = command_extractRMSD_2atoms.replace("pos2", atom2)
        filename="pos1-pos2.xvg"
        filename=filename.replace("pos1", atom1)
        filename=filename.replace("pos2", atom2)
        if os.path.exists(filename):
            print(f"file {filename} already exists")
        else:
            os.system(command_extractRMSD_2atoms)
        with open(filename, 'r') as rmsd_xvg:
            xvg = []
            for line in rmsd_xvg.readlines():
                if '000' in line:
                    xvg.append(re.sub(r'^.*?000 ', '', line))
                
            # print(xvg)
            nums = [float(i.strip()) for i in xvg]
            max_value.append(max(nums))
            
        print(f"the maximum distance between the atom {atom1} and the atom {atom2} is: {max(nums)} nm")
        with open("distance_between_3Nterm.dat", 'a') as _3Nterm:
            _3Nterm.write(f"the maximum distance between the atom {atom1} and the atom {atom2} is: {max(nums)} nm\n")

    print(f"the maximum distance between any Nterm is: {max(max_value)} nm")
    with open("distance_between_3Nterm.dat", 'a') as _3Nterm:
        _3Nterm.write(f"the maximum distance between any Nterm is: {max(max_value)} nm")
            
    #---------------------------------FOLDX_CALCULATIONS---------------------------------

    t = [5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000]
    occurrences=["interaction between A and B", "interaction between A and C", "interaction between B and C"]
    final_avg=[]
    command_path2 = "/data/tests/ieo7224_gianluca/GROMACS_workingdirectory/foldx_Linux_VM/foldx --command=AnalyseComplex --pdb-dir={path} --pdb=md_out_${}ns.pdb > {path}foldx-${}.out"
    # command_path2 = "~/Downloads/foldxMacC11_0/foldx --command=AnalyseComplex --pdb-dir={path} --pdb=test_pdb${}.pdb > {path}foldx-${}.out"
    command_path2 = command_path2.replace("{path}", path_u)

    for i in t:
        i2=str(i)
        
        command1= "echo \"1\" |  gmx trjconv -f md_out_centered.xtc -s md_out.tpr -o md_out_${}ns.pdb -dump ${}"
        command1 = command1.replace("${}", i2)
        PDB="md_out_${}ns.pdb"
        PDB=PDB.replace("${}", i2)
        if os.path.exists(PDB):
            print(f"file {PDB} already exists, remove it or check if it is ok")
        else: 
            os.system(command1)
        
        command_path2bis = command_path2.replace("${}", i2)
        foldxout="foldx-${}.out"
        foldxout=foldxout.replace("${}", i2)
        if os.path.exists(foldxout):
            print(f"file {foldxout} already exists, remove it or check if it is ok")
        else: 
            os.system(command_path2bis)
        
        foldx_path=f"{path_u}{foldxout}"
        with open(foldx_path, 'r') as txt:
            lines = txt.readlines()
            match=[]
            match_strip=[]
            for j in occurrences:
                condition=False
                string="Total          ="
                for row in lines:
                    if j in row:
                        condition = True
                    if string in row and condition:
                        match.append(row.strip()) 
                        condition=False
            for k in match:
                result=float(re.sub(r'[^0-9,.,-]', '', k))
                match_strip.append(result)
            final_avg.append(sum(match_strip))

    os.system("rm Indiv*")
    os.system("rm Interaction*")
    os.system("rm Interface*")
    os.system("rm Summary*")

    tot = sum(final_avg)
    avg = tot/len(final_avg)
    avg_str=str(avg)
    os.system("touch foldx.dat")
    foldx_dat_path=f"{path_u}foldx.dat"
    with open(foldx_dat_path, 'w') as dat:
        dat.write(f"The average complex interaction energy is: {avg_str}")
        dat.close()
    print(f"The average complex interaction energy is: {avg_str}")
