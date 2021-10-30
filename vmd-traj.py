from os import listdir
from os.path import isfile, join
import numpy as np
import math as mt
import re

frame_dir='./frames'            # Path to the directory of the frame trajectories
reservoir_pos=[-2,-2,-2]        # Coordinates of the particle reservoir
folded = True                  # Applies periodic boundary conditions
nametraj = "vmd-trajectory.vtf" # Name of the output vmd-readable trajectory file
nametcl  = "visualization.tcl"  # Name of the script for vmd visualization
radcyl = 3                      # Width of the simulation box

# Create the list of files

files = [name for name in listdir(frame_dir) if isfile(join(frame_dir, name))]

files=sorted(files) 

# Order the trajectory files 

files_ordered={}

for frame in files:
    
    num=re.findall(r'\d+', frame)
    files_ordered[int(num[0])]=frame


# Read file data

files_data=[]

unitcell={}
atom={}
bond={}
coord={}

N_frame=0
type_list=[]

minimum_index=min(files_ordered.keys())
maximum_index=max(files_ordered.keys())

for frame in range(minimum_index,maximum_index+1): # 'file' is a builtin type, 'frame' is a less-ambiguous variable name.
    
    file_name=files_ordered[frame]
    path=frame_dir+"/"+file_name

    print("Loading frame " + path)
    with open(path) as f:
        
        frame_atom_list=[]
        frame_bond_list=[]
        frame_coord_list=[]

        for line in f:

            line_clean=line.split() 
            
            if line_clean: # Non-empty line
                
                header=line_clean[0]

                if (header == "unitcell"):

                    box=line_clean[1:4]
                    unitcell[N_frame]=box

                elif (header == "atom"):
                    
                    id_atom=line_clean[1]
                    r_atom=line_clean[3]
                    name_atom=line_clean[5]
                    type_atom=line_clean[7]
                    data_atom=[id_atom,r_atom,name_atom,type_atom]
                    frame_atom_list.append(data_atom)

                    if type_atom not in type_list:

                        type_list.append(type_atom)

                elif (header == "bond"):

                    bond_info=line_clean[1]
                    frame_bond_list.append(bond_info)

                elif (header.isnumeric()):
                    
                    id_part=line_clean[0]
                    coord_part=line_clean[1:]
                    frame_coord_list.append([id_part,coord_part])

    atom[N_frame]=frame_atom_list
    bond[N_frame]=frame_bond_list
    coord[N_frame]=frame_coord_list
    N_frame+=1

# Count the maximum number of particles of each type

N_type={}

for typ in type_list:

    N_type[typ]=0

N_type_frame={}
radi_types={}
name_types={}

for frame in range(N_frame):

    for typ in type_list:

        N_type_frame[typ]=0
    
    for at in atom[frame]:

        at_type=at[-1]
        N_type_frame[at_type]+=1
        
        if at_type not in radi_types.keys():

            radi_types[at_type]=at[1]
        
        if at_type not in name_types.keys():

            name_types[at_type]=at[2]
        
    for typ in type_list:

        if (N_type_frame[typ] > N_type[typ]):

             N_type[typ]=N_type_frame[typ]

N_part=sum(N_type.values())

# Create a new list of particle ids for each type

id0=0
ids_type={}

for typ in type_list:

    ids_list=[]

    for id_at in range(id0,N_type[typ]+id0):
        
        ids_list.append(id_at)

    id0=id_at+1
    ids_type[typ]=ids_list

#write the vmd-readable trajectory file

with open(nametraj, "w") as f_vmd:

    # write the topology

    unit_cell_frame=unitcell[0]
    f_vmd.write("unitcell")
    for cell_coor in unit_cell_frame:
        f_vmd.write(" "+str(cell_coor))
    f_vmd.write("\n")
    for id_at in range(N_part):
        f_vmd.write("atom " + str(id_at))
       
        for id_list in ids_type.values():

            if id_at in id_list:

                key_list = list(ids_type.keys())
                val_list = list(ids_type.values())
                
                position = val_list.index(id_list)
                typ=key_list[position]

        f_vmd.write(" radius " + radi_types[typ])
        f_vmd.write(" name " + name_types[typ])
        f_vmd.write(" type " + str(typ) +"\n")

    # write each frame coordinate

    for frame in range(N_frame):

        f_vmd.write("\n")
        f_vmd.write("timestep indexed")
        f_vmd.write("\n")

        if (unitcell[frame] != unit_cell_frame):

            f_vmd.write("unitcell")
            
            for cell_coor in unit_cell_frame:
                f_vmd.write(" "+str(cell_coor))
            
            f_vmd.write("\n")

        frame_atom_list=atom[frame]
        frame_coord_list=coord[frame]

        frame_ids_type=ids_type.copy()
        
        for id_at in range(N_part):

            if (id_at < len(frame_coord_list)):

                at_list=frame_atom_list[id_at]
                at_type=int(at_list[-1])
                at_coord=frame_coord_list[id_at][1]
                poss_id=frame_ids_type[str(at_type)].copy()
                new_id=poss_id[0]
                f_vmd.write(str(new_id)+" ")
                
                n_cor=0
                for cor in at_coord:

                    if folded:

                        size_box=float(unit_cell_frame[n_cor])
                        cor=float(cor)
                        pbc_cor=cor-mt.floor(cor/size_box)*size_box
                        f_vmd.write(str(pbc_cor)+" ")
                        n_cor+=1
                    else:

                        f_vmd.write(str(cor)+" ")
                
                f_vmd.write("\n") 
                poss_id.remove(new_id)
                frame_ids_type[str(at_type)]=poss_id
               
        for id_list in frame_ids_type.values():

            if id_list:

                for id_at in id_list:

                    f_vmd.write(str(id_at)+" ")
            
                    for cor in reservoir_pos:

                        f_vmd.write(str(cor)+" ")

                    f_vmd.write("\n")

# Write an script to produce a good-looking visualization

with open(nametcl, "w") as f_visual:

    f_visual.write("mol delete top")
    f_visual.write("\n")
    f_visual.write("mol load vtf " + nametraj)
    f_visual.write("\n")
    f_visual.write("mol delrep 0 top")
    f_visual.write("\n")
    f_visual.write("display resetview")
    f_visual.write("\n")

    f_visual.write("pbc box_draw -color purple -width " + str(radcyl) + "\n")

    color=0

    for typ in type_list:

        var="typ"+typ
        f_visual.write("set " + var + ' [atomselect top "name ' + typ + ' "]')
        f_visual.write("\n")
        f_visual.write( "$"+ var + " set radius " +str(radi_types[typ]))
        f_visual.write("\n")
        f_visual.write("mol representation VDW 1.000000 16.000000")
        f_visual.write("\n")
        f_visual.write("mol selection name " + typ)
        f_visual.write("\n")
        f_visual.write("mol material Opaque")
        f_visual.write("\n")
        f_visual.write("mol color ColorID " + str(color))
        f_visual.write("\n")
        f_visual.write("mol addrep top")
        f_visual.write("\n \n")
        color=color+1

    f_visual.write("animate goto 0")
    f_visual.write("\n")
    f_visual.write("color Display Background white")
    f_visual.write("\n")
    f_visual.write("axes location off")
    f_visual.write("\n")

