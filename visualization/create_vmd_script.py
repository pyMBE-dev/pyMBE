import glob
import re
import argparse
import os
from glob import glob
import math as mt

def check_files_extension_consistency (extension_list):
    ''' 
    Check if all the provided file names have the same extension
    
    Args:
        extension_list (list): list that contains all the file extensions
    
    Returns:
        no returns

    Raises:
        ValueError: If not the extensions are not all the same
    
    '''

    for item in extension_list:
        if extension_list[0] != item :
            raise ValueError (f'\n The provided coordinate files have different extensions, found {extension_list[0]} and {item}.\
                \n Please provide file names of one type only')
    return None

def organize_numbered_files (input_file_names):
    ''' 
    Organizes numbered files     
    Inputs:
    - input_file_names (list): list that contains all the files 
    Outputs: 
    - files_ordered (dict): a dictionarty containing the number of the frame and the file information
    - ValueError: if there is one file that does not match the file format. 
  
    '''
    input_file_names = sorted (input_file_names) 
    files_ordered = {}

    for frame in input_file_names:
        num = re.findall(r'\d+', frame)
        if num == []:
            raise ValueError("\nfilename " + frame + " not compatible with the expected filename format filename_{number}.{extension}\n")
        files_ordered [int(num[0])]=frame

    return files_ordered 

def not_numbered_files (input_file_names):
    ''' 
    Organizes the file to have the same format as the numbered files     
    Inputs:
    - input_file_names (list): list that contains all the files 
    Outputs: 
    - files_ordered (dict): a dictionarty containing the number of the frame and the file information
    '''
    n_file = 0 
    files_ordered = {} 
    for frame in input_file_names:
        n_file += 1
        files_ordered[n_file] = frame 
    return files_ordered 

def read_vtf_files (input_file_names): 
    '''
    Reads the information on the vtf files 

    Inputs:
    - input_file_names (list): list that contains all the VTF files
    - args.numbered_files (args parse boolean): for numbered files expects --numbered_files as a shell argument else
                                            if the argument is not given consideres that the input_file_names contains only one file. 
    Outputs: 
    data_vtf_frames a dictionary that contains:
    - frame_atom (dict): a dictionarty which has as a key the frame number an as value a list with [id_atom, r_atom,name_atom,type_atom] 
    - frame_coord (dict): a dictionarty which has as a key the frame number an as value a a list with the coordinates of each particle 
    - N_frame (int) : total number of frames
    - unitcell (dict): information of the box size 
    - type_list (list): list containing the all the atom types 

    '''
    
    type_list = []
    unitcell = {}
    frame_coord = {}
    frame_atom = {}
    frame_bond ={}   
    N_frame=0

    if args.numbered_files:
        files_ordered = organize_numbered_files (input_file_names) 
    else: 
        
        files_ordered = not_numbered_files (input_file_names)
 
    minimum_index = min(files_ordered.keys())
    maximum_index = max(files_ordered.keys())

    for frame in range (minimum_index, maximum_index+1): 
        
        filename = files_ordered [frame]

        with open(filename) as f:   

            frame_atom_list = []
            frame_coord_list = []
            frame_bond_list=[]

            for line in f:

                line_clean=line.split() 
                
                if line_clean: # Non-empty line
                    
                    header=line_clean[0]

                    if header == "unitcell":

                        box = line_clean[1:4]
                        unitcell[N_frame] = box

                    elif header == "atom":
                        
                        id_atom = line_clean[1]
                        r_atom = line_clean[3]
                        name_atom = line_clean[5]
                        type_atom = line_clean[7]
                        data_atom = [id_atom,r_atom,name_atom,type_atom]
                        frame_atom_list.append(data_atom)

                        if type_atom not in type_list:

                            type_list.append(type_atom)

                    elif header == "bond":

                        bond_info = line_clean[1]
                        frame_bond_list.append(bond_info)

                    elif header.isnumeric():
                        
                        id_part = line_clean[0]
                        coord_part = line_clean[1:]
                        frame_coord_list.append([id_part,coord_part])

        frame_atom [N_frame] = frame_atom_list
        frame_bond [N_frame] = frame_bond_list
        frame_coord [N_frame]= frame_coord_list
        
        N_frame +=1

    data_vtf_frames = {'frame_atom' : frame_atom, 'frame_coord':frame_coord, \
        'N_frame': N_frame,'unitcell':unitcell,'type_list':type_list }

    return data_vtf_frames

def count_particles_types_from_vtf_frames (data_vtf_frames):
    '''
    Count the particles on the vtf files 

    Inputs:
    - data_vtf_frames (dict): dictionary that contains: frame_atom, frame_coord, N_frame, unitcell and type_list
   
    Outputs: 
    particles_count A dictionary that contains: 
    - n_part:  total number of particles  
    - n_type (dict): total number of particles of each type 
    - n_type_frame (dict): as a key the number of frame an as a value total number of particles in each frame 
    - name_type (dict): as a key the atom_type and as a value the name of the particle 
    - radius_type (dict): as a key the atom_type and value the radius of the particle
    '''

    n_type = {}
    n_type_frame = {}
   
    name_types = {} 
    radius_types = {} 


    frame_atom = data_vtf_frames['frame_atom']
    type_list = data_vtf_frames['type_list']

    for type in data_vtf_frames['type_list']:
        n_type[type]=0

    for frame in range(data_vtf_frames['N_frame']):
        
        for typ in type_list:
            n_type_frame[typ]=0

        for atom in frame_atom [frame]:
            atom_type = atom[-1]
            n_type_frame[atom_type]+=1

            if atom_type not in name_types.keys():
                name_types[atom_type]=atom[2]
            
            if atom_type not in radius_types.keys():
                radius_types[atom_type]=atom[1]
            
        for type in type_list:
            if n_type_frame[type] > n_type[type]:
                n_type[type] = n_type_frame[type]
    
    n_part = sum(n_type.values())

    particles_count = {'n_part': n_part, 'n_type':n_type, 'n_type_frame':n_type_frame, 'name_types':name_types,\
         'radius_types':radius_types} 

    return particles_count

def read_xyz_files (input_file_names, box = 20.0):
    '''
    Reads the information on the XYZ files 

    Inputs:
    - input_file_names (list): list that contains all the XYZ files
    - args.numbered_files (args parse boolean): for numbered files expects --numbered_files as a shell argument else
                                            if the argument is not given consideres that the input_file_names contains only one file. 
    Outputs: 
    data_xyz a dictionary that contains:
    - frame_list (list): a list that contains the number of frames
    - coord_list (list): a list with the coordinates of each particle 
    - type_atom (list) : a list with the type of each particle 
    - type_list (list): list containing the all the atom types 
    - unitcell (dict): information of the box size 

    '''

    type_list = []
    frame_list = []
    type_atom = []
    x_coord=[]
    y_coord=[]
    z_coord=[]
    coord_list = []
    comment=[]

    if args.numbered_files:
        files_ordered = organize_numbered_files (input_file_names) 
    else: 
        files_ordered = not_numbered_files (input_file_names)
 
    minimum_index = min(files_ordered.keys())
    maximum_index = max(files_ordered.keys())

    for frame in range (minimum_index, maximum_index+1): 
        
        filename = files_ordered [frame]

        with open (filename,'r') as file: 

            for line in file:

                if len(line) < 10:

                    frame_list.append(int(line))

                elif line[0] != "#":

                    type_atom.append(line.split()[0])
                    x_coord = line.split()[1]
                    y_coord = line.split()[2]
                    z_coord = line.split()[3]

                    coord_list.append([x_coord,y_coord,z_coord])

                else:
                    comment.append(line)


            for type in type_atom:
                if type not in type_list:
                    type_list.append(type)

            unitcell = {0: [box,box,box]} 


            data_xyz = {'frame_list' : frame_list ,'type_atom':type_atom,'coord_list':coord_list,\
                'type_list':type_list,'unitcell' :unitcell } 

            return data_xyz

def count_particles_types_from_xyz_file (data_xyz) :

    '''
    Count the particles on the vtf files 

    Inputs:
    - data_vtf_frames (dict): dictionary that contains: frame_atom, frame_coord, N_frame, unitcell and type_list
   
    Outputs: 
    particles_count A dictionary that contains: 
    - n_part:  total number of particles  
    - n_type (dict): total number of particles of each type 
    - n_type_frame (dict): as a key the number of frame an as a value total number of particles in each frame 
    - name_type (dict): as a key the atom_type and as a value the name of the particle 
    - radius_type (dict): as a key the atom_type and value the radius of the particle
    '''

    n_type = {}
    n_type_frame = {}
    name_types = {} 
    radius_types = {} 

    N_frame = 0 
    num_frame = 0

    frame_list = data_xyz ['frame_list']
    type_atom = data_xyz ['type_atom']
    type_list = data_xyz ['type_list']
    for type in type_list:
        n_type [type] = 0 

    for frame in frame_list:

        frame_line = 0 

        for type in type_list:
            n_type_frame [type] = 0

        for line in range (frame):
        
            frame_line = num_frame*frame + line 
            atom_type = type_atom [frame_line]
            n_type_frame [atom_type] += 1   #### amount of types per frame 

            if atom_type not in name_types.keys():
                name_types[atom_type] = type_atom [frame_line]

            if atom_type not in radius_types.keys():
                radius_types[atom_type] = 1

        for type in type_list:

            if n_type_frame[type] > n_type[type]:

                n_type[type] = n_type_frame[type] ### maximum of each type 

        num_frame = num_frame + 1
        N_frame +=1 

    n_part =sum(n_type.values()) 

    particles_count = {'n_part' : n_part ,'n_type' : n_type ,'n_type_frame' : n_type_frame , \
        'name_types' : name_types,'radius_types':radius_types,'N_frame': N_frame} 

    return particles_count

def reorganize_xyz_data (data_xyz,particles_count) :
    ''' 
    reorganize xyz data in dictionary that separates the information in frames 

    Inputs:
    - data_xyz (dict)
    - particles_count (dict )

    Outputs: 
    data_frames A dictionary that contains: 
    - frame_atom (dict): a dictionarty which has as a key the frame number an as value a list with [id_atom, r_atom,name_atom,type_atom] 
    - frame_coord (dict): a dictionarty which has as a key the frame number an as value a a list with the coordinates of each particle 
    - N_frame (int) : total number of frames
    - unitcell (dict): information of the box size 
    - n_type_frame (dict): as a key the number of frame an as a value total number of particles in each frame 
    '''
    nframe = 0

    unitcell = data_xyz ['unitcell']
    type_list = data_xyz['type_list'] 
    n_type_frame = particles_count ['n_type_frame']
    frame_atom = {}
    frame_coord = {}

    for frame in data_xyz['frame_list']: 

        frame_line = 0 

        id_atom = 0

        for type in type_list:

            for line in range (frame):

                frame_line = nframe*frame + line 

                if type == data_xyz['type_atom'][frame_line]:

                    atom_type = data_xyz['type_atom'] [frame_line]
                    n_type_frame [atom_type] += 1  

                    r_atom = 1
                    id_atom +=1 

                    data_atom = [str(id_atom),r_atom,atom_type, atom_type]
                    coord_part = [str(atom_type),data_xyz['coord_list'][frame_line]]
               
                    if nframe in frame_atom :
                        frame_atom[nframe].append(data_atom)
                        frame_coord[nframe].append(coord_part)
                    else:
                        frame_atom[nframe] = [data_atom]
                        frame_coord[nframe] = [coord_part]


        nframe = nframe + 1


    data_frames = {'frame_atom': frame_atom,'frame_coord':frame_coord, 'n_type_frame':n_type_frame,'unitcell' :unitcell, \
        'N_frame': nframe }

    return data_frames

def create_new_particle_list (type_list, particles_count) :

    '''
    create a new particle id list      
    Input:
    - type_list: list of particles types
    - particle_count 

    Output: 
    new_particles_list dictionary that contains: 
    - ids_type (dict):  a dictionarty with a type as a key and the id as a value 
    - ids_list (list): a list with the id of each particle
    - new_id : the new id 
    '''
  
    new_id = 0
    ids_type={}

    n_type = particles_count['n_type']

    for type in type_list:
        
        ids_list=[]

        for id_at in range(new_id,n_type[type]+new_id):
            ids_list.append(id_at)

        new_id = id_at + 1
        ids_type[type]=ids_list

    new_particles_list = {'ids_type':ids_type, 'ids_list':ids_list, 'new_id':new_id }

    return new_particles_list

def create_output_vtf_trajectory_for_vmd (output_trajectory,data_frames,new_particles_list,particles_count) :
    '''
    Creates the output as a vtf file that is vmd compatible

    Input:
    - data_frames: dictionary with the information of particles, coordinates, types of the input files
    - new_particles_list: dictionary with the information of ids
    - particles_count: a dictionary with the information of n_part, names_types,radius_types
    - output_trajectory: the name of the output file 

    Output: 
    - output_trajectory file 
    
    '''
    reservoir_pos=[-2,-2,-2]        # Coordinates of the particle reservoir
    folded = True                   # Applies periodic boundary conditions

    ids_type = new_particles_list ['ids_type']
    n_part = particles_count['n_part']
    name_types = particles_count ['name_types']
    radius_types = particles_count ['radius_types']

    frame_atom = data_frames['frame_atom']
    frame_coord = data_frames['frame_coord']

    unitcell = data_frames['unitcell']

    print ( f'Output trajectory filename: {output_trajectory}')  

    with open(output_trajectory, "w") as output_vmd:

        unit_cell_frame=unitcell[0]    
        output_vmd.write("unitcell")
    
        for cell_coor in unit_cell_frame:
            output_vmd.write(" "+str(cell_coor))
        output_vmd.write("\n")

        for id_at in range(n_part ):

            output_vmd.write(f'atom {str(id_at)}')
        
            for id_list in ids_type.values():

                if id_at in id_list:

                    key_list = list(ids_type.keys())
                    val_list = list(ids_type.values())
                    
                    position = val_list.index(id_list)
                    typ=key_list[position]

            output_vmd.write(f' radius {radius_types[typ]}')
            output_vmd.write(f' name  {name_types[typ]}')
            output_vmd.write(f' typ {str(typ)} \n')

        # write the hidding particle
        hidden_type='hid'
        hidden_name='hid'
        hidden_id=0
        n_part +=1

        output_vmd.write (f'atom {str(hidden_id)}')
        output_vmd.write (' radius 1')
        output_vmd.write (f' name {hidden_name}')
        output_vmd.write (f' typ {hidden_type} \n')

        # write each frame coordinate

        for frame in range(data_frames['N_frame']):
           
            output_vmd.write('\ntimestep indexed\n')
            
            frame_atom_list  = frame_atom[frame]
            frame_coord_list = frame_coord[frame]

            frame_ids_type = ids_type.copy()
            
            for id_at in range(n_part ):

                if id_at < len(frame_coord_list):
                    
                    at_list=frame_atom_list[id_at]
                    at_type=at_list[-1]
                    at_coord=frame_coord_list[id_at][1]
                    poss_id=frame_ids_type[str(at_type)].copy()
                    new_id=poss_id[0]

                    output_vmd.write(str(new_id)+" ")
                    
                    n_cor=0
                    
                    for cor in at_coord:

                        if folded:
                            size_box=float(unit_cell_frame[n_cor])
                            cor=float(cor)
                            pbc_cor=cor-mt.floor(cor/size_box)*size_box
                    
                            output_vmd.write(str(pbc_cor)+" ")
                    
                            n_cor+=1
                        else:
                    
                            output_vmd.write(str(cor)+" ")
                    
                    output_vmd.write("\n") 
                    poss_id.remove(new_id)
                    frame_ids_type[str(at_type)]=poss_id
                
            for id_list in frame_ids_type.values():

                if id_list:

                    for id_at in id_list:

                        output_vmd.write(str(id_at)+" ")
                
                        for cor in reservoir_pos:

                            output_vmd.write(str(cor)+" ")

                        output_vmd.write("\n")

            # Write the hidding particle

            output_vmd.write(f'{str(hidden_id)} ')

            for cor in reservoir_pos:
                output_vmd.write(f'{str(cor)} ')

    return  0

def create_output_visualization_tcl (output_trajectory,output_visualization, type_list,  color_code ):
    '''
    Creates a good-looking visualization of the output trajectory to run in VMD 
    
    Input:
    - output_trajectory: name of the output vtf file 
    - type_list (list): list containin all the particles types
    - color_code (args.parse): if --color_code with the keyword 'protein' process the information with the specific color choice for proteins/peptides
    - output_visualization : name of the visualization file 
    Output:
    - output_visualization file 
    '''
    print ( f'Output visualization filename: {output_visualization}')  

    width_box_line = 3.0

    with open(output_visualization, "w") as output_file:

        output_file.write('mol delete top\n')
        output_file.write('display depthcue off\n')
        output_file.write(f'mol load vtf {os.path.basename(output_trajectory)}\n')
        output_file.write('mol delrep 0 top\n')
        output_file.write('display resetview\n')
        output_file.write('mol representation CPK 1.000000 0.000000\n')
        output_file.write('mol selelection all\n')
        output_file.write('animate goto 0\n')
        output_file.write('color Display Background white\n')
        output_file.write('axes location off\n')
        output_file.write(f'pbc box_draw -color gray -width {str(width_box_line)}\n')

        color=0

        if color_code == 'protein':

            acidic_charged_groups =['D','E','Y','C','c']
            basic_charged_groups = ['HH','KH','RH','nH']

            for typ in type_list:

                if typ in acidic_charged_groups:
                    
                    var="typ"+typ
                    
                    output_file.write(f'set {var} [atomselect top "name {typ} "]\n')
                    output_file.write(f'${var} set radius 1.0\n')
                    output_file.write(f'mol selection name {typ}\n')
                    output_file.write('mol material Opaque\n')
                    output_file.write('mol color ColorID 1 \n') 
                    output_file.write('mol addrep top\n\n')

                elif typ in basic_charged_groups:

                    var="typ"+typ
                    
                    output_file.write(f'set {var} [atomselect top "name {typ} "]\n')
                    output_file.write(f'${var} set radius 1.0\n')
                    output_file.write(f'mol selection name {typ}\n')
                    output_file.write('mol material Opaque\n')
                    output_file.write('mol color ColorID 0 \n') 
                    output_file.write('mol addrep top\n\n')

                elif typ == 'CA':

                    var="typ"+typ
                    
                    output_file.write(f'set {var} [atomselect top "name {typ} "]\n')
                    output_file.write(f'${var} set radius 1.0\n')
                    output_file.write(f'mol selection name {typ}\n')
                    output_file.write('mol material Opaque\n')
                    output_file.write('mol color ColorID 4 \n') 
                    output_file.write('mol addrep top\n\n')

                elif typ == 'Na':

                    var="typ"+typ
                    
                    output_file.write(f'set {var} [atomselect top "name {typ} "]\n')
                    output_file.write(f'${var} set radius 1.0\n')
                    output_file.write(f'mol selection name {typ}\n')
                    output_file.write('mol material Opaque\n')
                    output_file.write('mol color ColorID 10 \n') 
                    output_file.write('mol addrep top\n\n')

                elif typ == 'Cl':

                    var="typ"+typ
                    
                    output_file.write(f'set {var} [atomselect top "name {typ} "]\n')
                    output_file.write(f'${var} set radius 1.0\n')
                    output_file.write(f'mol selection name {typ}\n')
                    output_file.write('mol material Opaque\n')
                    output_file.write('mol color ColorID 9 \n') 
                    output_file.write('mol addrep top\n\n')

                else:

                    var="typ"+typ
                    
                    output_file.write(f'set {var} [atomselect top "name {typ} "]\n')
                    output_file.write(f'${var} set radius 1.0\n')
                    output_file.write(f'mol selection name {typ}\n')
                    output_file.write('mol material Opaque\n')
                    output_file.write('mol color ColorID 7 \n') 
                    output_file.write('mol addrep top\n\n')

        else:  

            for typ in type_list:

                var="typ"+typ
                
                output_file.write(f'set {var} [atomselect top "name {typ} "]\n')
                output_file.write(f'${var} set radius 1.0\n')
                output_file.write(f'mol selection name {typ}\n')
                output_file.write('mol material Opaque\n')
                output_file.write(f'mol color ColorID {str(color)}\n')
                output_file.write('mol addrep top\n\n')
                color=color+1


        hidden_type='hid'
        var="typ"+hidden_type
        output_file.write(f'set {var} [atomselect top "name {hidden_type} "]\n')
        output_file.write(f'${var} set radius 1\n')
        output_file.write('mol representation CPK 1.000000 16.000000\n')
        output_file.write(f'mol selection name {hidden_type}\n')
        output_file.write('mol material Goodsell\n')
        output_file.write(f'mol color ColorID {str(8)}\n')
        output_file.write('mol addrep top\n\n')

    return 0

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Script to create merge the trajectory frames to a single vtf trajectory for VMD and produce a tcl script for visualization of this trajectory.')
    parser.add_argument('file_names', nargs ='*',help='Expected the input files')
    parser.add_argument('--print_input_file_names', help='print the input file names on the screen while loading', default=False, action='store_true') 
    parser.add_argument('--separator', type = str, help = 'the character that separates the base name and the frame number', default=r'\d+' )
    parser.add_argument('--numbered_files', help = 'work with numbered files', default=False, action= 'store_true')
    parser.add_argument('--color_code', type = str, help = 'for protein color_code add protein as an argument', default='' )

    args = parser.parse_args()

    print_input_file_names = args.print_input_file_names

    input_file_names = []
    count_files = 0
    extension_list = []
    for file_name in args.file_names: 
        input_file_names += glob(file_name) 
        extension = os.path.splitext((file_name))[1]
        extension_list.append(extension)
        if count_files > 1 and not args.numbered_files:
            raise ValueError ('\nSeems like multiple frames have been provided. To work with multiple trajectory files add --numbered_files\n')
    
    check_files_extension_consistency (extension_list) 

    for file_name in input_file_names:        
        extension = os.path.splitext((file_name))[1]
        basename = re.split(f'{args.separator}',(os.path.splitext((file_name))[0].split('/'))[-1])[0]
        
        if print_input_file_names:
            print (f'Format: {extension}    Loading: {args.file_names}')

    output_trajectory  = os.path.join(os.path.dirname(file_name), f'{basename}-frames.vtf' ) 
    output_visualization = os.path.join(os.path.dirname(file_name), 'visualization.tcl')

    # FIXME PK: the code below is somewhat repetitive. Each of these functions should take the extension (file type) as argument and internally decide which sub-function to call
    if extension == '.vtf':
        
        print (f'Input files extension: {extension}')

        data_frames = read_vtf_files(input_file_names)

        particles_count = count_particles_types_from_vtf_frames (data_frames)
        
        type_list = data_frames['type_list']

        new_particles_list = create_new_particle_list (type_list,particles_count)
        

    elif extension == '.xyz':

        print (f'Input files extension: {extension}')
    
        data_xyz = read_xyz_files (input_file_names)
        
        particles_count = count_particles_types_from_xyz_file (data_xyz)
        
        data_frames = reorganize_xyz_data (data_xyz, particles_count)

        type_list = data_xyz['type_list']

        new_particles_list = create_new_particle_list (type_list, particles_count)
        

create_output_vtf_trajectory_for_vmd (output_trajectory, data_frames, new_particles_list, particles_count)

create_output_visualization_tcl (output_trajectory= output_trajectory, output_visualization = output_visualization,
                                 type_list = type_list, color_code = args.color_code )

print ('Finished')
