import numpy as np
import pandas as pd 
import argparse
from biopandas.pdb import PandasPdb
from genericpath import isfile
pd.options.mode.chained_assignment = None 

def verbose_print(**kwargs):
    """
    Print message if verbose = True, otherwise do nothing

    Args:
        message (str): Text to print
        verbose (bool): Switch that controls if message is printed
    """
    if kwargs.pop('verbose'):
        print(kwargs.pop('message'))

def create_biopandas_df():
    """
    Creates a BioPandas dataframe using either an input pdb file or downloading a pdf from RCSB

    Returns:
        biopandas_df (cls): BioPandas dataframe with the input pdb data stored
    """
    
    if args.filename:
        if  isfile (args.filename): 
            verbose_print(message=f'Loading {args.filename}', verbose=args.verbose)
            biopandas_df = PandasPdb()                 
            biopandas_df.read_pdb(args.filename)  
            verbose_print(message=f'Creating BioPandas data frame of: {args.filename} ', verbose=args.verbose)
        else:
            raise ValueError (f'{args.pdb_code} does not exist in this folder, STOP \n \
                To download the pdb code from RCSB PDB add --download as argument in the command line ')
    elif args.pdb_code:
        verbose_print (message=f'Downloading {args.pdb_code} from RCSB PDB', verbose=args.verbose)
        biopandas_df = PandasPdb().fetch_pdb(args.pdb_code) 
        verbose_print (message=f'Creating BioPandas data frame of: {args.pdb_code} ', verbose=args.verbose)

    return biopandas_df 

def create_pandas_df_of_the_protein_and_metal_ions (biopandas_df): 

    """
    Creates a data frame that contains only the information of aminoacids and metal ions (if exist) in the protein  
    (this is necessary because biopandas does not include all the functionalities of pandas)
    Args:
        biopandas_df (cls): biopandas data frame with all the pdb information
    Return:
        pdb_df (cls): pandas df with the protein information.
    """

    pdb_df_protein = pd.DataFrame(biopandas_df.df['ATOM']) 

    pdb_df_metal_ions   = pd.DataFrame(biopandas_df.df['HETATM']) 

    pdb_df_metal_ions   = pdb_df_metal_ions.loc[pdb_df_metal_ions['residue_name'] != 'HOH']
    if pdb_df_metal_ions.empty:
        verbose_print (message=f'It does not contain metal ions', verbose=args.verbose)

    pdb_df_metal_ions ['residue_name'] = pdb_df_metal_ions['residue_name'].replace (['CA'],'Ca')
    pdb_df_metal_ions ['atom_name'] = pdb_df_metal_ions['atom_name'].replace (['CA'],'Ca')
    pdb_df = pd.concat ( [pdb_df_protein,pdb_df_metal_ions ], ignore_index = True )

    return pdb_df

def drop_unnecesary_columns (pdb_df):
    """
    Drop columns that are empty or unnecesary in the creation of the coarse grain model
    Args:
        pdb_df (cls): pandas df with the protein information.
    """

    drop_columns = ['record_name','blank_1','blank_2','blank_3','blank_4','alt_loc','insertion','segment_id','charge','line_idx','occupancy','b_factor']
    pdb_df.drop(drop_columns,axis=1,inplace=True)

    return pdb_df 

def change_aminoacids_name_to_one_letter_name (pdb_df): 
    """
    Changes the three letter aminoacid name to the one letter name

    Args:
        pdb_df (cls): pandas df with the protein information.
    """
    amino3to1 = {
        "ASH": "A", "ALA": "A", "CYX": "C", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
        "GLY": "G", "HIS": "H", "HID": "H", "HIE": "H", "HIP": "H", "ILE": "I", "LYS": "K",
        "LEU": "L", "MET": "M", "MSE": "M", "ASN": "N", "PYL": "O", "HYP": "P", "PRO": "P",
        "GLN": "Q", "ARG": "R", "SER": "S", "THR": "T", "SEL": "U", "VAL": "V", "TRP": "W", 
        "TYR": "Y", }

    pdb_df ['resname_one_letter'] = pdb_df ['residue_name']
    
    pdb_df = pdb_df.replace({'resname_one_letter': amino3to1},inplace=True)

    return pdb_df

def create_coarse_grain_df (pdb_df):
    """
    Creates a df with the coarse grain data of the protein

    Args:
        pdb_df (cls): pandas df with the protein information.
    Return:
        coarse_grain_df (df): df with information of the coarse grain model.
    """
    coarse_grain_per_chain = {}  
    coarse_grain_df = pd.DataFrame()

    if args.chain_id is not None:
        verbose_print (message=f'Coarse grain model for chain: {args.chain_id}',verbose=args.verbose)
        pdb_df = pdb_df.loc[pdb_df['chain_id']== args.chain_id]
        alpha_carbon_and_terminus  = create_alpha_carbons_and_terminus_beads (pdb_df)
        sidechain =  create_sidechain_beads  (pdb_df)
        coarse_grain_df= merge_alpha_carbon_and_sidechain (alpha_carbon_and_terminus,sidechain)
        
    else:
        verbose_print (message=f'Coarse grain model for all the sub-units in the PDB', verbose=args.verbose)
        for chain in pdb_df.groupby(['chain_id']).groups.keys():
            pdb_df_per_chain = pdb_df.loc[pdb_df['chain_id']== chain]
            alpha_carbons_and_terminus  = create_alpha_carbons_and_terminus_beads (pdb_df_per_chain)
            sidechain =  create_sidechain_beads  (pdb_df_per_chain)
            coarse_grain_per_chain [chain] = merge_alpha_carbon_and_sidechain (alpha_carbons_and_terminus,sidechain)
            coarse_grain_df = pd.concat ([coarse_grain_df,coarse_grain_per_chain[chain]],ignore_index=0,axis=0)
        
        coarse_grain_df.index = np.arange(1,len(coarse_grain_df)+1)

    return coarse_grain_df

def create_alpha_carbons_and_terminus_beads (pdb_df):

    """
    Creates the alpha carbons beads and terminus beads position from the original protein coordinates.

    Args:
        pdb_df (cls): pandas df with the protein information.
    Return:
        alpha_carbons_and_terminus (dict): dictionary that contains: atom_numbers, position, resname,resname_one_letter,radius,chain
        
    """

    residue_number_list = pdb_df.residue_number.to_list()
    min_residue_number = min (residue_number_list)

    for atom in pdb_df.atom_name.keys():

        if pdb_df.residue_number[atom] == min_residue_number and pdb_df.atom_name[atom] == 'N':

            atom_numbers = pd.Series (pdb_df.residue_number[atom])
            x_coord = pd.Series (pdb_df.x_coord[atom])
            y_coord = pd.Series (pdb_df.y_coord[atom])
            z_coord = pd.Series (pdb_df.z_coord[atom])

            resname = 'n' 
            resname = pd.Series (resname)
            resname_one_letter =  pd.Series (resname)

            resid = pd.Series(pdb_df.residue_number[atom])

            atom_radius = 0.7
            radius = pd.Series(atom_radius)
            chain = pd.Series(pdb_df.chain_id[atom])

        elif pdb_df.atom_name[atom] == 'CA':
        
            atom_number_pdb_ca= pd.Series(pdb_df.residue_number[atom])
            x_coord_ca = pd.Series(pdb_df.x_coord[atom])
            y_coord_ca = pd.Series(pdb_df.y_coord[atom])
            z_coord_ca = pd.Series(pdb_df.z_coord[atom])
            resname_ca = pd.Series(pdb_df.residue_name[atom])
            
            resname_one_letter_ca =  pd.Series ('CA')
            resid_ca = pd.Series(pdb_df.residue_number[atom])

            atom_radius_ca = 0.7
            radius_ca = pd.Series(atom_radius_ca)
            chain_ca = pd.Series(pdb_df.chain_id[atom])

            atom_numbers = pd.concat ([atom_numbers,atom_number_pdb_ca], axis = 0, ignore_index = True)
            x_coord = pd.concat ([x_coord,x_coord_ca], axis=0, ignore_index = True)
            y_coord = pd.concat ([y_coord,y_coord_ca], axis=0, ignore_index = True)
            z_coord = pd.concat ([z_coord,z_coord_ca], axis=0, ignore_index = True)
            resname  = pd.concat ([resname,resname_ca], axis=0, ignore_index=True)
            resname_one_letter  = pd.concat ([resname_one_letter,resname_one_letter_ca], axis=0, ignore_index=True)
            radius = pd.concat ([radius,radius_ca], axis = 0, ignore_index = True )
            chain = pd.concat ([chain,chain_ca],axis =0, ignore_index = 0)

            resid = pd.concat ([resid,resid_ca],axis =0, ignore_index = 0)


        elif pdb_df.residue_number [atom] == pdb_df.residue_number.iloc[-2] and 'C' == pdb_df.atom_name[atom]:

            atom_number_pdb_c= pd.Series(pdb_df.residue_number[atom])
            x_coord_c = pd.Series(pdb_df.x_coord[atom])
            y_coord_c = pd.Series(pdb_df.y_coord[atom])
            z_coord_c = pd.Series(pdb_df.z_coord[atom])

            resname_c = 'c'
            resname_c = pd.Series(resname_c)
            resname_one_letter_c =  pd.Series (resname_c)
            resid_c = pd.Series(pdb_df.residue_number[atom])

            atom_radius_c = 0.7
            radius_c = pd.Series(atom_radius_c)
            chain_c = pd.Series (pdb_df.chain_id[atom])

            atom_numbers = pd.concat([atom_numbers,atom_number_pdb_c], axis=0 , ignore_index = True)
            x_coord = pd.concat([x_coord,x_coord_c], axis=0, ignore_index = True)
            y_coord = pd.concat([y_coord,y_coord_c], axis=0, ignore_index = True)
            z_coord = pd.concat([z_coord,z_coord_c], axis=0, ignore_index = True)
            resname  = pd.concat([resname,resname_c], axis=0, ignore_index = True)
            resname_one_letter  = pd.concat ([resname_one_letter,resname_one_letter_c], axis=0, ignore_index=True)
            radius = pd.concat ([radius,radius_c], axis = 0, ignore_index = True )
            chain = pd.concat ([chain,chain_c],axis =0, ignore_index = 0)
            resid = pd.concat ([resid,resid_c],axis =0, ignore_index = 0)

        elif pdb_df.residue_name[atom] == 'ca' : 

                atom_number_pdb_metal= pd.Series(pdb_df.residue_number[atom])
                x_coord_metal = pd.Series(pdb_df.x_coord[atom])
                y_coord_metal = pd.Series(pdb_df.y_coord[atom])
                z_coord_metal = pd.Series(pdb_df.z_coord[atom])
                
                resname_metal = 'ca'

                resname_metal = pd.Series(resname_metal)
                resname_one_letter_metal =  pd.Series (resname_metal)
                resid_metal = pd.Series(pdb_df.residue_number[atom])

                atom_radius_metal = 2.3
                radius_metal = pd.Series (atom_radius_metal)
                chain_metal = pd.Series(pdb_df.chain_id[atom])

                atom_numbers = pd.concat([atom_numbers,atom_number_pdb_metal], axis=0, ignore_index = True)
                x_coord = pd.concat ([x_coord,x_coord_metal],axis=0, ignore_index = True)
                y_coord = pd.concat ([y_coord,y_coord_metal],axis=0, ignore_index = True)
                z_coord = pd.concat ([z_coord,z_coord_metal],axis=0, ignore_index = True)
                resname  = pd.concat ([resname,resname_metal], axis=0, ignore_index = True)
                resname_one_letter  = pd.concat ([resname_one_letter,resname_one_letter_metal], axis=0, ignore_index=True)

                radius = pd.concat ([radius,radius_metal], axis = 0, ignore_index = True )
                chain = pd.concat ([chain,chain_metal],axis =0, ignore_index = 0)
                resid = pd.concat ([resid,resid_metal],axis =0, ignore_index = 0)

    alpha_carbons_and_terminus = {'atom_numbers': atom_numbers,'x_coord':x_coord,'y_coord':y_coord,'z_coord':z_coord,\
        'resname':resname,'resname_one_letter':resname_one_letter ,'radius': radius, 'chain': chain, 'resid': resid}

    return alpha_carbons_and_terminus

def create_sidechain_beads  (pdb_df) :
        """
        Creates the sidechain beads position, which is created considering the original position of all the atoms in the sidechain 
        and then calculating the center of mass. 

        Args:
            pdb_df (cls): pandas df with the protein information.
        Return:
            residues_bead (dict): dictionary with the information of atom_numbers, coordinates , resname, chain and radius.
        
        """

        pdb_df.loc[pdb_df['element_symbol'] == 'C', 'mass'] = 12.0
        pdb_df.loc[pdb_df['element_symbol'] == 'N', 'mass'] = 14.0
        pdb_df.loc[pdb_df['element_symbol'] == 'O', 'mass'] = 16.0
        pdb_df.loc[pdb_df['element_symbol'] == 'S', 'mass'] = 32.0
        pdb_df.loc[pdb_df['element_symbol'] == 'CA', 'mass'] = 40.0

        x_cm = pdb_df.x_coord*pdb_df.mass
        y_cm = pdb_df.y_coord*pdb_df.mass
        z_cm = pdb_df.z_coord*pdb_df.mass

        pdb_df ['x_cm'] = x_cm
        pdb_df ['y_cm'] = y_cm
        pdb_df ['z_cm'] = z_cm

        pdb_df = pdb_df.loc[(pdb_df.atom_name != 'O') & (pdb_df.atom_name != 'C') & \
                    (pdb_df.atom_name != 'N') & (pdb_df.atom_name != 'CA') & (pdb_df.residue_name != 'GLY') & (pdb_df.atom_name != 'ca')  ]

        pdb_df['sum_x_cm'] = pdb_df.groupby(['residue_number'])['x_cm'].transform(sum)
        pdb_df['sum_y_cm'] = pdb_df.groupby(['residue_number'])['y_cm'].transform(sum)
        pdb_df['sum_z_cm'] = pdb_df.groupby(['residue_number'])['z_cm'].transform(sum)
        pdb_df['sum_mass'] = pdb_df.groupby(['residue_number'])['mass'].transform(sum)
        pdb_df['sum_res'] = pdb_df.groupby('residue_number')['residue_number'].transform('count')

        pdb_df['x_pos'] = pdb_df ['sum_x_cm']/ pdb_df['sum_mass']
        pdb_df['y_pos'] = pdb_df ['sum_y_cm']/ pdb_df['sum_mass']
        pdb_df['z_pos'] = pdb_df ['sum_z_cm']/ pdb_df['sum_mass']

        ''' sidechain radius of gyration  '''
        
        pdb_df['delta_x'] = (pdb_df['x_coord'] - pdb_df['x_pos'])**2
        pdb_df['delta_y'] = (pdb_df['y_coord'] - pdb_df['y_pos'])**2
        pdb_df['delta_z'] = (pdb_df['z_coord'] - pdb_df['z_pos'])**2

        cols_to_sum = ['delta_x','delta_y','delta_z']

        pdb_df['sum_rg'] = pdb_df[cols_to_sum].sum(axis=1)
        pdb_df['sum_rg'] = pdb_df.groupby(['residue_number'])['sum_rg'].transform(sum)
        pdb_df['radius_r'] = np.sqrt(pdb_df['sum_rg']/pdb_df['sum_res'],dtype='float').round(4)

        pdb_df = pdb_df.drop_duplicates(subset='residue_number')

        check_sidechain_radius (pdb_df)
        verbose_print (message='It will be placed a bead with an estimated radius', verbose=args.verbose)

        x_coord_r = pdb_df['x_pos'] 
        y_coord_r = pdb_df['y_pos'] 
        z_coord_r = pdb_df['z_pos'] 
        
        pdb_df['radius_mean'] = pdb_df.groupby(['residue_name'])['radius_r'].transform (np.mean).round(4)

        atom_number_pdb_r = pd.Series(pdb_df.residue_number)
        resname_r = pd.Series(pdb_df.residue_name)
        resname_one_letter_r=  pd.Series (pdb_df.resname_one_letter)
        resid_r = pd.Series(pdb_df.residue_number)

        radius_r = pd.Series(pdb_df['radius_mean'])
        chain_r = pd.Series(pdb_df.chain_id)
                
        residues_bead = {'atom_numbers_r': atom_number_pdb_r,\
             'x_coord_r':x_coord_r,'y_coord_r':y_coord_r,'z_coord_r':z_coord_r, \
                'resname_r':resname_r,'resname_one_letter_r':resname_one_letter_r , 'radius_r':radius_r, 'chain_r': chain_r, 'resid_r': resid_r}

        return  residues_bead

def check_sidechain_radius (pdb_df):
    """
    Checks if the sidechain radius is calculated correctly or if there are missing atoms in the aminoacid.
    Args:
        pdb_df (cls): pandas df with the protein information.
    """
    
    for bead in pdb_df.residue_name.keys():
        if pdb_df.sum_res[bead] == 1:
            if pdb_df.residue_name[bead] == 'ALA':
                pdb_df.radius_r [bead] = 1.0 
            else:
                verbose_print (message=f'{pdb_df.residue_name[bead]} from chain {pdb_df.chain_id[bead]} has missing atoms.', verbose=args.verbose)    
                pdb_df.radius_r [bead] = 1.0    
    return 0

def merge_alpha_carbon_and_sidechain (alpha_carbons_and_terminals,residues_bead):
    """
    Merges the information on alpha_carbons_and_terminals and residues_bead to create the coarse grain df with
    all the data of the protein coarse grain model.

    Args:
        alpha_carbons_and_terminus (dict): dictionary with the information of atom_numbers, coordinates , resname, chain and radius.
        residues_bead (dict): dictionary with the information of atom_numbers, coordinates , resname, chain and radius.

    Return:
        coarse_grain_df (cls): pandas df with the information of the coarse grain beads. 
    """

    if args.model == "2beadAA":

        atom_numbers = pd.concat([alpha_carbons_and_terminals['atom_numbers'],residues_bead ['atom_numbers_r']], axis=0 , ignore_index = True).astype('int')

        x_coord = pd.concat( [alpha_carbons_and_terminals ['x_coord'], residues_bead ['x_coord_r']], axis = 0, ignore_index = True)
        y_coord = pd.concat( [alpha_carbons_and_terminals ['y_coord'], residues_bead ['y_coord_r']], axis = 0, ignore_index = True)
        z_coord = pd.concat( [alpha_carbons_and_terminals ['z_coord'], residues_bead ['z_coord_r']], axis = 0, ignore_index = True)

        resname = pd.concat( [alpha_carbons_and_terminals ['resname'], residues_bead ['resname_r']], axis = 0, ignore_index = True)
        resname_one_letter = pd.concat( [alpha_carbons_and_terminals ['resname_one_letter'],residues_bead ['resname_one_letter_r']], axis = 0, ignore_index = True)

        radius = pd.concat ([alpha_carbons_and_terminals ['radius'], residues_bead ['radius_r']], axis = 0, ignore_index = True)
        chain =  pd.concat ([alpha_carbons_and_terminals ['chain'], residues_bead ['chain_r']], axis = 0, ignore_index = True)

        resid = pd.concat ([alpha_carbons_and_terminals ['resid'], residues_bead ['resid_r']], axis = 0, ignore_index = True)

    elif args.model == "1beadAA":

        atom_numbers = pd.concat([residues_bead['atom_numbers_r']], axis=0 , ignore_index = True).astype('int')

        x_coord = pd.concat( [residues_bead ['x_coord_r']], axis = 0, ignore_index = True)
        y_coord = pd.concat( [residues_bead ['y_coord_r']], axis = 0, ignore_index = True)
        z_coord = pd.concat( [residues_bead ['z_coord_r']], axis = 0, ignore_index = True)

        resname = pd.concat([residues_bead ['resname_r']], axis = 0, ignore_index = True)
        resname_one_letter = pd.concat([residues_bead ['resname_one_letter_r']], axis = 0, ignore_index = True)

        radius = pd.concat ([residues_bead['radius_r']], axis = 0, ignore_index = True)
        chain =  pd.concat ([residues_bead['chain_r']], axis = 0, ignore_index = True)
        resid = pd.concat ([residues_bead ['resid_r']], axis = 0, ignore_index = True)


    coarse_grain_df = pd.DataFrame()
    coarse_grain_df = pd.concat ( [atom_numbers,resname,resname_one_letter, resid, x_coord,y_coord,z_coord,radius,chain], axis = 1, ignore_index = True )
    
    coarse_grain_df.columns = ['atom_numbers','resname','resname_one_letter','resid','x_coord','y_coord','z_coord','radius','chain'] 
    coarse_grain_df.sort_values (by = 'atom_numbers', inplace = True )
    coarse_grain_df.index = np.arange( 1 , len(coarse_grain_df) + 1)

    return coarse_grain_df


def create_list_of_particles_bond (coarse_grain_df):
    """
    Creates a list of all the bonds in the coarse grain model of the protein
    Args:
        coarse_grain_df (cls): pandas df with the information of the coarse grain beads. 
    Return:
        particles_bond_list (list): list with all the particles bond, with the format 'XX:YY' 
    """
    atom_numbers_list = []
    bond_dict = {}
    particles_bond_list = []

    for key in coarse_grain_df.groupby(['chain']).groups.keys():

        coarse_grain_chain = coarse_grain_df.loc[coarse_grain_df['chain']==key]
        particle_id = coarse_grain_chain.loc[(coarse_grain_chain.resname_one_letter == 'CA')]
        particle_id ['indice'] = particle_id.index 
        alpha_carbons_id_list = particle_id['indice'].to_list()

        ''' adds to a list all the alpha carbon bonds '''

        for previous_ca , current_ca in zip (alpha_carbons_id_list, alpha_carbons_id_list[1:]):
            alpha_carbons_bond = f'{previous_ca}:{current_ca}'
            particles_bond_list.append((str (alpha_carbons_bond)))

        ''' adds to a list all the alpha carbons and sidechain bonds '''

        minimun = min (coarse_grain_chain.atom_numbers.keys())
        maximum = max (coarse_grain_chain.atom_numbers.keys())

        for i in  range (minimun, maximum+1): 
            if coarse_grain_chain.atom_numbers[i] not in atom_numbers_list: 
                atom_numbers_list.append (coarse_grain_chain.atom_numbers[i] )

        for residues in atom_numbers_list: 
            bond_dict [residues] = list ()

        for i in  range (minimun,maximum+1): 
            if i < (maximum):
                if coarse_grain_chain.atom_numbers[i] == coarse_grain_chain.atom_numbers[i+1]:

                    bond = f'{i}:{i+1}'
                    bond_dict [coarse_grain_chain.atom_numbers[i]].append(bond)

        for i in bond_dict.keys ():
            if len (bond_dict[i]) > 1:
                particles_bond_list.append(str (bond_dict[i][0]))
                particles_bond_list.append(str(bond_dict[i][1]))
            elif len (bond_dict[i]) == 1:
                particles_bond_list.append(str(bond_dict[i][0]))

    return particles_bond_list

def create_output_coarse_grain_model_as_vtf_file  (coarse_grain,beads_bond, identifier):

    """
    Creates the output as a VTF file that contains the protein coarse grained model information and the bonds
    in the coarse grain model 
    Args:
        coarse_grain (cls): pandas df with the coarse graine information of the protein
        beads_bond (list): list with all the bond between particles 
        identifier (str): protein identifier
    """

    output_file = f'coarse_grain_model_of_{identifier}.vtf'

    verbose_print (message=f'\nCreate output file in VTF format: {output_file} ', verbose=args.verbose)

    coarse_grain = coarse_grain.drop (['atom_numbers'], axis = 1 )

    with open(f'{output_file}','w+') as output_vtf_file:

        for atom in coarse_grain.resname.keys ():
            output_vtf_file.write( f'atom {atom} name {coarse_grain.resname_one_letter[atom]} resname {coarse_grain.resname[atom]} resid {coarse_grain.resid[atom]} chain {coarse_grain.chain[atom]} radius {coarse_grain.radius[atom]} \n')

        for i in (beads_bond):
            output_vtf_file.write(f'bond {i}\n')

        output_vtf_file.write (f'\ntimestep indexed \n')
        cols_to_drop = ['resname','resname_one_letter','radius','chain','resid']
        coarse_grain_df = coarse_grain.drop (cols_to_drop, axis = 1 )
        output_vtf_file.write(coarse_grain_df.to_string (header = False))
        
    output_vtf_file.close()

    return 0

if __name__ == '__main__':
    
    """Creates a coarse-grained model from a pdb structure of a protein. 
    The input can either be a pdb file or a pdb code. 
    If a pdb code is given, the structure of the corresponding protein is downloaded from RCSB.
    By default, it constructs a 2-bead model of the protein with a bead in the alpha carbon and another in the side chain of each aminoacid

    Args:
        --filaname (str,opt): Path to the input pdb file.
        --download_pdb (str,opt): PDB code of the protein. The structure will be directly downloaded from RCSB.
        --model (str): Coarse-grained model to be used. Defaults to 2bead, which considers a bead for the alpha carbon and another for the side chain of each aminoacid
        --chain_id (str): Chain id to be filtered from the pdb structure. Only this selected chain id will be coarse-grained

    Returns:
        coarse_grain_model_of_{pdb_code}.vtf: output file with the coarse-grained structure of the protein in vtf format

    Examples:
        Using  a local pdb file and keeping only chain id A
        `python3 create_coarse_grain_from_pdb.py  --filename tests/1f6s.pdb --chain_id A`
        Downloading a pdb file from RCSB and doing a one bead model (one bead per aminoacid)
        `python3 create_coarse_grain_from_pdb.py  --download_pdb 1f6r --model 1beadAA`
    """

    # Define argparse arguments to the script and parse them
    parser = argparse.ArgumentParser(description='Creates a coarse-grained model from a protein structure given in PDB format')    
    parser.add_argument('--filename', dest='filename', help='\nPath to the PDB file\n')
    parser.add_argument('--download_pdb', dest='pdb_code', help='Downloads the corresponding PDB from RCSB and coarse-grains it') 
    parser.add_argument('--model', dest='model', default='2bead', type=str , help='\nCoarse-grained model to be used\n')
    parser.add_argument('--chain_id', type=str , help='\nSpecific chaid_id to coarse-grain\n') 
    parser.add_argument('--verbose', dest='verbose', action='store_true')
    parser.add_argument('--no-verbose', dest='verbose', action='store_false')
    parser.set_defaults(verbose=True)
    args = parser.parse_args()

    valid_keys_models=["1beadAA","2beadAA"]

    # Check the arguments given

    if args.filename is None and args.pdb_code is None:
        verbose_print(message="WARNING: no inputs were provided, nothing was done", verbose=args.verbose)
        exit()
    elif args.filename is not None and args.pdb_code is not None:
        verbose_print(message="ERROR: --filename and --download_pdb argparse modes cannot be active at the same time, please choose one functionality", verbose=args.verbose)
        exit()
    if args.model not in valid_keys_models:
        verbose_print(message="ERROR: --model only takes " + str(valid_keys_models), verbose=args.verbose)
        exit()
    
    biopandas_df = create_biopandas_df()
    pdb_df = create_pandas_df_of_the_protein_and_metal_ions(biopandas_df)
    drop_unnecesary_columns(pdb_df)
    change_aminoacids_name_to_one_letter_name(pdb_df)
    coarse_grain_df = create_coarse_grain_df(pdb_df)
    beads_bond = create_list_of_particles_bond(coarse_grain_df)
    if args.filename:
        identifier=args.filename[-8:-4]
    elif args.pdb_code:
        identifier=args.pdb_code
    create_output_coarse_grain_model_as_vtf_file(coarse_grain_df,beads_bond,identifier)
    verbose_print(message=f'Finished coarse grain model', verbose=args.verbose) 
