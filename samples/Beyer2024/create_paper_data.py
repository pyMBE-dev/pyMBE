# Import pyMBE and other libraries
import pyMBE
from lib import analysis
import os
import sys
import numpy as np
import argparse 
import subprocess

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(SEED=42)

valid_fig_labels=["7a", "7b", "7c", "8a", "8b", "9"]
valid_modes=["short-run","long-run", "test"]

parser = argparse.ArgumentParser(description='Script to create the data from Beyer2024')
parser.add_argument('--fig_label', 
                    type=str, 
                    required= True,  
                    help=f'Label of the corresponding figure in Beyer2024, currently supported: {valid_fig_labels}')
parser.add_argument('--mode', 
                    type=str, 
                    default= "long-run",  
                    help='Sets for how long the simulation runs, valid modes are {valid_modes}')
parser.add_argument('--plot', action='store_true', help="Switch to plot the data")
args = parser.parse_args()

# Inputs
fig_label=args.fig_label
mode=args.mode
plot=args.plot
# Sanity checks
if fig_label not in valid_fig_labels:
    raise ValueError(f"The figure label {fig_label} is not supported. Supported figure labels are {valid_fig_labels}")


if mode not in valid_modes:
    raise ValueError(f"Mode {mode} is not currently supported, valid modes are {valid_modes}")

## Peptide plots (Fig. 7)
labels_fig7=["7a", "7b", "7c"]

if fig_label in labels_fig7:
    script_path=pmb.get_resource("samples/Beyer2024/peptide.py")
    if fig_label == "7a":
        sequence="K"*5+"D"*5
    elif fig_label == "7b":
        sequence="E"*5+"H"*5
    elif fig_label == "7c":
        sequence="nDSHAKRHHGYKRKFHEKHHSHRGYc"
    else:
        raise RuntimeError()
    pH_range = np.linspace(2, 12, num=21)
    for pH in pH_range:
        run_command=[sys.executable, script_path, "--sequence", sequence, "--pH", str(pH), "--mode", mode]
        print(subprocess.list2cmdline(run_command))
        subprocess.check_output(run_command)



## Protein plots (Fig. 8)

labels_fig8=["8a", "8b"]

if fig_label in labels_fig8:
    
    script_path=pmb.get_resource("samples/Beyer2024/globular_protein.py")
    pH_range = np.linspace(2, 7, num=11)
    run_command_common=[sys.executable, script_path, "--mode", mode, "--no_verbose"]

    if fig_label == "8a":
        
        protein_pdb = "1f6s"
        path_to_cg = f"parameters/globular_proteins/{protein_pdb}.vtf"
        for pH in pH_range:
            
            run_command=run_command_common + ["--pH", str(pH),"--pdb", protein_pdb, "--path_to_cg", path_to_cg, "--metal_ion_name", "Ca", "--metal_ion_charge", str(2)]
            print(subprocess.list2cmdline(run_command))
            subprocess.check_output(run_command)

    elif fig_label == "8b":
        protein_pdb = "1beb"
        path_to_cg = f"parameters/globular_proteins/{protein_pdb}.vtf"
        for pH in pH_range:
            run_command=run_command_common + ["--pH", str(pH),"--pdb", protein_pdb, "--path_to_cg", path_to_cg]
            print(subprocess.list2cmdline(run_command))
            subprocess.check_output(run_command)
    else:
        raise RuntimeError()


## Weak polyelectrolyte dialysis plot (Fig. 9)
if fig_label == "9" and not plot: 
    script_path=pmb.get_resource("samples/Beyer2024/weak_polyelectrolyte_dialysis.py")
    pH_range = np.linspace(1, 13, num=13)
    c_salt_res = 0.01 * pmb.units.mol/pmb.units.L
    for pH in pH_range:
        run_command=[sys.executable, script_path, "--c_salt_res", str(0.01), "--c_mon_sys", str(0.435), "--pH_res", str(pH), "--pKa_value", str(4.0), "--mode", mode]
        print(subprocess.list2cmdline(run_command))
        subprocess.check_output(run_command)

# Analyze all time series
if fig_label in labels_fig7:
    time_series_folder_path=pmb.get_resource("samples/Beyer2024/time_series/peptides")
    
if fig_label in labels_fig8:    
    time_series_folder_path=pmb.get_resource("samples/Beyer2024/time_series/globular_protein")
    
if fig_label == "9":
    time_series_folder_path=pmb.get_resource("samples/Beyer2024/time_series/grxmc")

data=analysis.analyze_time_series(path_to_datafolder=time_series_folder_path)

# Store mean values and other statistics
data_path=pmb.get_resource("samples/Beyer2024/")+"data"
if not os.path.exists(data_path):
    os.makedirs(data_path)
data.to_csv(f"{data_path}/fig{fig_label}.csv")

if plot:
    # Plot the data
    # Import matplotlib for plotting
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r"\usepackage{mathptmx}")
    plt.rcParams["font.family"] = "serif"
    plt.tight_layout()
    mpl.rc('axes', linewidth=1)
    mpl.rcParams['lines.markersize'] = 5
    mpl.rcParams['lines.linewidth'] = 1.0

    width = 3.033
    height = 0.65*width
    plt.figure(figsize=(width,height), dpi=600)
    plt.grid(which='major', 
             color='#CCCCCC', 
             linestyle='--', 
             linewidth=0.6)

    # Set labels for the axes
    if fig_label in labels_fig7+labels_fig8:
        plt.ylabel(r"Net charge $Z/e$")
        plt.xlabel(r"pH in the solution")
    elif fig_label == "9":
        plt.ylabel(r"Degree of ionization $\alpha$")
        plt.xlabel(r"pH in the reservoir")
    else:
        raise RuntimeError()

    # Load pka set
    if fig_label in ["7a","7b"]:
        pka_path=pmb.get_resource("parameters/pka_sets/CRC1991.json")
        pmb.load_pka_set (filename=pka_path)
    elif fig_label in ["7c", "8a", "8b"]:
        pka_path=pmb.get_resource("parameters/pka_sets/Nozaki1967.json")
        pmb.load_pka_set (filename=pka_path)
        # FIXME: this is only necessary due to an undesired feature in calculate_HH
        # that forces to have all particles defined in pyMBE
        par_path=pmb.get_resource("parameters/peptides/Blanco2021.json")
        pmb.load_interaction_parameters(par_path)

    # Load ref data    
    if fig_label == "7a":
        ref_path=pmb.get_resource("testsuite/data/Lunkad2021a.csv")
    elif fig_label == "7b":
        ref_path=pmb.get_resource("testsuite/data/Lunkad2021b.csv")
    elif fig_label == "7c":
        ref_path=pmb.get_resource("testsuite/data/Blanco2020a.csv")
    elif fig_label == "8a":
        ref_path=pmb.get_resource("testsuite/data/Torres2022.csv")
    elif fig_label == "8b":
        ref_path=pmb.get_resource("testsuite/data/Torres2017.csv")
    elif fig_label == "9":
        ref_path=pmb.get_resource("testsuite/data/Landsgesell2020a.csv")
    else:
        raise RuntimeError()

    ref_data=analysis.read_csv_file(path=ref_path)

    # Calculate and plot Henderson-Hasselbalch (HH)
    if fig_label in labels_fig7:
        
        pmb.define_peptide (name=sequence, 
                            sequence=sequence,
                            model="1beadAA")
        pH_range_HH = np.linspace(2, 12, num=1000)

        Z_HH = pmb.calculate_HH(molecule_name=sequence,
                                pH_list=pH_range_HH)

        # Plot HH
        plt.plot(pH_range_HH,
                   Z_HH, 
                   label=r"HH", 
                   color="black")

    elif fig_label in labels_fig8:
        
        if fig_label == "8a":
            protein_pdb = '1f6s'
        elif fig_label == "8b":
            protein_pdb = '1beb'
    
        path_to_cg=pmb.get_resource(f"parameters/globular_proteins/{protein_pdb}.vtf")
        topology_dict = pmb.read_protein_vtf_in_df (filename=path_to_cg)
    
        pmb.define_protein (name=protein_pdb, 
                            topology_dict=topology_dict, 
                            model = '2beadAA')
            
        pH_range_HH = np.linspace(2, 7, num=1000)
        
        Z_HH = pmb.calculate_HH(molecule_name=protein_pdb,
                                pH_list=pH_range_HH)

        print (Z_HH)
         # Plot HH
        plt.plot(pH_range_HH,
               Z_HH, 
               label=r"HH", 
               color="black")
    

    elif fig_label == "9":
        c_salt_res = 0.01 * pmb.units.mol/pmb.units.L
        pmb.define_particle(name='A', acidity='acidic', sigma=1*pmb.units('reduced_length'), epsilon=1*pmb.units('reduced_energy'), pka=4.0)
        pmb.define_residue(name='rA', central_bead="A", side_chains=[])
        pmb.define_molecule(name='polyacid', residue_list=['rA'])

        pH_range = np.linspace(1.0, 13.0, num=1000)
        Z_HH = pmb.calculate_HH(molecule_name='polyacid', pH_list=pH_range)
        alpha_HH = np.abs(np.asarray(Z_HH))

        HH_Donnan_charge_dict = pmb.calculate_HH_Donnan(
                c_macro={'polyacid': 0.435*pmb.units.mol/pmb.units.L},
                c_salt=c_salt_res, 
                pH_list=pH_range)
        Z_HH_Donnan = HH_Donnan_charge_dict["charges_dict"]
        alpha_HH_Donnan = np.abs(np.asarray(Z_HH_Donnan['polyacid']))

        # Plot HH
        plt.plot(pH_range,
                   alpha_HH, 
                   label=r"HH", 
                   color="black")

        plt.plot(pH_range,
                   alpha_HH_Donnan, 
                   label=r"HH+Don", 
                   color="#666666", 
                   linestyle="dashed")

    # Plot Ref data
    ref_data=ref_data.sort_values(by="pH",ascending=True)
    
    if fig_label in ["7a","7b"]:
        style={"linestyle":"none", 
                "marker":"s", 
                "label":"Lunkad  et al.", 
                "color":"C0"}
    elif fig_label == "7c":
        style={"linestyle":"none", 
                "marker":"^", 
                "label":"Blanco  et al.", 
                "color":"green",  
                "markeredgewidth":1.5}
        
    elif fig_label in labels_fig8:
        
        style={"linestyle":"none", 
            "marker":"s", 
            "label":"Torres  et al.", 
            "color":"C0",
            "markeredgewidth":1.5}
    
        
    elif fig_label == "9":
        style={"linestyle":"none", 
                "marker":"s", 
                "label":"Landsgesell  et al.", 
                "color":"C0",  
                "markeredgewidth":1.5}

    if fig_label in labels_fig7+labels_fig8:
        plt.errorbar(ref_data["pH"], 
                       ref_data["charge"], 
                       ref_data["charge_error"], 
                        **style)

    elif fig_label == "9":
        ref_data = ref_data.loc[np.isclose(ref_data["cs_bulk"], (c_salt_res.to('mol/reduced_length**3')*pmb.N_A).magnitude, rtol=1e-03) & np.isclose(ref_data["c_poly_[mol/l]"], 0.435/50, rtol=1e-03) & np.isclose(ref_data["Kcideal_in_mol_per_l"], 10**(-4), rtol=1e-03)]
        plt.errorbar(ref_data["pH"], 
                        ref_data["degree_of_dissociation"], 
                        yerr=ref_data["err_degree_of_dissociation"], 
                        **style)


    # Plot data produced by pyMBE
    if fig_label in ["7a", "7b", "7c"]:
        data=data.astype({("pH","value"): 'float'}).sort_values(by=("pH","value"))
        data=data[data.sequence.value == f"{sequence}"]
        plt.errorbar(data["pH"]["value"], 
                       data["mean","charge"], 
                       yerr=data["err_mean","charge"], 
                       linestyle="none", 
                       marker="o", 
                       label="pyMBE", 
                       color="C1", 
                       fillstyle="none",
                       markeredgewidth=1.5)
        plt.xticks([2,4,6,8,10,12])
        
    elif fig_label in labels_fig8:   

        data=data.astype({("pH","value"): 'float'}).sort_values(by=("pH","value"))
        data=data[data.pdb.value == f'{protein_pdb}']


        plt.errorbar(data["pH"]["value"], 
                   data["mean","charge"], 
                   yerr=data["err_mean","charge"], 
                   linestyle="none", 
                   marker="o", 
                   label="pyMBE", 
                   color="C1", 
                   fillstyle="none",
                   markeredgewidth=1.5)
        plt.xticks([2,3,4,5,6,7])
        
    elif fig_label == "9":
        data=data.astype({("pH","value"): 'float'}).sort_values(by=("pH","value"))
        plt.errorbar(data["pH"]["value"], 
                       data["mean","alpha"], 
                       yerr=data["err_mean","alpha"], 
                       linestyle="none", 
                       marker="o", 
                       label="pyMBE", 
                       color="C1", 
                       fillstyle="none",
                       markeredgewidth=1.5)
        plt.xticks([2,4,6,8,10,12])

    # Save plot
    fig_path=pmb.get_resource("samples/Beyer2024")+"/figs"
    if not os.path.exists(fig_path):
        os.makedirs(fig_path)
    plt.legend(frameon=False, loc="lower left", fontsize=9, bbox_to_anchor=(0,1.02,1,0.2), mode="expand", borderaxespad=0, ncol=2)
    plt.savefig(f"{fig_path}/{fig_label}.pdf", 
                bbox_inches='tight')
    plt.close()
