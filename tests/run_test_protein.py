from matplotlib import colors
from matplotlib.axis import YAxis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse 
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import StrMethodFormatter

import os
import sys
import inspect

# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library()

parser = argparse.ArgumentParser(description='Script to reproduce the numerical data for globular proteins from Torres')
parser.add_argument('--run_command', type=str, required= True,  help='Run command for globular_protein.py')
parser.add_argument('--pdb_code', type=str, required= True,  help='PDB code of the protein')
args = parser.parse_args ()

pdb=args.pdb_code

#Run the main script for each pH value
for pH_value in np.arange(2.0, 7.5, 0.5):
    print (f'{args.run_command} --pH {pH_value}')
    os.system(f'{args.run_command} --pH {pH_value}')
    

if not os.path.exists(f'{pyMBE_path}/tests/observables_results'):
    os.system(f'mkdir {pyMBE_path}/tests/observables_results')

os.system(f'mv pH*.csv {pyMBE_path}/tests/observables_results')
os.system(f'python3 {pyMBE_path}/handy_scripts/data_analysis.py {pyMBE_path}/tests/observables_results/')

# Here you can adjust the width of the panda columns displayed when running the code 
pmb.pd.options.display.max_colwidth = 10

#Load the pmb.df 
pmb.df = pmb.read_pmb_df(filename='df.csv')
protein_sequence = pmb.df.loc[pmb.df['name']== pdb].sequence.values[0]
titratables_AA_df = pmb.df[['name','pka','acidity']].drop_duplicates().dropna()

# Get the pka_set
pka_set = {}
for index in titratables_AA_df.name.keys():
    name = titratables_AA_df.name[index]
    pka_value = titratables_AA_df.pka[index]
    acidity = titratables_AA_df.acidity[index]   
    pka_set[name] = {'pka_value':pka_value,'acidity':acidity}

# Calculate Henderson-Hasselbach titration curve 
pH_range = np.linspace(2, 7, num=31)
Z_HH = pmb.calculate_HH (object_name=pdb,
                        pH_list = pH_range, 
                        pka_set=pka_set )

# Here we have to add +2 for the Calcium in the protein charge by HH
if pdb == '1f6s':
    for index in range (len (Z_HH)):
        Z_HH[index] = Z_HH[index] +2 

# Read the reference data from Torres2022 

ref_data_torres = f'{pyMBE_path}/reference_data/{pdb}-10mM-torres.dat'
pH_list = []
znet_ref = []
sigma_ref = []

with open (ref_data_torres,'r') as file:
    for line in file: 
        line_split = line.split ()
        pH = float (line_split[0])
        pH_list.append(pH)
        znet_ref.append (float(line_split[1]))
        sigma_ref.append(float(line_split[2]))

# Read the results from espresso simulation using pyMBE 
espresso_data = 'analyzed_observables.csv'
full_data = pd.read_csv(espresso_data)
columns = full_data.columns.to_list()
full_data.sort_values('pH',ascending=True, inplace=True)

znet_espresso = full_data['Znet'].to_list()

numerical_comparison = pd.DataFrame()
numerical_comparison['pH'] = pH_list 
numerical_comparison['ref_torres'] = znet_ref
numerical_comparison['espresso'] = znet_espresso

numerical_comparison['error %'] = abs(( (numerical_comparison['espresso']) -  (numerical_comparison['ref_torres'])) /  (numerical_comparison['ref_torres'])) *100 

#Save `numerical_comparison` to a csv file

path_to_tests=pmb.get_resource("tests") 
numerical_comparison.to_csv(f'{path_to_tests}/{pdb}-numerical_comparison.csv',index = True)

#Plot results  

fig = plt.figure(figsize = (20,20))
ax1 = fig.add_subplot (111)

# Plots the HH equation
ax1.plot(
    pH_range,
    Z_HH,
    linewidth = 3,
    color = 'black',
    label = 'HH' 
    )

# Plot the ref data from Torres2022
ax1.errorbar (
    pH_list,
    znet_ref,
    yerr=sigma_ref,
    linewidth = 2,
    elinewidth = 3,
    marker = 'o',
    markersize = 10,
    markerfacecolor = 'none',
    alpha = 0.8, #changes the line opacity
    color = 'red',
    label = 'Ref Torres et al.')

#Plots the resuls from espresso

ax1.errorbar (
    full_data['pH'],
    full_data['Znet'],
    yerr = full_data['ZnetErr'],
    linewidth = 2,
    elinewidth = 3,
    marker = 's',
    markersize = 10,
    markerfacecolor = 'none',
    alpha = 0.8, #changes the line opacity
    color = 'blue',
    label = 'ESPResSo')

# Add axes information 
ax1.set_title ('Net charge vs pH. $c_{salt}$ = 0.01 M',fontsize ='40',pad = 30)
ax1.set_xlabel('$\it{pH}$',size = 45)
ax1.set_ylabel('$\it{Z}$', size = 45)

ax1.hlines(y=0,xmin=2,xmax=7,lw = 3, colors='grey',linestyles='dashed')

ax1.invert_yaxis ()
ax1.set_xlim ([2.0,7.0])
ax1.set_ylim ([-10.0,20.0])

y = np.linspace(-10,20,7)
minor_ticks = np.linspace(-10,20,31)

ax1.set_xticks(pH_list, minor=True)
ax1.set_yticks (y,minor=True)
ax1.set_yticks (minor_ticks,minor=True)

ax1.tick_params(axis='x',which='major',length=30,direction='in',width=3,colors='black',pad=15,labelsize=40)
ax1.tick_params(axis='x',which='minor',length=15,direction='in',width=3)

ax1.tick_params(axis='y',which='major',length=30,direction='in',width=3,colors='black',pad=10,labelsize=40) 
ax1.tick_params(axis='y',which='minor',length=15,direction='in',width=3,colors='black')

ax1.spines['left'].set_color('black')
ax1.spines['left'].set_lw(3)
ax1.spines['top'].set_lw(3)
ax1.spines['right'].set_lw(3)
ax1.spines['bottom'].set_lw(3)

ax1.legend(frameon =False)

plt.legend(prop={'size': 35})
pdf_name = f'{path_to_tests}/observables_results/{pdb}-analyzed_observables.pdf'
plt.savefig(pdf_name)
plt.show()


