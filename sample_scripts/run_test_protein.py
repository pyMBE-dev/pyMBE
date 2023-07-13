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

parser = argparse.ArgumentParser(description='Script to run all the pH of the globular protein')
parser.add_argument('-pdb', type=str, required= True,  help='Expected PDB code of the protein')
args = parser.parse_args ()

protein_name = args.pdb

#Run the main script for each pH value
for pH_value in np.arange(2.0, 7.5, 0.5):
    print (f'Currently running {pH_value}')
    os.system(f'../pypresso globular_protein.py -pdb {protein_name} -pH {pH_value}')

if not os.path.exists('./observables_results'):
    os.system('mkdir observables_results')

os.system('mv pH*.csv observables_results')
os.system('python3 ../handy_scripts/data_analysis.py observables_results/')

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library()

# Here you can adjust the width of the panda columns displayed when running the code 
pmb.pd.options.display.max_colwidth = 10

#This line allows you to see the complete amount of rows in the dataframe
# pmb.pd.set_option('display.max_rows', None)


#Load the pmb.df 
pmb.df = pd.read_csv ('df.csv', header=[0],index_col=[0])
protein_sequence = pmb.df.loc[pmb.df['name']== protein_name].sequence.values[0]

print ("The protein sequence is:", protein_sequence)

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
Z_HH = pmb.calculate_HH (sequence = protein_sequence, pH_list = pH_range, pka_set=pka_set )

# Here we have to add +2 for the Calcium in the protein charge by HH
if protein_name == '1f6s':
    for index in range (len (Z_HH)):
        Z_HH[index] = Z_HH[index] +2 

# Read the reference data from Torres2022 

ref_data_torres = f'../reference_data/{protein_name}-10mM-torres.dat'
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

znet_espresso = full_data['Znet'].to_list ()

numerical_comparison = pd.DataFrame()
numerical_comparison['pH'] = pH_list 
numerical_comparison['ref_torres'] = znet_ref
numerical_comparison['espresso'] = znet_espresso

numerical_comparison['error %'] = abs(( (numerical_comparison['espresso']) -  (numerical_comparison['ref_torres'])) /  (numerical_comparison['ref_torres'])) *100 

#Save `numerical_comparison` to a csv file 
numerical_comparison.to_csv(f'{protein_name}-numerical_comparison.csv',index = True)

print (numerical_comparison)

#Plot results  

fig = plt.figure(figsize = (20,20))
ax1 = fig.add_subplot (111)

# Plots the HH equation
ax1.plot(
    pH_range,
    Z_HH,
    linewidth = 3,
    #marker = 'o',
    #markersize = 10,
    # markerfacecolor = 'none',
    #alpha = 0.8,
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

# ax1.set_xticks(5)

ax1.tick_params(axis='x',which='major',length=30,direction='in',width=3,colors='black',pad=15,labelsize=40)
ax1.tick_params(axis='x',which='minor',length=15,direction='in',width=3)

ax1.tick_params(axis='y',which='major',length=30,direction='in',width=3,colors='black',pad=10,labelsize=40) 
ax1.tick_params(axis='y',which='minor',length=15,direction='in',width=3,colors='black')

# # ax2.tick_params(axis='y',which='major',length=20,direction='in',width=3,colors='black',pad=10,labelsize=30)
# # ax2.tick_params(axis='y',which='minor',length=10,direction='in',width=3)

# ax1.yaxis.label.set_color('black')
ax1.spines['left'].set_color('black')
ax1.spines['left'].set_lw(3)
ax1.spines['top'].set_lw(3)
ax1.spines['right'].set_lw(3)
ax1.spines['bottom'].set_lw(3)

ax1.legend(frameon =False)

plt.legend(prop={'size': 35})
pdf_name = f'{protein_name}-analyzed_observables.pdf'
plt.savefig(pdf_name)
plt.show()


