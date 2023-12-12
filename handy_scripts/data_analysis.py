import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 
from pathlib import Path
import argparse
import sys

parser = argparse.ArgumentParser(description='Data analysis')
parser.add_argument('folder', type=str, help='Folder to read the files', default=sys.path[0])
parser.add_argument('--plot',help='Folder to read the files',default=False,action='store_true')
parser.add_argument('--columns_to_plot', type=str, help='Folder to read the files', default=None)

args = parser.parse_args ()

folder = args.folder 
plot = args.plot
columns_to_plot = args.columns_to_plot 


def plot_time_trace(full_data,filename, columns="all", mean_values = None, save_pdf = False):
    if(columns == "all"):
        columns = full_data.columns.to_list()
    for x in ['time', 'step', 'sweep']:
        if(x in columns):
            columns.remove(x)
    if(save_pdf):
            directory = "analysis_outputs"
            if(os.path.isdir(directory)==False):
                os.mkdir(directory)
    n_cols = len(columns)
    if 'time' in full_data.columns.to_list():
        x_data = full_data['time']
    elif 'MC sweep' in full_data.columns.to_list():
        x_data = full_data['MC sweep']
    else:
        raise ValueError("columns: "+str( full_data.columns.to_list() ) )
    for col in columns:
        fig = plt.figure(figsize = (20,7))
        plt.plot(
            x_data,
            full_data[col],
            linewidth = 1,
            marker = 'o',
            markersize = 3,
            alpha = 0.8,
            color = 'blue',
            label = col
        )
        plt.xlabel('time [internal units]')
        plt.ylabel(col)
        plt.title(filename)
        if(save_pdf):
            pdf_name = filename.name.replace("observables.dat", "_"+col+".pdf").replace("csv","pdf")
            plt.savefig(os.path.join(directory,pdf_name))
        plt.show()

def get_params_from_file_name(file_name):
    system_name = file_name.parts[-1].replace('_observables.dat', '')
    #print(system_name)
    entries = system_name.split('_')
    #print(entries)
    params = {}
    for entry in entries:
        sp_entry = entry.split('-')
        if(sp_entry[0]=='N'):
            params[sp_entry[0]] = int(sp_entry[-1])
        else:
            params[sp_entry[0]] = sp_entry[-1]
    return params

def rename_columns(data_frame, suffix, skip = ['time','step', "MC sweep"]):
    new_index = []
    for index in data_frame.index:
        if(index in skip):
            new_index.append(index)
        else:
            new_index.append(index+suffix)
    data_frame.index =  new_index
    return data_frame 

def get_dt(data):
    if 'time' in data:
        time = data['time']
    elif 'MC sweep' in data:
        time = data['MC sweep']
    else:
        raise ValueError("neither 'time' nor 'MC sweep' column found in data, got " + str(data))
    imax = data.shape[0]
    dt_init = time[1] - time[0]
    warn_lines = [];
    for i in range(1,imax):
        dt = time[i] - time[i-1]
        if(np.abs((dt_init - dt)/dt) > 0.01 ):
            warn_lines.append(f"Row {i} dt = {dt} = {time[i]} - {time[i-1]} not equal to dt_init = {dt_init}")
    if(len(warn_lines) > 20):
        print("\n")
        for line in warn_lines:
            print(line)
    return dt

def block_analyze(full_data,filename, n_blocks=16, equil=0.1):
    params = get_params_from_file_name(filename)    
    dt = get_dt(full_data) # check that the data was stored with the same time interval dt
    drop_rows = int(full_data.shape[0]*equil) # calculate how many rows should be dropped as equilibration
    # drop the rows that will be discarded as equlibration
    data = full_data.drop(range(0,drop_rows))
    # drop the columns step, time and MC sweep
    for column_name in ['time', 'step', 'MC sweep']:
        if(column_name in data.columns):
            data = data.drop(columns=column_name)
    # first, do the required operations on the remaining data
    n_samples = data.shape[0] # number of samples to be analyzed
    block_size = n_samples/n_blocks # mean block size
    mean = data.mean() # calculate the mean values of each column
    var_all = data.var() # calculate the variance of each column
    params['Blocks'] = n_blocks
    params['B_size'] = block_size
    print("b:", n_blocks, "k:", block_size)
    
    # calculate the mean per each block    
    blocks = np.array_split(data,n_blocks) # split the data array into blocks of (almost) equal size
    block_means = [] # empty list that we will use to store the means per each block
    for block in blocks:
        block_mean = block.mean() # mean values of each column within a given block
        block_means.append(block_mean)            
    block_means = pd.concat(block_means, axis=1).transpose()
    
    # perform calculations using averages or individual data blocks
    var_mean = block_means.var() # variance of the block averages = variance of the mean
    err_mean = np.sqrt(var_mean) # standard error of the mean
    n_eff = n_blocks*var_all/var_mean # effective number of samples in the whole simulation using eq.(28) and eq.(38) from Janke
    tau_int = dt*n_samples/(2*n_eff) # autocorrelation time using eq.(28) from Janke
    
    # modify the column names of the temporary results
    err_mean = rename_columns(err_mean,"Err")
    n_eff    = rename_columns(n_eff,"Nef")
    tau_int  = rename_columns(tau_int,"Tau")
    # first, concatenate the observables and alphabetically sort them by column names
    result = pd.concat( [ mean, err_mean, n_eff, tau_int ] ).sort_index(axis=0)
    # next, concatenate the results with simulation parameters to produce a human-readable output
    result = pd.concat( [ pd.Series(params), result] )
    return result

########################
# Key parameters

n_blocks = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048] # uncomment this line if you want to see how n_blocks affects the result
n_blocks = [16] # default number of blocks in the block analysis
equil = 0.1 # fraction of of all data that should be discarded as equilibration

filenames = list(Path(folder).glob('*.csv')) # file names to be analyzed

columns_to_plot = [columns_to_plot]

results = [] # empty list that we will use to store the results
print("filenames: ", filenames)

for n_b in n_blocks:
    for filename in filenames:
        print("filename: ", filename)
        full_data = pd.read_csv(filename) # read the full trajectory as pandas dataframe 

        if(plot):
            plot_time_trace(full_data,filename, columns=columns_to_plot, save_pdf = True) 
        tmp_result = block_analyze(
            full_data = full_data,
            filename = filename,
            n_blocks = n_b, 
            equil = equil, 
        )
        results.append(tmp_result)

results = pd.concat(results, axis=1).transpose() # re-shape the matrix with data to prepare it for writing to the csv file
print("final_results:\n", results)
results.to_csv(open("analyzed_observables.csv","w")) # save the results csv file that can be imported to spreadsheet calculators
print("###\nFinished\n###")