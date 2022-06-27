import numpy as np

def get_system_name_cph(seq, pH, c_acid, c_salt, n_chains, chain_length ):
    return "Seq-{:s}_pH-{:.2f}_ca-{:.4f}_cs-{:.4f}_chains-{:.0f}_N-{:.0f}".format(seq, 
            pH, 
            c_acid.to("mol/L").magnitude, 
            c_salt.to('mol/L').magnitude,
            n_chains,
            chain_length
            )

def get_system_name_peptide(seq, model, pH, c_pep, c_salt, n_chains):
    return "Seq-{:s}_mod-{:s}_pH-{:.2f}_cp-{:.4f}_cs-{:.4f}_chains-{:.0f}".format(seq,
            model,
            pH, 
            c_pep.to("mol/L").magnitude, 
            c_salt.to('mol/L').magnitude,
            n_chains,
            )

def get_params_from_file_name(filename):
    return_dict = {}
    items = filename.rstrip('.csv').split('_')
    for item in items:
        values = item.split('-')
        return_dict[ values[0] ] = values[1]
    return return_dict

# the reference data from Henderson-Hasselbalch equation
def ideal_alpha(pH, pKa):
    return 1. / (1 + 10**(pKa - pH))

# statistical analysis of the results
def block_analyze(input_data, n_blocks=16):
    data = np.array(input_data)
    block = 0
    # this number of blocks is recommended by Janke as a reasonable compromise
    # between the conflicting requirements on block size and number of blocks
    block_size = int(data.shape[0] / n_blocks)
    print("block_size:", block_size)
    # initialize the array of per-block averages
    block_average = np.zeros((n_blocks, data.shape[0]))
    # calculate averages per each block
    for block in range(0, n_blocks):
        block_average[block] = np.average(data[block * block_size: (block + 1) * block_size])
    # calculate the average and average of the square
    av_data = np.average(data)
    av2_data = np.average(data * data)
    # calculate the variance of the block averages
    block_var = np.var(block_average)
    # calculate standard error of the mean
    err_data = np.sqrt(block_var / (n_blocks - 1))
    # estimate autocorrelation time using the formula given by Janke
    # this assumes that the errors have been correctly estimated
    tau_data = np.zeros(av_data.shape)
    if av_data == 0:
        # unphysical value marks a failure to compute tau
        tau_data = -1.0
    else:
        tau_data = 0.5 * block_size * n_blocks / (n_blocks - 1) * block_var / (av2_data - av_data * av_data)
    return av_data, err_data, tau_data, block_size

