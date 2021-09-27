import moana
import numpy as np
import pandas as pd

def corner(dir_mcmc, dir_out)
if __name__ == '__main__':

    # Path, file names and MCMC parameters
    path_mcmc = "./"
    run_name = 'mcmc_fit'
    N_burn = 200

    # --- Load MCMC file ---
    # Define column names
    colnames = np.array(['X1', 'Y1', 'X2', 'Y2', 'SEP', 'FRATIO', 'FTOTAL', 'F1', 'F2', 'chi2'])
    
    # Add flux for each dataset
    #datasets = ['k']
    #flux = np.array([f"fs_{a} fb_{a}".split(' ') for a in datasets])
    #colnames = np.append(colnames, flux)

    # Number of columns
    cols = range(len(colnames))
    print(cols)

    # Load MCMC output
    fname = f'{path_mcmc}/{run_name}.dat'
    print(f'Loading {fname}...')
    data = pd.read_table(fname, sep='\s+', names=colnames, usecols=cols, 
        dtype=np.float64, skiprows=N_burn)

    # --- Preparation of the Dataframe ---
    # Remove column if parameter is constant
    #data = data.loc[:, (data != data.iloc[0]).any()]

    # Add some missing quantities
    data['dchi2'] = data['chi2'] - np.min(data['chi2'])
    #data['tE'] = 1 / data['invtE']
    #data['q'] = moana.dbc.mass_fration_to_mass_ratio(data['eps1'])
    #data['qrescaled'] = 1e2 * data['q']

    # Reset Dataframe index
    #data.reset_index(drop=True, inplace=True)

    # Display some basic information
    print('Preview of first few lines...')
    print(data.head(3))
    print(f"Minimum chi-square: {data['chi2'].min():.2f}")
    print(f"Sample size: {len(data)}")

    # --- Create object to get statistical info, including scatter plot  ---
    # Choose here which parameters you want to include in the corner plot
    labels = ['X1', 'Y1', 'X2', 'Y2', 'SEP', 'FRATIO', 'FTOTAL', 'F1', 'F2']

    # Create the object
    posterior = moana.SampledPosterior(data, labels=labels)

    # --- Plot options ---

    # Define label names using LaTeX
    # If not defined, label names will be parameter name.
    column_labels = dict()
    column_labels.update({
        'sep': 's',
        'tE': r't_E\ (days)', 
        'u0': r'u_0',
        'q': 'q',
        'qrescaled': r'q/10^2',
    })

    # Range, major and minor ticks locators
    # -- Example 1:
    # {'tE': [[11.7, 17.1], [10.1, 19], [2, 4], [1, 2]]}
    # means for the subplots:
    # * tE_min = 11.7, tE_max = 17.1 when tE is represented on x-axis
    # * tE_min = 10.1, tE_max = 19 when tE is represented on y-axis
    # * Major ticks will be mutiples of 2 when tE is represented on x-axis
    # * Major ticks will be mutiples of 1 when tE is represented on x-axis
    # * There will be 4 minor ticks between each major ticks when tE is represented on x-axis
    # * There will be 2 minor ticks between each major ticks when tE is represented on y-axis
    # -- Example 2:
    # {'tE': [None, [10.1, 19], [2, 4], None]}
    # means that the choice of x range and y ticks will be automated.
    # TIP --> It is highly recommanded to start without specifying any ax_opts, because
    # if you ask for too many ticks, then the code may crash. If you do not
    # specify any ax_opts, then everything will be automated.
    #ax_opts = dict()
    #ax_opts.update({
    #    'u0': [[0.7, 1.25], [0.7, 1.25], [0.2, 4], [0.2, 4]],
    #    'tE': [[11.7, 17.1], [11.7, 17.1], [2, 4], [2, 4]],
    #    'sep': [[0.4, 0.53], [0.4, 0.53], [0.05, 5], [0.05, 5]],
    #    'qrescale': [[0.42, 0.75], [0.42, 0.75], [0.1, 4], [0.1, 4]],
    #})

    # Bins for 2-dimensional contours
    # -- Example 1:
    # {'u0tE': [100, 80]}
    # means for the subplot where u0 is on the x axis, and tE on the y axis:
    # * 100 bins will be used in the x axis direction
    # * 80 bins will be used in the y axis direction
    # TIP --> It is highly recommanded to start without specifying any bins.
    # If you do not specify any bins, then 80 bins in both x and y axes will
    # be used.
    bins = dict()
    bins.update({
        'u0sep': [1000, 1000],
        'u0tE': [1000, 1000],
        'tEsep': [100, 80],
        'sepqrescale': [100, 100]
    })

    # --- Plot the beautiful corner plot: EXAMPLE 1 ---
    print("mcmc_cdf...")
    fname = 'MCMC_CDF.pdf'
    posterior.corner_plot(filename=fname, labels=column_labels,
        rotation=30)

    # --- Plot the beautiful corner plot: EXAMPLE 4 ---
    # We want to see the samples and hide the 1-3 sigmas filled contours
    print("mcmc_chisq...")
    fname = 'MCMC_CHISQ.png'
    posterior.corner_plot(filename=fname, labels=column_labels,
        diagonal='chi2', show_samples=True, credible_intervals=False,
        rotation=30)
















