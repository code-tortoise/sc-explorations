import numpy as np # 1.20.3
import pandas as pd # 1.2.4
import scanpy as sc # 1.7.2
import os
import sys
from pathlib import Path
import scipy # 1.6.3
from scipy.sparse import csr_matrix
import scrublet as scr # 0.2.3
from statsmodels import robust # 0.12.2
import matplotlib.pyplot as plt # 3.4.2
import random
import datetime
import time

# Note These functions have not been tested with different versions and are not guarenteed to work if any of the package versions are different 

#to set unique date for outputs
calc_date = datetime.datetime.now()
date=calc_date.strftime('%Y-%m-%d')
#date=calc_date.strftime('%d-%m-%Y')
date = date.replace('-', '')



##########################################################################################################################################################

# Type checks

def string_check(x):
    if not isinstance(x,str):
        raise TypeError(f'Error: expected string arguemnt but "{x}" was inputted which is {type(x)}')

def numeric_check(x):
    if not isinstance(x,int):
        if not isinstance(x,float):
            if not isinstance(x,complex):
                raise TypeError(f'Error: expected numeric arguemnt but "{x}" was inputted which is {type(x)}')
                
def sequence_check(x):
    if not isinstance(x,list):
        if not isinstance(x,tuple):
            if not isinstance(x,range):
                raise TypeError(f'Error: expected sequence arguemnt but "{x}" was inputted which is {type(x)}')
                              
def mapping_check(x):
    if not isinstance(x,dict):
        raise TypeError(f'Error: expected mapping arguemnt but "{x}" was inputted which is {type(x)}')
                              
def set_check(x):
    if not isinstance(x,set):
        if not isinstance(x,frozenset):
             raise TypeError(f'Error: expected set arguemnt but "{x}" was inputted which is {type(x)}')
                              
def boolean_check(x):
    if not isinstance(x,bool):
        raise TypeError(f'Error: expected boolean arguemnt but "{x}" was inputted which is {type(x)}')
                              
def binary_check(x):
    if not isinstance(x,bytes):
        if not isinstance(x,bytearray):
            if not isinstance(x,memoryview):
                raise TypeError(f'Error: expected binary arguemnt but "{x}" was inputted which is {type(x)}')

##########################################################################################################################################################

# Data checks and formatting
                
def describe_basic(O):
    print(f'Overall shape: {O.shape}')
    print(f'Min count: {O.X.min()}')
    print(f'Max count: {O.X.max()}')

def describe_advanced(O):
    print(f'Overall shape: {O.shape}')
    print(f'Min count: {O.X.min()}')
    print(f'Max count: {O.X.max()}')
    print("Number of rows in obs: ", len(O.obs))
    print("Number of columns in obs: ", len(O.obs.columns))
    print("Number of rows in var: ", len(O.var))
    print("Number of columns in var: ", len(O.var.columns))
    if scipy.sparse.issparse(O.X):
        data = O.X.todense()
    else:
        data = O.X
    data = pd.DataFrame(data)
    data.columns = O.var.index
    data.index = O.obs.index
    plt.hist(data.loc[:,'FTL'])
    plt.title('FTL')
    plt.show()   

def simple_subset(data, var=False, subset_col = '', restrict = [''], restrict_index = False):
    string_check(subset_col)
    sequence_check(restrict)
    boolean_check(restrict_index)
    
    data_restrict = data
    
    if var == False:
        if restrict_index == True:
            data_restrict = data[data.obs.index.isin(restrict),:][:]
        else:
            data_restrict = data[data.obs[subset_col].isin(restrict),:][:]
    else:
        if restrict_index == True:
            data_restrict = data[:,data.obs.index.isin(restrict)][:]
        else:
            data_restrict = data[:,data.obs[subset_col].isin(restrict)][:]
    return data_restrict
   
##########################################################################################################################################################

# Directory functions

def save_path(save_path):
    if os.path.exists(save_path):
        print('File path provided exists')
        print('To save data to this location use the inputed variable at the start of the file path')
    else:
        print('File path provided does not exist!')
        save_path_answer = input(f'Do you want to create a file path at: "{save_path}"? \nIf so input "yes":') # Default of input is to return string
        print("")
        if (save_path_answer.lower() == 'yes') or (save_path_answer.lower() == 'y'):
            print(f'Selected to make a file path to: {save_path}')
            os.makedirs(save_path)
            print('To save data to this location use the inputed variable at the start of the file path')
        else:
            print('You did not select "yes". save_path not set or made')


def load(path, Method='', Transpose=False, cache = False):
    print('Beginning to load data:')
    string_check(path)
    print(f'Method chosen for loading is {Method} out of the options: h5ad , h5 or mtx')
    print()
    
    if (Path(path).is_file()):
        if 'h5ad' in Method:
            if (Transpose == False):
                object_to_load = sc.read(path)
            else:
                object_to_load = sc.read(path).T
        
        elif 'h5' in Method:
            if (Transpose == False):
                object_to_load = sc.read_10x_h5(path)
            else:
                object_to_load = sc.read_10x_h5(path).T
        
        elif 'mtx' in Method:
            raise Exception(f'Error in loading!\nSelected method mtx but provided path to a file! Please provide path to directory instead.')
        
        print(f'Data loaded in successfully')
        
    elif (os.path.isdir(path)):
        if 'mtx' in Method:
            if (Transpose == False):
                object_to_load = sc.read_10x_mtx(path, var_names = 'gene_symbols', cache=cache)
            else:
                object_to_load = sc.read_10x_mtx(path, var_names = 'gene_symbols', cache=cache)
        
        elif ('h5ad' in Method) or ('h5' in Method):
            raise Exception(f'Error in loading!\nSelected method {Method} but provided path to a directory not a file! Please provide path to directory instead.')
        
        print(f'Data loaded in successfully')
        
    else:
        raise Exception(f'Error in loading!\nMake sure you are using one of the following in Method: h5ad , h5 or mtx\nIf trying to load using methods h5ad or h5 object make sure path provided is the whole directory\nIf trying to load using method mtx make sure path is to the directory holding a matrix.mtx , features.tsv and barcodes.tsv')
        
    describe_basic(object_to_load)
    print("")
    return object_to_load


##########################################################################################################################################################

# QC functions

# To actually calculate scrublet using defined threshold cutoff
def simple_scrublet_run(data, column = '', cutoff = 3, save_path = save_path, date = date):

    string_check(column)
    numeric_check(cutoff)
    Cutoff = str(cutoff)
    
    if (os.path.exists(save_path+'scrublet/') == True):
        pass
    else:
        print(f'File path does not exists, making scrublet folder at the end of: {save_path}')
        os.makedirs(save_path+'scrublet/')
    
    if (os.path.exists(save_path+'scrublet/'+'Srub_run_val_' + Cutoff) == True):
        pass
    else:
        print(f'Making path for scrublet run output using cutoff: {Cutoff}')
        os.makedirs(save_path+'scrublet/'+'Srub_run_val_' + Cutoff)
        print('File path made')
    
    meta_10x_channels = column
    RUNs, DSs, CELLs, THRs, MEDs, MADs, CUTs, no_thr = [], [], [], [], [], [], [], []
    orig_stdout = sys.stdout
    sys.stdout = open(save_path+ "scrublet/" + 'Srub_run_val_' + Cutoff + '/' + 'scrublet_output_file_mad_' + date + '.txt', 'w')
    for run in data.obs[meta_10x_channels].unique():
        print(run)
        ad = data[data.obs[meta_10x_channels] == run, :]
        x = ad.X
        scrub = scr.Scrublet(x)
        ds, prd = scrub.scrub_doublets()
        RUNs.append(run)
        DSs.append(ds)
        CELLs.append(ad.obs_names)
        # MAD calculation of threshold:
        MED = np.median(ds)
        MAD = robust.mad(ds)
        CUT = (MED + (cutoff * MAD))
        MEDs.append(MED)
        MADs.append(MAD)
        CUTs.append(CUT)
        try:  # not always can calculate automatic threshold
            THRs.append(scrub.threshold_)
            print('Threshold found by scrublet')
        except:
            THRs.append(0.4)
            no_thr.append(run)
            print('No threshold found, assigning 0.4 to', run)
            scrub.call_doublets(threshold=0.4) # so that it can make the plot
        fig = scrub.plot_histogram()
        fig[0].savefig(save_path+ "scrublet/" + 'Srub_run_val_' + Cutoff + '/' + run + '_' + date + '.png')
        # Alternative histogram for MAD-based cutoff
        scrub.call_doublets(threshold=CUT)
        fig = scrub.plot_histogram()
        fig[0].savefig(save_path+ "scrublet/" + 'Srub_run_val_' + Cutoff + '/' + run + '_mad_' + date + '.png')
        plt.close('all')
        print()
        print()
    print()
    print('The following sample(s) do not have automatic threshold:')
    print(no_thr)
    sys.stdout.close()
    sys.stdout = orig_stdout
    ns = np.array(list(map(len, DSs)))
    tbl = pd.DataFrame({
        'run': np.repeat(RUNs, ns),
        'ds': np.concatenate(DSs),
        'thr': np.repeat(THRs, ns),
        'mad_MED': np.repeat(MEDs, ns),
        'mad_MAD': np.repeat(MADs, ns),
        'mad_thr': np.repeat(CUTs, ns),
        }, index=np.concatenate(CELLs))
    tbl['auto_prd'] = tbl['ds'] > tbl['thr']
    tbl['mad_prd'] = tbl['ds'] > tbl['mad_thr']
    tbl.to_csv(save_path + 'scrublet/' + 'Srub_run_val_' + Cutoff + '/' + 'doublets_score_mad_' + date + '.csv', header=True, index=True)
    #mad_prd = 'cutoff_' + Cutoff + '_mad_prd'
    #auto_prd = 'cutoff_' + Cutoff + '_auto_prd'
    is_doublet = 'cutoff_' + Cutoff + '_is_doublet'
    #data.obs[mad_prd] = tbl['mad_prd'].astype("str")
    #data.obs[auto_prd] = tbl['auto_prd'].astype("str")
    #data.obs[is_doublet] = data.obs[mad_prd]
    data.obs[is_doublet] = tbl['mad_prd'].astype("str")
    doublets_count = str(data.obs[is_doublet].value_counts()['True'])
    not_doublets_count = str(data.obs[is_doublet].value_counts()['False'])
    print(f'The number of doublets using column {column} at cutoff {Cutoff}:')
    print(f'Doublets: {doublets_count}')
    print(f'Not Doublets: {not_doublets_count}')
    print()
    print(f'Scrublet run using column {column} at cutoff {Cutoff} complete')
    print()


# Uses the simple_scrublet_run def on loop and outputs cleaner to understand
def scrublet(data, column = '', single_cutoff = False, multi_cutoffs = False, save_path = save_path, date = date):
    
    print('Starting scrublet run')
    
    if not single_cutoff is False:
        print(f'Single cutoff selected at value: {single_cutoff}')
        numeric_check(single_cutoff)
        simple_scrublet_run(data, column = column, cutoff = single_cutoff, save_path = save_path, date = date)
        
    elif not multi_cutoffs is False:
        print('Multiple cutoffs selected, proceeding to cycle through different cutoffs inputted')
        sequence_check(multi_cutoffs)
        print('Cutoff values inputted:')
        for i in multi_cutoffs:
            print(i)
            numeric_check(i)
        print()
        for i in multi_cutoffs:
            print(f'Starting scrublet run for value {i}')
            simple_scrublet_run(data, column = column, cutoff = i, save_path = save_path, date = date)
            print('Scrublet run finished')
            print()
            print()
    else:
        raise Exception('No cutoff selected!')
        

# To calculate impact of mito cutoff on data
def mito_calc(data, data_name = '', column='', subset_cat = '',mito_column = '', mito_cutoff=20, Whole_data = True):
    if (column == 'all'):
        lanes = len(data.obs)
    else:
        lanes = len(data.obs[column].unique())
        data = data[data.obs[column].isin([subset_cat])]
    before = len(data.obs)
    mito_val = mito_cutoff * 0.01
    after = len(data.obs[data.obs[mito_column] < mito_val])
    df_append = pd.DataFrame(columns=['Data','Subset_col','Number_unique_categories','Category_of_interest','Number_cells_before_cutoff','Percent_cuttoff','Number_cells_after_cutoff','Total_number_cells_removed','Percentage_loss'])
    df_append['Data'] = data_name
    if (Whole_data == True):
        df_append['Subset_col'] = ['None']
    else:
        df_append['Subset_col'] = [column]
    df_append['Number_unique_categories'] = [lanes]
    df_append['Category_of_interest'] = [subset_cat]
    df_append['Number_cells_before_cutoff'] = [before]
    df_append['Percent_cuttoff'] = [mito_cutoff]
    df_append['Number_cells_after_cutoff'] = [after]
    df_append['Total_number_cells_removed'] = [before-after]
    df_append['Percentage_loss'] = [round((100-((after/before)*100)),2)]
    return df_append

# To use summarise output of mito_calc in a easy to understand table format 
def mito_calc_summary(data, data_name = [''], columns=[''], mito_column = '', mito_cutoff=20, Whole_data = True):
    mito_df = pd.DataFrame(columns=['Data','Subset_col','Number_unique_categories','Category_of_interest','Number_cells_before_cutoff','Percent_cuttoff','Number_cells_after_cutoff','Total_number_cells_removed','Percentage_loss'])
    mito_append = mito_calc(data, data_name=data_name, column='all', subset_cat='None', mito_column=mito_column, mito_cutoff=mito_cutoff, Whole_data = Whole_data)
    mito_df = pd.concat([mito_df, mito_append], axis=0)
    Whole_data = False
    for col in columns:
        for x in data.obs[col].unique():
            mito_append = mito_calc(data, data_name = data_name, column = col, subset_cat=x, mito_column = mito_column, mito_cutoff=mito_cutoff, Whole_data = Whole_data)
            mito_df = pd.concat([mito_df, mito_append], axis=0)
    return mito_df


##########################################################################################################################################################

# Data transformations

def gene_removal(data, genes_to_remove):
    data = data[:, ~data.var_names.isin(genes_to_remove)].copy()
    return data


def transform_data(data, how='',  gene_removal_list = [] ,count_layer=False, hvg_min_mean = 0.0125, hvg_max_mean = 3, batch_key_input=None , subset_hvg = False, scale_zero_center=True, sacle_max_value=10, warn_check=False):
    
    if warn_check == True:
        for c in how:
            if (c in 'frnlhs') != True:
                raise Exception('Error: inputted character not expected. \nUse one of the following: f,n,l,h,s')
    
    Transformations_ran = []
    
    start_time = time.time()
    
    if "f" in how:
        sc.pp.filter_genes(data, min_cells=5)
        Transformations_ran.append('Filtered low quality genes')
        
    if "r" in how:
        if gene_removal_list:
            data = gene_removal(data,gene_removal_list)
            Transformations_ran.append('Removed genes from list provided by user')
        else:
            raise Exception("Error: trying to remove genes without passing a list of genes to remove")
    
    if count_layer == True:
        data.layers['raw_counts'] = data.X.copy()
        Transformations_ran.append('Made a copy of .X in layers called raw_counts')
        
    if "n" in how:
        sc.pp.normalize_per_cell(data, counts_per_cell_after=1e4)
        Transformations_ran.append('Normalised per cell to 10000')
        
    if "l" in how:
        sc.pp.log1p(data)
        Transformations_ran.append('Natural Logged')
        
    if "h" in how:
        if batch_key_input == None:
            if subset_hvg == True:
                sc.pp.highly_variable_genes(data, min_mean = hvg_min_mean, max_mean = hvg_max_mean, subset=True)
                Transformations_ran.append('Gene feature selection and data subsetted to gene features')
            else:
                sc.pp.highly_variable_genes(data, min_mean = hvg_min_mean, max_mean = hvg_max_mean)
                Transformations_ran.append('Gene feature selection')
        else:
            if subset_hvg == True:
                sc.pp.highly_variable_genes(data, min_mean = hvg_min_mean, max_mean = hvg_max_mean, batch_key = batch_key_input, subset=True)
                Transformations_ran.append('Gene feature selection by batch and data subsetted to gene features')
            else:
                sc.pp.highly_variable_genes(data, min_mean = hvg_min_mean, max_mean = hvg_max_mean, batch_key = batch_key_input)
                Transformations_ran.append('Gene feature selection by batch')
                
    if "s" in how:
        sc.pp.scale(data, zero_center=scale_zero_center, max_value=sacle_max_value)
        Transformations_ran.append('Scaled')
        
    finish_time = time.time() - start_time
    print(f'Data transformation complete. \nRuntime: {str(finish_time)}s')
    print()
    print(f'Options performed in order: {Transformations_ran}')
    return data
    
    
    
#def transform_data(data, how='', count_layer=False, zero_center=True, max_value=10, batch_key=None ,warn_check=False):
#    
#    if warn_check == True:
#        for c in how:
#            if (c in 'fnlhs') != True:
#                raise Exception('Error: inputted character not expected. \nUse one of the following: f,n,l,h,s')
#    
#    Transformations_ran = []
#    
#    start_time = time.time()
#    if "f" in how:
#        sc.pp.filter_genes(data, min_cells=5)
#        Transformations_ran.append('Filtered low quality genes')
#    if count_layer == True:
#        data.layers['raw_counts'] = data.X.copy()
#        Transformations_ran.append('Made a copy of .X in layers called raw_counts')
#    if "n" in how:
#        sc.pp.normalize_per_cell(data, counts_per_cell_after=1e4)
#        Transformations_ran.append('Normalised per cell to 10000')
#    if "l" in how:
#        sc.pp.log1p(data)
#        Transformations_ran.append('Natural Logged')
#    if "h" in how:
#        if batch_key == None:
#            sc.pp.highly_variable_genes(data)
#            Transformations_ran.append('Gene feature selection')
#        else:
#            sc.pp.highly_variable_genes(data, batch_key = batch_key)
#            Transformations_ran.append('Gene feature selection by batch')
#    if "s" in how:
#        sc.pp.scale(data, zero_center=zero_center, max_value=max_value)
#        Transformations_ran.append('Scaled')
#    finish_time = time.time() - start_time
#    print(f'Data transformation complete. \nRuntime: {str(finish_time)}s')
#    print()
#    print(f'Options performed in order: {Transformations_ran}')