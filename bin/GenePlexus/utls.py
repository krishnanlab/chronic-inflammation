import numpy as np
import pandas as pd
import pickle
from scipy.stats import hypergeom
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import average_precision_score
import time
from scipy.spatial.distance import cosine
import os

'''
This script contains all the function needed to excute the geneplexus pipeline
'''

################################################################################################################################


def intial_ID_convert(input_genes,file_loc='local'):
    '''
    This function takes in the input genes and then sees which ones can convert to Entrez IDs
    '''
    #load all the possible conversion dictionaries 
    convert_types = ['ENSG','Symbol','ENSP','ENST']
    all_convert_dict = {}
    for anIDtype in convert_types:
        convert_tmp = load_dict('to_Entrez',file_loc,anIDtype_=anIDtype)
        all_convert_dict[anIDtype] = convert_tmp
            
    # make some place holder arrays
    convert_IDs = [] # This will be a flat list for Entrez IDs to use as positives
    convert_out = [] # This will be a list of lists that will be used to tell user the conversions made
    for agene in input_genes:
        try:
            agene_int = int(agene)
            convert_out.append([agene_int,agene_int])
            convert_IDs.append(agene_int)
        except ValueError:
            for idx, anIDtype in enumerate(convert_types):
                if agene in all_convert_dict[anIDtype]:
                    convert_IDs = convert_IDs + all_convert_dict[anIDtype][agene]
                    convert_out.append([agene,', '.join(all_convert_dict[anIDtype][agene])])
                    break
                elif idx == len(convert_types)-1:
                    convert_out.append([agene,'Could Not be mapped to Entrez'])
    df_convert_out = pd.DataFrame(convert_out,columns=['Original_ID','ID_converted_to_Entrez'])
    df_convert_out = df_convert_out.astype({'Original_ID':str,'ID_converted_to_Entrez':str})
    return convert_IDs, df_convert_out
    
def make_validation_df(df_convert_out,file_loc='local'):
    '''
    This function make a dataframe of how the genes were converted for each network
    The ouptput is not used in the inflammation project
    '''
    table_summary = []
    #num_converted_to_Entrez = df_convert_out[~(df_convert_out['ID_converted_to_Entrez']=='Could Not be mapped to Entrez')].shape[0]
    input_count = df_convert_out.shape[0]
    converted_genes = df_convert_out['ID_converted_to_Entrez'].to_numpy()
    for anet in ['BioGRID','STRING','STRING-EXP','GIANT-TN']:
        net_genes = load_txtfile('net_genes',file_loc,net_type_=anet)
        df_tmp = df_convert_out[df_convert_out['ID_converted_to_Entrez'].isin(net_genes)]
        pos_genes_in_net = np.intersect1d(converted_genes,net_genes)
        table_row = {'Network': anet, 'NetworkGenes': len(net_genes), 'PositiveGenes': len(pos_genes_in_net)}
        table_summary.append(dict(table_row))
        tmp_ins = np.full(len(converted_genes),'N',dtype=str)
        tmp_ins[df_tmp.index.to_numpy()] = 'Y'
        df_convert_out['In %s?'%anet] = tmp_ins

    df_convert_out = df_convert_out.rename(columns = {'Original_ID': 'Original ID', 'ID_converted_to_Entrez': 'Entrez ID'})
    return df_convert_out, table_summary, input_count
        
def get_genes_in_network(convert_IDs,net_type,file_loc='local'):
    '''
    This function gets the genes in the choosen network and finds
    out which input_genes are in the network and lables thme as positives
    '''
    net_genes = load_txtfile('net_genes',file_loc,net_type_=net_type)
    pos_genes_in_net = np.intersect1d(np.array(convert_IDs),net_genes)
    genes_not_in_net = np.setdiff1d(np.array(convert_IDs),net_genes)
    return pos_genes_in_net, genes_not_in_net, net_genes
    
def get_negatives(pos_genes_in_net,net_type,GSC,file_loc='local'):
    '''
    This function finds which genes should be used as negatvies by
    1. Find all genes in the geneset collection that was choosen
    2. Remove genes that are considered postives
    3. Remove any genes in sets that are "too similar" to the input gene list
    '''
    uni_genes = load_txtfile('uni_genes',file_loc,net_type_=net_type,GSC_=GSC)
    good_sets = load_dict('good_sets',file_loc,GSC_=GSC,net_type_=net_type)
    M = len(uni_genes)
    N = len(pos_genes_in_net)
    genes_to_remove = pos_genes_in_net
    for akey in good_sets:
        n = len(good_sets[akey]['Genes'])
        k = len(np.intersect1d(pos_genes_in_net,good_sets[akey]['Genes']))
        pval = hypergeom.sf(k-1, M, n, N)
        if pval < 0.05:
            genes_to_remove = np.union1d(genes_to_remove,good_sets[akey]['Genes'])
    negative_genes = np.setdiff1d(uni_genes,genes_to_remove)
    return negative_genes
    
def run_SL(pos_genes_in_net,negative_genes,net_genes,net_type,features,file_loc='local'):
    '''
    1) Train the machine learning model and make predictions on the full set of genes
    2) Evaluate model performance by doing cross valaidation
    '''
    pos_inds = [np.where(net_genes==agene)[0][0] for agene in pos_genes_in_net]
    neg_inds = [np.where(net_genes==agene)[0][0] for agene in negative_genes]
    data = load_npyfile('data',file_loc,features_=features,net_type_=net_type)
    
    std_scale = StandardScaler().fit(data)
    data   = std_scale.transform(data)
    Xdata = data[pos_inds+neg_inds,:]
    ydata = np.array([1]*len(pos_inds) + [0]*len(neg_inds))
    clf = LogisticRegression(max_iter=10000,solver='lbfgs',penalty='l2',C=1.0)
    clf.fit(Xdata,ydata)
    mdl_weights = np.squeeze(clf.coef_)
    probs = clf.predict_proba(data)[:,1]
    
    if len(pos_genes_in_net) < 15:
        avgps = [-10, -10, -10]
    else:
        avgps = []
        n_folds = 3
        skf= StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=None)
        for trn_inds, tst_inds in skf.split(Xdata,ydata):
            clf_cv = LogisticRegression(max_iter=10000,solver='lbfgs',penalty='l2',C=1.0)
            clf_cv.fit(Xdata[trn_inds],ydata[trn_inds])
            probs_cv = clf_cv.predict_proba(Xdata[tst_inds])[:,1]
            avgp = average_precision_score(ydata[tst_inds],probs_cv)
            num_tst_pos = np.sum(ydata[tst_inds])
            prior = num_tst_pos/Xdata[tst_inds].shape[0]
            log2_prior = np.log2(avgp/prior)
            avgps.append(log2_prior)
        # avgp = '{0:.2f}'.format(np.median(avgps)) # used in webserver but not for inflamation work
    return mdl_weights, probs, avgps
    
def make_prob_df(net_genes,probs,pos_genes_in_net,negative_genes,file_loc='local'):
    '''
    Format the predictions from the machine learning model into a more
    convienent form
    '''
    Entrez_to_Symbol = load_dict('Entrez_to_Symbol',file_loc)
    Entrez_to_Name = load_dict('Entrez_to_Name',file_loc)
    prob_results = []
    for idx in range(len(net_genes)):
        if net_genes[idx] in pos_genes_in_net:
            class_label = 'P'
        elif net_genes[idx] in negative_genes:
            class_label = 'N'
        else:
            class_label = 'U'
        try:
            syms_tmp = '/'.join(Entrez_to_Symbol[net_genes[idx]]) #allows for multimapping
        except KeyError:
            syms_tmp = 'N/A'
        try:
            name_tmp = '/'.join(Entrez_to_Name[net_genes[idx]]) #allows for multimapping
        except KeyError:
            name_tmp = 'N/A'
        prob_results.append([net_genes[idx],syms_tmp,name_tmp,probs[idx],class_label])
    df_probs = pd.DataFrame(prob_results,columns=['Entrez','Symbol','Name','Probability','Class-Label'])
    df_probs = df_probs.astype({'Entrez':str,'Probability':float})
    df_probs = df_probs.sort_values(by=['Probability'],ascending=False)
    return df_probs, Entrez_to_Symbol
    
def save_files(fp_save,jobname,df_probs,avgps):
    'Save the predictions and the cross validation results'
    if not os.path.exists(fp_save):
        os.makedirs(fp_save)
    df_probs.to_csv(fp_save+jobname+'--predictions.tsv',sep='\t',header=True,index=False)
    np.savetxt(fp_save+jobname+'--CVvalues.txt',avgps,header='CVs (log2p)')
    
        
################################################################################################################################

'''
Below is code to load the various data files
'''
fp_HPCC = '/mnt/research/compbio/krishnanlab/tmp/to_alex_from_Chris/inflammation/chronic-inflammation/data_Zenodo/GenePlexus/'
def load_txtfile(file_type,file_loc,dtype_=str,net_type_=None,GSC_=None,target_set_=None):
    if file_type == 'net_genes':
        if file_loc == 'local':
            output_txt = np.loadtxt('../../data_Zenodo/GenePlexus/Node_Orders/%s_nodelist.txt'%net_type_,dtype=dtype_)
        elif file_loc == 'HPCC':
            output_txt = np.loadtxt(fp_HPCC + 'Node_Orders/%s_nodelist.txt'%net_type_,dtype=dtype_)
    elif file_type == 'uni_genes':
        if file_loc == 'local':
            output_txt = np.loadtxt('../../data_Zenodo/GenePlexus/GSCs/%s_%s_universe.txt'%(GSC_,net_type_),dtype=dtype_)
        elif file_loc == 'HPCC':
            output_txt = np.loadtxt(fp_HPCC + 'GSCs/%s_%s_universe.txt'%(GSC_,net_type_),dtype=dtype_)
    return output_txt

def load_npyfile(file_type,file_loc,features_=None,net_type_=None,GSC_=None,target_set_=None):
    if file_type == 'data':
        if file_loc == 'local':
            output_npy = np.load('../../data_Zenodo/GenePlexus/%s/%s_data.npy'%(features_,net_type_))
        elif file_loc == 'HPCC':
            output_npy = np.load(fp_HPCC + '%s/%s_data.npy'%(features_,net_type_))
    return output_npy
        
def load_dict(file_type,file_loc,anIDtype_=None,GSC_=None,net_type_=None,target_set_=None,features_=None):
    if file_type == 'to_Entrez':
        if file_loc == 'local':
            with open('../../data_Zenodo/GenePlexus/ID_conversion/Homo_sapiens__%s-to-Entrez__All-Mappings.pickle'%anIDtype_,'rb') as handle:
                output_dict = pickle.load(handle)
        elif file_loc == 'HPCC':
            with open(fp_HPCC + 'ID_conversion/Homo_sapiens__%s-to-Entrez__All-Mappings.pickle'%anIDtype_,'rb') as handle:
                output_dict = pickle.load(handle)
    elif file_type == 'good_sets':
        if file_loc == 'local':
            with open('../../data_Zenodo/GenePlexus/GSCs/%s_%s_GoodSets.pickle'%(GSC_,net_type_),'rb') as handle:
                output_dict = pickle.load(handle)
        elif file_loc == 'HPCC':
            with open(fp_HPCC + 'GSCs/%s_%s_GoodSets.pickle'%(GSC_,net_type_),'rb') as handle:
                output_dict = pickle.load(handle)
    elif file_type == 'Entrez_to_Symbol':
        if file_loc == 'local':
            with open('../../data_Zenodo/GenePlexus/ID_conversion/Homo_sapiens__Entrez-to-Symbol__All-Mappings.pickle','rb') as handle:
                output_dict = pickle.load(handle)
        elif file_loc == 'HPCC':
            with open(fp_HPCC + 'ID_conversion/Homo_sapiens__Entrez-to-Symbol__All-Mappings.pickle','rb') as handle:
                output_dict = pickle.load(handle)
    elif file_type == 'Entrez_to_Name':
        if file_loc == 'local':
            with open('../../data_Zenodo/GenePlexus/ID_conversion/Homo_sapiens__Entrez-to-Name__All-Mappings.pickle','rb') as handle:
                output_dict = pickle.load(handle)
        elif file_loc == 'HPCC':
            with open(fp_HPCC + 'ID_conversion/Homo_sapiens__Entrez-to-Name__All-Mappings.pickle','rb') as handle:
                output_dict = pickle.load(handle)          
    return output_dict


    
    
    
