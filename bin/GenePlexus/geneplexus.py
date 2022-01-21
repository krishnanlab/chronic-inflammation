import argparse
import utls
import numpy as np
import pandas as pd

'''
This function calls function in utls.py to step through the geneplexus pipeline
'''

    
def read_input_file(file_path,sep=', '):
    '''
    This function is specific to file formats used in the inflammation project
    '''
    input_genes = pd.read_csv(file_path,sep='\t')
    input_genes = input_genes['Gene'].tolist()
    input_genes = [str(item) for item in input_genes]
    print('The length of input genes is',len(input_genes))
    return input_genes

    
def run_model(input_genes, flie_loc, net_type, GSC, features, jobname, fp_save):

    print('1. Validating the input genes')
    convert_IDs, df_convert_out = utls.intial_ID_convert(input_genes,file_loc=flie_loc)
    df_convert_out, table_summary, input_count = utls.make_validation_df(df_convert_out,file_loc=flie_loc)
    
    print('2. Geting the genes in the choosen network')
    pos_genes_in_net, genes_not_in_net, net_genes = utls.get_genes_in_network(convert_IDs,net_type,file_loc=flie_loc)

    print('3. Finding negatives genes to use in the machine learning model')
    negative_genes = utls.get_negatives(pos_genes_in_net,net_type,GSC,file_loc=flie_loc)

    print('4. Training the machine learning model')
    mdl_weights, probs, avgps = utls.run_SL(pos_genes_in_net,negative_genes,net_genes,
                                            net_type,features,file_loc=flie_loc)

    print('5. Making a dataframe of the predictions')
    df_probs, Entrez_to_Symbol = utls.make_prob_df(net_genes,probs,pos_genes_in_net,negative_genes,file_loc=flie_loc)
    print('6. Saving the results')
    utls.save_files(fp_save,jobname,df_probs,avgps)

    

    
