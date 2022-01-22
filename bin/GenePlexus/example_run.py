import geneplexus
import numpy as np
import argparse


'''
This script allows a user to enter command line arguments to run the geneplexus pipeline
'''

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',
                    default = '../../data/disease_gene_files/Alzheimers_Disease.txt',
                    type = str,
                    help = 'file path to the gene functions')
parser.add_argument('-fl','--file_loc',
                    default = 'local',
                    type = str,
                    help = 'options local, HPCC, cloud')
parser.add_argument('-n','--net_type',
                    default = 'STRING',
                    type = str,
                    help = 'options are BioGRID, STRING-EXP, STRING, GIANT-TN')
parser.add_argument('-f','--features',
                    default = 'Adjacency',
                    type = str,
                    help = 'options are Embedding, Adjacency, Influence')
parser.add_argument('-g','--GSC',
                    default = 'GO',
                    type = str,
                    help = 'options are GO, DisGeNet')
parser.add_argument('-j','--jobname',
                    default = 'mynameargparse_AZ',
                    type = str,
                    help = 'sting to use for jobname')
parser.add_argument('-s','--fp_save',
                    default = '../../results/GenePlexus_output/',
                    type = str,
                    help = 'The path to save the file to (needs the / in it for now)')
args = parser.parse_args()

# functions to run geneplexus script
input_genes = geneplexus.read_input_file(args.input)
geneplexus.run_model(input_genes,args.file_loc,args.net_type,args.GSC,args.features,args.jobname,args.fp_save)
