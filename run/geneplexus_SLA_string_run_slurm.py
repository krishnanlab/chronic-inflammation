import numpy as np
import pandas as pd
import os
import pickle
import time
import glob

# these can be changed
net_type = 'STRING'
features = 'Adjacency'
save_path = '/mnt/gs18/scratch/users/mckimale/chronic/chronic-inflammation/results/GenePlexus_output/'
if not os.path.exists(save_path):
    os.makedirs(save_path)
# dir to put the slurm files
slurm_dir = '/mnt/gs18/scratch/users/mckimale/chronic/chronic-inflammation/run/Geneplexus_SLA_string/' 
if not os.path.exists(slurm_dir):
    os.makedirs(slurm_dir)

# get all the files and then get the part we care about for file name
files_to_do = glob.glob('/mnt/gs18/scratch/users/mckimale/chronic/chronic-inflammation/data/disease_gene_files/*.txt')
# add how to do background selection
# first choose a amin background source, then add files to switch to the other
bg_main = 'DisGeNet'
other_bg_files = ['chronic_inflammation_go']
# automatically set the other background below
if bg_main == 'DisGeNet':
    bg_other = 'GO'
elif bg_main == 'GO':
    bg_other = 'DisGeNet'
# set the job names and type of background to do
backgrounds = [bg_main if item.strip().split('/')[-1].split('.t')[0] not in other_bg_files else bg_other for item in files_to_do]
jobnames = [item.strip().split('/')[-1].split('.t')[0]+'--%s--%s--%s'%(net_type,features,backgrounds[idx]) for idx, item in enumerate(files_to_do)]
print('The number of files to do is',len(files_to_do))


for idx, afile in enumerate(files_to_do):

    mylist = ['#!/bin/bash']
    mylist.append('### define resources needed:')
    mylist.append('#SBATCH --time=00:30:00')
    mylist.append('#SBATCH --nodes=1')
    mylist.append('#SBATCH --mem=50G')
    mylist.append('#SBATCH --cpus-per-task=1')
    mylist.append('#SBATCH --job-name=%s'%jobnames[idx])
    mylist.append('#SBATCH --account=wang-krishnan')
    mylist.append('#SBATCH --output=%sslurm-%%x-%%j.out'%slurm_dir)
    mylist.append('module load powertools')
    mylist.append('umask g+rw')
    mylist.append('export PATH="/mnt/home/mckimale/anaconda3/bin:$PATH"') #### change this
    mylist.append('cd /mnt/gs18/scratch/users/mckimale/chronic/chronic-inflammation/bin/GeneplexusPublic')
    mylist.append('python example_run.py -i %s -j %s -n %s -f %s -g %s -s %s'%(afile,jobnames[idx],
                                                                               net_type,features,
                                                                               backgrounds[idx],save_path))

    with open(slurm_dir + 'Geneplexus-%s.sb'%jobnames[idx], 'w') as thefile:
        for item in mylist:
            thefile.write("%s\n" % item)
    os.system('sbatch ' + slurm_dir + 'Geneplexus-%s.sb'%jobnames[idx])
