# A network-based approach for isolating the chronic inflammation gene signatures underlying complex diseases towards finding new treatment opportunities
This GitHub repository contains all code used for reproducing results from
the manuscript "**A network-based approach for isolating the chronic inflammation gene signatures underlying complex diseases towards finding new treatment opportunities**", which can be found [here] (Get link)

This markdown documents will provide instructions on how to run the code for a sample disease and how to recreate the results.

### General requirements

1. unix or unix-like OS
2. Anaconda3 distribution
3. R version >=4.0.0 
4. Slurm workload manager (To recreate project results)

R libraries needed:

1. tidyverse 1.3.1
2. parallel 
3. mccf1 1.1
4. grid
5. org.Hs.eg.db 3.12.0
6. igraph 1.2.9
7. topGO 2.42.0

### Zenodo Download
Required data, as well as copies of pertinent results, are included on our [Zenodo record](https://zenodo.org/record/5919086/)

Note that this record contains ~45 GB.

These can be downloaded with the script `get_data.sh`. Run this script in
the repo so the Zenodo folder appears in this default filepath.

This folder, `data_Zenodo`, will include:

- `GenePlexus`: A local instance of GenePlexus
- `GenePlexus_parameter_checks`: Output for running GenePlexus with different arguments
- `GenePlexus_String_Adjacency`: Record of output from pipeline and analysis
- `pascal_out`: Pascal output files for UK BioBank traits we used
- `clinical_trials`: Clinical trial data for analysis
- `biogrid`: Biogrid network
- `ConsensusPathDB`: ConsensusPathDB network
- `string`: String network
- `string-exp`: String-exp network
- `prediction_clusters_same_graph`: Compendium of results to create paper figures

### data directory
The 'data' directory has some required data. This includes:

- `disease_gene_files`: Seed gene lists, generated in the pipeline,
- `drugcentral`: Data from DrugCentral and scripts used to format it.
- `dgidb`: Data from DGIdb and scripts used to download/format it.


### src directory
Scripts used to run the pipeline are located in `src`. 

Contains ```chronic_inflammation_functions.R```, which contains utility functions used by most scripts in the pipeline.

### run directory
`run` contains slurm submission scripts that were used to do our analysis.
This readme has instructions for how to run each script without use of slurm
for a sample disease

### figures directory
`figures` contains Markdown notebooks used to analyze final results. These results
can be found in the Zenodo record

## Pipeline instructions
The instructions in this readme will be how to run the pipeline with one sample
disease. Running with all traits used in this project requires the use of slurm
workload manager and is impractical without it.

### Format seed genes
__Script__:  
`prep_disease_gene_dfs.R`

__Purpose__:
Creates the seed gene files in `data/disease_gene_files` for specified diseases from
Disgenet.

__Arguments__:

1. Text file with one column (no header) containing the disease ids of interest from disgenet
2. Directory where "disease_gene_files" folder will go

__Run__:
```bash
Rscript prep_disease_gene_dfs.R \
 ../data/chronic_inflammation_diseases_non-ovlp_cuid.txt \
 ../data
```
### Getting inflammation genes from human genome
__Script__:  
`getInflammationGenesFrom_org.HS.eg.db.R`

__Purpose__:
Takes the human genome and gets genes from inflammation related GO terms

__Arguments__:
N/A

__Run__:
```bash
Rscript getInflammationGenesFrom_org.HS.eg.db.R
```

### Creating network edgelists
__Script__:  
`prepEdgelist.R`

__Purpose__:
Formats and creates an edgelist and Rdata object for a given network

__Arguments__:

1. Path to tab delimited edgelist
2. Path to output dir
3. Network name
4. True/False, keep edge weights or not

__Run__:
```bash
Rscript prepForClusterSaverunner.R \
 ../data_Zenodo/biogrid/biogrid_entrez_edgelist.txt \
 ../data_Zenodo/biogrid/ \
 bioGRID \
 TRUE

```

### Getting UK Biobank seed genes
__Script__:  
`getNegativeControls.R`

__Purpose__:
Creates the seed gene files in `data/disease_gene_files` for the UK BioBank
traits used in this project.

__Arguments__:

1. File from Zenodo of UK BioBank traits of interest
2. Location of Pascal output in Zenodo for each trait
3. location of `disease_gene_files`, where files will be output

__Run__:
```bash
Rscript getNegativeControls.R \
 ../data_Zenodo/our_ukbb_traits_description.tsv \
 ../data_Zenodo/pascal_out \
 ../data/
```

### Running Geneplexus
__Script__:  
`bin/GenePlexus/example_run.py`

__Purpose__:
Runs GenePlexus on a trait of interest and output the results. This project
used `ConsensusPathDB`,`Adjacency`, and `DisGeNet` for its final results

__Arguments__:  

-i : Disease seed genes \
-j : Job name \
-n : Network, options are BioGRID, STRING-EXP, STRING, ConsensusPathDB \
-f : Features, options are Embedding, Adjacency, Influence \
-g : GSC type, options are GO or DisGeNet \
-s : Output directory \
-fl : Option for how to run, always use `local` for this project \

__Run__:
```bash
python example_run.py \
 -i ../../data/disease_gene_files/Chronic_Obstructive_Airway_Disease.txt \
 -j Chronic_Obstructive_Airway_Disease--ConsensusPathDB--Adjacency--DisGeNet \
 -n ConsensusPathDB \
 -f Adjacency \
 -g DisGeNet \
 -s ../../results/GenePlexus_output/ \
 -fl local
```

### Summarize Geneplexus output
__Script__:  
`summarizeGeneplexusPredictions.R`

__Purpose__:
Returns multiple figures showing results for the network combination,
along with a summarized Rdata file that has pertinent disease results used in
later parts of the pipeline

__Arguments__:

1. Path to directory with GenePlexus predictions
2. Output directory
3. Average cv threshold

__Run__:
```bash
Rscript summarizeGeneplexusPredictions.R \
  ../results/GenePlexus_output/ \
  ../results/GenePlexus_parameters \
  1.0
```


### Clustering Geneplexus results
__Script__:  
`filterAndClusterGeneplexusPredictions.R`

__Purpose__:
Takes the GenePlexus predictions and assign genes to clusters
for each disease.

__Arguments__:

1. GenePlexus prediction path
2. Prediction threshold, either `mccf1` or a number < 1
3. Path to igraph object containing network for clustering
4. Leiden algorithm partition type
5. Resolution parameter
6. GenePlexus results path
7. True/False, Is the network weighted?

__Run__:
```bash
Rscript filterAndClusterGeneplexusPredications.R \
 ../results/GenePlexus_output/Chronic_Obstructive_Airway_Disease--bioGRID--Adjacency--DisGeNet--predictions.tsv \
 0.8 \
 ../data_Zenodo/biogrid/BioGrid_igraph.Rdata \
 ModularityVertexPartition \
 0.1 \
 ../results/prediction_clusters_same_graph
 FALSE

```


### Clustering inflammation genes
__Script__:  
`clusterInflammationGenes.R`

__Purpose__:
Clustering the inflammation genes

__Arguments__:

1. Path to inflammation genes
2. Path to igraph object that has network for clustering
3. Partition type
4. Resolution parameter
5. Results path
6. True/False, Is the network weighted?

__Run__:
```bash
Rscript clusterInflammationGenes.R \
 ../data/disease_gene_files/inflammation_genes/chronic_inflammatory_response_GO2ALLEGS.txt \
 ../data_Zenodo/biogrid/BioGrid_igraph.Rdata \
 ModularityVertexPartition \
 0.1
 ../results/prediction_clusters_same_graph/
 FALSE

```


### Clustering random genes
__Script__:  
`clusterRandomGenes.R`

__Purpose__:
Takes the 5000 fake traits that were randomly generated from a
disease and assigns the genes to clusters

__Arguments__:

1. Path to data containing all fake traits generated
2. Disease of interest
3. Path to igraph object containing network for clustering
4. Partition type
5. Resolution parameter
6. Results path
7. True/False, Is the network weighted?

__Run__:
```bash
Rscript clusterRandomGenes.R \
 ../data_Zenodo/5000Expandedfaketraits_biogrid.tsv \
 Chronic_Obstructive_Airway_Disease \
 ../data_Zenodo/biogrid/BioGrid_igraph.Rdata \
 ModularityVertexPartition \
 0.1
 ../results/prediction_clusters_same_graph/clusters/predicted_withBioGRID--clustered_on_BioGRID
 FALSE

```


### Get GOBP enriched clusters from GenePlexus resuls
__Script__:  
`find_GOBP_enriched_clusters_GenePlexus.R`

__Purpose__:
Finds GOBPs that are enriched in each cluster of a disease

__Arguments__:

1. Path to cluster file
2. Background genes from network
3. Output directory

__Run__:
```bash
Rscript find_GOBP_enriched_clusters_GenePlexus.R \
 ../data_Zenodo/prediction_clusters_same_graph/clusters/predicted_withConsensusPathDB--clustered_on_ConsensusPathDB/Malignant_neoplasm_of_pancreas--threshold--0.8--PredictionGraph--ConsensusPathDB--ClusterGraph--ConsensusPathDB_clusters.csv \
 ../data_Zenodo/ConsensusPathDB/ConsensusPathDB_genes.csv \
 ../results/prediction_clusters_same_graph/GOBP_enrichment

```
### Get GOBP enriched clusters for inflammation clusters
__Script__:  
`find_GOBP_enriched_inflammation_clusters.R`

__Purpose__:
Finds GOBPs that are enriched in each cluster of a disease

__Arguments__:

1. Path to cluster file
2. Background genes from network
3. Output directory

__Run__:
```bash
Rscript find_GOBP_enriched_clusters_GenePlexus.R \
 ../results/prediction_clusters_same_graph/clusters/predicted_withConsensusPathDB--clustered_on_ConsensusPathDB/Malignant_neoplasm_of_pancreas--threshold--0.8--PredictionGraph--ConsensusPathDB--ClusterGraph--ConsensusPathDB_clusters.csv \
 ../data_Zenodo/ConsensusPathDB/ConsensusPathDB_genes.csv  \
 ../results/prediction_clusters_same_graph/GOBP_enrichment

```


### Cluster overlap scores
__Script__:  
`scoreClusterOverlaps_GenePlexus.R`

__Purpose__:
For a disease, outputs a file with the overlap score between all
real and fake trait clusters that have >=5 genes with chronic inflammation genes. 

Also outputs a file with the shared genes between all real and fake trait clusters
with chronic inflammation genes

__Arguments__:

1. Path to folder containing leiden cluster output files
2. Path to chronic inflammation prediction file
3. Path to output directory
4. Disease of interest
5. Chronic inflammation prediction threshold
6. List of genes in network the disease genes were clustered on 

__Run__:
```bash
Rscript scoreClusterOverlaps_GenePlexus.R \
 ../results/prediction_clusters_same_graph/clusters/predicted_withConsensusPathDB--clustered_on_ConsensusPathDB \
 ../results/GenePlexus_output/chronic_inflammation_gene_shot--ConsensusPathDB--Adjacency--GO--predictions.tsv \
 ../results/prediction_clusters_same_graph \
 Chronic_Obstructive_Airway_Disease \
 0.8 \
 ../data_Zenodo/ConsensusPathDB/ConsensusPathDB_genes.csv

```


### Retrieve significant overlaps
__Script__:  
`filterSignificantOverlaps_GenePlexus.R`

__Purpose__:
Filters the FDRs for significant values, outputting the significant clusters and
real and fake cluster assignments

__Arguments__:

1. overlap_results.Rdata location
2. FDR cutoff
3. Output directory
4. Path to gene cluster assignment files 

__Run__:
```bash
Rscript filterSignificantOverlaps_GenePlexus.R \
 ../results/prediction_clusters_same_graph/scores/chronic_inflammatory_response_GO2ALLEGS_thresh=0.8_clusteredOn_STRING-EXP_overlap_results.Rdata \
 .05 \
 ../results/prediction_clusters_same_graph/ \
 ../results/prediction_clusters_same_graph/clusters/predicted_withSTRING-EXP--clustered_on_STRING-EXP

```

### Running SAveRUNNER to compare clusters with STRING-EXP as the interactome
__Script__:  
`prepForClusterSaverunner.R`

__Purpose__:
Sets up files for running SAveRUNNER with cluster genes. This instance is stored in `data_Zenodo/prediction_clusters_same_graph/SAveRUNNER`

__Arguments__:

1. Saverunner input directory
2. path to interactome edgelist
3. path to "final for alex" file with significant clusters
4. path to gene cluster assigments

__Run__:
```bash
Rscript prepForClusterSaverunner.R \
 ../data_Zenodo/prediction_clusters_same_graph/SAveRUNNER/code/input_files \
 ../data_Zenodo/prediction_clusters_same_graph/SAveRUNNER/code/input_files/interactome.txt \
 ../data_Zenodo/prediction_clusters_same_graph/chronic_inflammation_gene_shot_pubs_greater10--predictedWith--ConsensusPathDB--clusteredOn--ConsensusPathDB_final_for_alex.csv \
 ../data_Zenodo/prediction_clusters_same_graph/chronic_inflammation_gene_shot_pubs_greater10--predictedWith--ConsensusPathDB--clusteredOn--ConsensusPathDB_relevant_gene_cluster_assigments.csv

```



### Finding drugs
Drugs were obtained using the SAveRUNNER software, located at https://github.com/sportingCode/SAveRUNNER.
The instance is stored in `data_Zenodo/drugs/SAveRUNNER`



### Analyses and visualizations
The `figures` has R markdown scripts that will recreate figures (including supplemental) in the paper


### Drug data downloads
The relevant files from DrugCentral are included in `data/drugcentral`. They came from a local
PostgreSQL instance of the DrugCentral database, which can be obtained from
DrugCentral.

Relevant files from DGIdb are included in `data/dgidb`.

#### DrugCentral Entrez
__Script__:  
`data/drugcentral/getDrugCentralEntrez.R`

__Purpose__:
This script takes tables from DrugCentral and returns a mapping of Drugs and
Entrez targets for humans. This output is already provided.

__Run__:
```bash
Rscript getDrugCentralEntrez.R
```

