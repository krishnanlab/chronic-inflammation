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
`figures` contains Markdown notebooks used to analyze final results. Output
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

1. A one column text file of disease(s), CUIDs, of interest
2. Output directory where `disease_gene_files` directory and results will go

__Run__:
```bash
Rscript prep_disease_gene_dfs.R \
 ../data/chronic_inflammation_diseases_non-ovlp_cuid.txt \
 ../data
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
3. location of `disease_gene_files`, where files will be outputted

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
used `STRING`,`Adjacency`, and `DisGeNet` for its final results

__Arguments__:  

-i : Disease seed genes \
-j : Job name \
-n : Network, options are BioGRID, STRING-EXP, STRING, GIANT-TN \
-f : Features, options are Embedding, Adjacency, Influence \
-g : GSC type, options are GO or DisGeNet \
-s : Output directory \
-fl : Option for how to run, always use `local` for this project \

__Run__:
```bash
python example_run.py \
 -i ../../data/disease_gene_files/Chronic_Obstructive_Airway_Disease.txt \
 -j Chronic_Obstructive_Airway_Disease--STRING--Adjacency--DisGeNet \
 -n STRING \
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
  1
```


### Clustering Geneplexus results
__Script__:  
`filterAndClusterGeneplexusPredictions.R`

__Purpose__:
Takes the GenePlexus predictions and assign genes to clusters
for each disease.

__Arguments__:

1. GenePlexus predictions
2. Prediction threshold, either `mccf1` or a number < 1
3. Path to igraph object containing network for clustering
4. Leiden algorithm partition type
5. Resolution parameter
6. GenePlexus results path, creates folder `clusters_threshold.80`

__Run__:
```bash
Rscript filterAndClusterGeneplexusPredications.R \
 ../results/GenePlexus_output/predictions_cv_greater1/Chronic_Obstructive_Airway_Disease--STRING--Adjacency--DisGeNet--predictions.tsv \
 0.8 \
 ../data_Zenodo/biogrid/BioGrid_igraph.Rdata \
 ModularityVertexPartition \
 0.1 \
 ../results/GenePlexus_output/

```


### Clustering random genes
__Script__:  
`clusterRandomGenes.R`

__Purpose__:
Takes the 5000 fake traits that were randomly generated from a
disease and assigns the genes to clusters

__Arguments__:

1. Path to dat containing all fake traits generated
2. Disease of interest
3. Path to igraph object containing network for clustering
4. Partition type
5. Resolution parameter
6. Results path

__Run__:
```bash
Rscript clusterRandomGenes.R \
 ../data_Zenodo/5000Expandedfaketraits_biogrid.tsv \
 Chronic_Obstructive_Airway_Disease \
 ../data_Zenodo/biogrid/BioGrid_igraph.Rdata \
 ModularityVertexPartition \
 0.1
 ../results/GenePlexus_output/clusters_threshold.80

```


### Get GOBP enriched clusters
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
 ../results/GenePlexus_output/clusters_threshold.80/chronic_inflammation_go--threshold--0.8--ClusterGraph--BioGrid_clusters.csv \
 ../data_Zenodo/biogrid/BioGrid_genes.csv \
 ../results/GenePlexus_output/GOBP_enrichment

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
 ../results/GenePlexus_output/clusters_threshold.80 \
 ../results/GenePlexus_output/predictions_cv_greater1/chronic_inflammation_go--STRING--Adjacency--GO--predictions.tsv \
 ../results/GenePlexus_String_Adjacency \
 Chronic_Obstructive_Airway_Disease \
 0.8 \
 ../data_Zenodo/biogrid/BioGrid_genes.csv

```

### FDR calculation
__Script__:  
`calculatePermutedFDR.R`

__Purpose__:
Calculates the FDR for a disease based on how the real clusters overlap with 
chronic inflammation genes relative to the 5000 fake trait clusters.

__Arguments__:

1. Path to folder with overlap score files from `scoreClusterOverlaps_GenePlexus.R`
2. Output directory
3. Boolean. Make `TRUE` if want FDR for only the real diseases/traits

__Run__:
```bash
Rscript calculatePermutedFDR.R \
 ../results/GenePlexus_String_Adjacency/scores \
 ../results/ \
 TRUE
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
4. Path to num_random_clustered folder (created in clusterRandomGenes.R)
5. Path to gene cluster assignment files 

__Run__:
```bash
Rscript filterSignificantOverlaps_GenePlexus.R \
 ../results/overlap_results_real_only.Rdata \
 .01 \
 ../results/ \
 ../results/num_random_clustered \
 ../results/GenePlexus_output/clusters_threshold.80

```


### Finding drugs
Drugs were obtained using the SAveRUNNER software, located at https://github.com/sportingCode/SAveRUNNER.



### Analyses and visualizations
`drug_repurposing_analysis.html` and `figures_and_tables_geneplexus0.80.html` can be recreated with the R notebooks in `figures`


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

