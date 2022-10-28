# Network_pleiotropy
This repo contains all scripts used for analysis in "Network topology can explain differences in the pleiotropy of cis- and trans-regulatory mutations."

The "Deletiondata.R" script contains data cleaning of the Kemmeren et al (2014) dataset, as well as construction of the network and analysis of cis- and trans-acting pleiotropy. 

The "fitness.R" script contains analysis of the MacLean et al (2017) dataset for comparisons of cis- and trans-acting deletion mutants fitness.

The "paramsweep.R" script contains analysis of the cis and trans pleiotropy comparisons using different significance and fold change cutoffs to construct the network from the Kemmeren dataset.

The "permutationtests.R" script contains the code used to permute edges in the network, both maintaining and destroying the network degree distribution.

The "Effectsizecis.R" script contains the code used to compare effect sizes between cis- and trans-acting deletion mutants on the expression of the focal gene. 
