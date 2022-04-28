`0.0_sampleNameMapping.csv`: Sample name mapping, used to map sample ID for Kallisto-related analysis.

`0_e*.txt`: Marker genes from previously published review papers, separated by processes.

`20210707_stringDB.cys`: Cytoscape session used to generate STRING subnetworks.
`20220324_GENIE3.cys`: Cytoscape session used to generate GENIE3 subnetworks.

`Mus_musculus_TF_cofactors.txt`: mouse transcription factor (TF) and co-TF obtained from AnimalTF database.
`mm10_TF_coTF_EnsemblID.txt`: mouse transcription factor (TF) and co-TF obtained from AnimalTF database, with gene stable ID.
`Mus_musculus_grcm38_coding_transcripts.txt`: mouse protein coding transcript, obtained from Ensembl Sep. 2019.

`TPM_allSamples.tsv`: TPM counts on transcript level for all samples, including outliers.
`estCounts_allSamples.tsv`: raw counts on transcript level for all samples, including outliers.
`e*.5geneLevelTPM.rda`: TPM counts on gene levels, separated by timepoints.

`combine-test-expression1.Rdata`: background data for PlacentaCellEnrich analysis.

`e*.5Group.txt`: timepoint hierarchical cluster groups.
`e*.5_e*.5vsE*.5_DEtransGenes.txt`: DE transcripts that have DE genes from pairwise-timepoint DE analysis.

`e*.5specific_TFensGenes.txt`: timepoint-specific transcription factor genes.
`e*.5specific_TSS.bed`: timepoint-specific genes' promoters, defined as 50bp up and downstream of the gene. For this publication, we only used the 4th column which corresponds to gene stable ID.
`e*.5specific_ensGenes.txt`: timepoint-specific genes with gene stable ID.
`e*.5specific_trans.txt`: timepoint-specific transcripts.

`tpmForClustering.txt`: TPM file used for clustering validation.

`trans_cluster_*.txt`: clustering groups from three different algorithms.

`lo.rda`: data object obtained from deconvolution analysis.

GENIE3/: directory containing files related to GENIE3 networks
- `GENIE3/e*/e*.5_0.9.txt`: original GENIE3 networks before subclustering.
- `GENIE3/e*/e*.5_*/e*.5_*BP_filtered.rda`: gene ontology results of the subnetwork.
- `GENIE3/e*/e*.5_*/e*.5_*_TSS.bed`: promoter regions of genes in subnetwork,
- `GENIE3/e*/e*.5_*/e*.5_*_hubGenes.rda`: hub genes of the subnetwork.
- `GENIE3/e*/e*.5_*/e*.5_*_nodeTable.csv`: network node tables with related metrics.

STRING/: directory containing files related to STRING networks
- `STRING/e*.5/e*.5_string_interactions.tsv`: original networks before subclustering.
- `STRING/e*.5/e*.5_string_mapping.tsv`: mapping from gene names to protein names, obtained from STRINGdb ver 11.0b
- `STRING/e*.5/largestComponent/e*.5_largestComponent_nodes.csv`: largest connected component from original networks.
- `STRING/e*.5/largestComponent/e*.5_*/`: directory containing files related to subnetworks. Structures are similar to these of GENIE3.

PlacentaOntology/: directory with PlacentaOntology results, obtained from Webgestal.
