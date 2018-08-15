# Research-project

#Explanatory file: 

#####SCRIPTS: 

#Visualising the data folder:
- gene_comparison_new.R = script for function to plot the distribution of 1/2/3 genes relative to all counts. Makes plots comparing histogram for all counts and 1/2/3 genes.
- 10percent_givengene.R = script for function that when given a gene name or names will identify the experiments with the top and bottom 10% of counts and make a wordcloud of the metadata of those experiments. 


#Transcripts folder:
- Transcripts.tsv = Gene transcripts downloaded from FlyMine 
- Transcripts.R = Code for making transcript exploratory histograms, calculating the average transcript length for every gene and making the log2 transcript length fold-change plot
- mean_gene_lengths.csv = contains all the mean gene lengths used. 
- Length_foldchange_transcripts.csv = transcript table but with logFC values added on 


#Ranking the data folder:
- ranking_FUNCTIONS_TPM.R = Most of the other scripts depend on this script being loaded prior using source(). Script loads in the TPM normalised counts, adds a mask as 1st column, has the function for generating a per-gene threshold and binarising the data (setthresh). Function also for global threshold (globalthresh), calculating the hamming distance and generating the ranked list (hamming), calculating hamming distance but don't return full binary sequence attached = quicker (ham2), multiple_set function (think this is obsolete now and no script uses it), mode function (how to assign a 1 or 0 when forming a consensus sequence), Majority function = creates consensus sequence. 

##Single query folder (within ranking the data):
- TPM_rankonce_rep.R = script for selecting 5 random genes and flipping locations and recording change in rank using single query and plotting 
- TPM_rankOnce_rep.Rdata = R object containing result from script above
- Assessment_singlequery.R = Using a gene from each of the gene test lists as a single query and finding ranks of other genes and plotting single query test list box plots
- Assess_single_query_testlists.Rdata = Results from above script

##Multiple query folder (within ranking the data):
- TPM_multiple_rankonce_same.Rdata = Script for making 2:11 copies of a gene and mutating and compressing into consensus and plotting
- TPM_multiple_rankOnce.Rdata = Results from above script
- Generate_random_numbers.R = Generating random numbers to be used for mutating random genes for both multiple query methods
- random_numbers.Rdata = random numbers generated from above script 
- Plotting_pergene_random_2to11.R = Plotting the results of mutating 2:11 random genes when using the per-gene method
- pergene_comparing_2to11_random.R = script for taking 2:11 random genes, creating a ranked list using the per-gene method and mutating and recording the new rank
- pergene_comparing_2to11_random_doPar.R = Same as above but using in parallel for HPC 
- pergene_2to11gene_random_dopar.Rdata = Results for above script 
- pergene_2to11_close.R = Script for taking 2:11 related genes (smallest hamming distance genes), creating a ranked list using the per-gene method and mutating and recording the new rank 
- pergene_2to11genes_close_dopar.Rdata = Results for above script
- mulq_genes_FUNCTIONS.R = script with functions for the multiple query scripts. Contains functions consensus_rank (create ranked list of genes using consensus method), pergene_rank (create ranked list of genes using pergene method), runsum_rank (runsum algorithm for pergene method), mutate_and_rank (gradually mutate query and find new ranking)
- Consensus_comparing_2to11_random_PG.R = script for taking 2:11 random genes, creating a ranked list using the consensus method and mutating and recording the new rank
- Consensus_2to11genes_random.Rdata = results for above script
- Consensus_comparing_2to11_close.R = Script for taking 2:11 related genes (smallest hamming distance genes), creating a ranked list using the consensus method and mutating and recording the new rank
- Consensus_2to11genes_close.Rdata = Results for above script

##Assessment of multiple query (within ranking the data): 
- Test_fly_lists folder = the test gene lists used (downloaded from FlyMine)
- Plotting_lists.R = Plotting the median, 25% quantile, 75% quantile and their standard deviations for each group, plot boxplots comparing the two multiple query methods 
- All the other .R files and .Rdata files are doing training/test 90/10 splitting 10 times for each of the lists and the results  


#Normalising folder:
	- TPM_normalised_look.R = script for making a histogram of the TPM normalised counts 
	- Making_log2FC_plot.R = script for making the log2 fold-change transcript length plot 
	- Plot_summed_data.R = script for plotting sum of TPM binarised counts and sum of discretised normalised counts
	- TPM_normalised.py = python script for normalising the counts - function is taken from Justin's Github


#Metadata folder:
	- WordCloud.R = R script for making the word clouds 


#Data:
	- Nonlogged_tpm_normalised_counts.csv = TPM normalised counts (not logged) 
	- LoadingData.R = exploring the data, making some plots of the metadata columns 
	- Intermine.R = Looking at using the intermine R package
