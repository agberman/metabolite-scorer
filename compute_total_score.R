# Adam Berman
# Princeton University
# 5 May 2017


setwd("/Users/adamberman/metabolite-scorer")

library(futile.matrix)
library(reshape2)
library(KEGGREST)
library(plyr)
library(ROCR)


#----------Compute mutational score----------------------------------------------------------------
#----------Compute mutational score----------------------------------------------------------------
#----------Compute mutational score----------------------------------------------------------------

# Read in metabolites and their associated genes
# The 'kegg_metabolite_table.txt' file is created by pipelining the output of 
# make_hmdb_metabolites.py into a text file, where make_hmdb_metabolites.py 
# (as well as kegg_metabolite_table.txt itself) as can be found at: https://github.com/agberman/metabolite-scorer.git
metabolites <- read.table('kegg_metabolite_table.txt', header = TRUE, sep = '\t', quote = "", comment.char = "")

# Read in genes and their gene lengths
# The amino_acid_gene_lengths file can be found in the Github repository at 
# the following link: https://github.com/agberman/metabolite-scorer.git
gene_lengths <- read.table('amino_acid_gene_lengths', header = TRUE, sep = '\t', quote = "", comment.char = "")

# Read in raw TCGA data
# BRCA_US_DATA.csv can be downloaded in the Github repo above, or at the 
# following link: https://dcc.icgc.org/api/v1/download?fn=/current/Projects/BRCA-US/donor.BRCA-US.tsv.gz
samples <- read.csv('BRCA_US_DATA.csv')

# Sum together total number of mutations to each gene
gene_scores <- aggregate(samples$X..mutations.at.gene, by=list(gene_symbol=samples$gene_symbol), FUN=sum)

# Remove empty string from top of gene_scores and rename columns
gene_scores = gene_scores[-1,]
colnames(gene_scores) <- c('gene_symbol', 'total_mutations')

# Match genes to gene lengths
gene_scores$gene_length <- gene_lengths$LENGTH[match(gene_scores$gene_symbol, gene_lengths$Hugo_Symbol)]

# Remove genes with no gene length
gene_scores <- na.omit(gene_scores)

# Compute mutational gene subscores by dividing total mutations by gene length
gene_scores$mutational_gene_subscore <- (gene_scores$total_mutations / gene_scores$gene_length)

# Sort mutational gene subscores
gene_scores_sorted <- gene_scores[with(gene_scores, order(-mutational_gene_subscore)), ]

# Plot a histogram of mutational gene subscores
hist(gene_scores$mutational_gene_subscore, xlab = 'Mutational Gene Subscore', main = 'Distribution of Mutational Gene Subscores')

# Match mutational gene subscores to genes in metabolite list
metabolites$mutational_gene_subscore <- gene_scores$mutational_gene_subscore[match(metabolites$Gene.Name, 
                                                                                    gene_scores$gene_symbol)]

# Remove rows where gene has no mutational gene subscore
metabolites <- na.omit(metabolites)

# Sort metabolites by decreasing mutational gene subscore
#metabolites_sorted <- metabolites[with(metabolites, order(-mutational_gene_subscore)), ]

# Find average of mutational gene subscores associated with each metabolite
metabolites_scored <- aggregate(mutational_gene_subscore~Metabolite, metabolites, FUN=mean)

# CLEAN UP: Rename columns and return KEGG IDs to each metabolite after aggregation, then reorder columns
colnames(metabolites_scored) <- c('metabolite', 'mutational_score')
metabolites_scored$KEGG_ID <- metabolites$KEGG[match(metabolites_scored$metabolite, metabolites$Metabolite)]
metabolites_scored <- metabolites_scored[,c(1,3,2)]

# Normalize mutational scores to be between 0 and 1
metabolites_scored$normalized_mutational_score <- ((metabolites_scored$mutational_score - min(metabolites_scored$mutational_score)) / (max(metabolites_scored$mutational_score) - min(metabolites_scored$mutational_score)))

# Sort metabolites by descending mutational score
metabolites_scored_sorted <- metabolites_scored[with(metabolites_scored, order(-mutational_score)), ]

# Copy results into new data frames
mutational_results <- metabolites_scored
mutational_results_sorted <- metabolites_scored_sorted

# Change column names
colnames(mutational_results) <- c("Metabolite", "KEGG_ID", "Mutation_Score", "Normalized_Score")
colnames(mutational_results_sorted) <- c("Metabolite", "KEGG_ID", "Mutation_Score", "Normalized_Score")

# Generate a histogram of the mutation scores
hist(mutational_results$normalized_mutational_score, xlab = 'Normalized Mutational Score', main = 'Distribution of Normalized Mutational Scores')





#----------Compute expressional score----------------------------------------------------------------
#----------Compute expressional score----------------------------------------------------------------
#----------Compute expressional score----------------------------------------------------------------

# Read raw data into a data frame, storing only the gene name, the expression value,
# and the indicator of whether a sample is from cancer tissue or healthy tissue
# "exp_array.BRCA-US.tsv", the raw expressional data, can be downloaded at the 
# following link: https://dcc.icgc.org/api/v1/download?fn=/current/Projects/BRCA-US/exp_array.BRCA-US.tsv.gz
df = read.delim("exp_array.BRCA-US.tsv", colClasses = c("NULL","NULL","NULL","NULL","NULL","NULL","NULL",NA,NA,
                                                        "NULL","NULL","NULL","NULL","NULL","NULL",NA,"NULL"))

# Convert raw data into a matrix with one sample per row and one gene per column, where the value 
# at a row-column intersection is the expression value of the column's gene in the row's sample
matrix <- acast(df, raw_data_accession ~ gene_id, value.var="normalized_expression_value")

row_names <- rownames(matrix)
abbrev_row_names <- data.frame(strsplit(row_names, "-"))[4,]
abbrev_row_names <- t(abbrev_row_names)
abbrev_row_names <- as.character(as.vector(abbrev_row_names[,1]))


# 0 indicates a cancer sample (row), 1 indicates a healthy sample (row)
bin_row_names <- as.numeric(substring(abbrev_row_names, 1, 1))
matrix_copy <- matrix
matrix_copy <- cbind(matrix_copy, bin_row_names)

# 1. Split matrix into cancer and healthy submatrices
cancer_matrix <- matrix_copy[matrix_copy[,17815] == 0, c(1:17814)]
healthy_matrix <- matrix_copy[matrix_copy[,17815] == 1, c(1:17814)]


# 2. Run t-test between each corresponding pair of rows in the two submatrices

# Define a function to run your t.test, grab the p value, and put them in a data.frame
f <- function(x,y){
  test <- t.test(x, y)
  #out <- data.frame(pval = test$p.value)
  out <- test$p.value
  return(out)
}

# Iterate over columns in cancer_matrix and healthy_matrix via sapply
p_vals <- sapply(seq(ncol(cancer_matrix)), function(x) f(cancer_matrix[,x], healthy_matrix[,x]))

# Adjust p-values for multiple hypothesis testing
p_vals <- p.adjust(p_vals, method = "bonferroni")

# Binarize p-values according to the threshold of 0.01
binarized_p_vals <- as.integer(p_vals < 0.01)

# Histogram of p_vals
hist(p_vals)
hist(binarized_p_vals)

bpv <- data.frame(colnames(matrix), binarized_p_vals)
colnames(bpv) <- c("gene_symbol", "binarized_p_val")

# Import metabolite data
metabolites <- read.table('kegg_metabolite_table.txt', header = TRUE, sep = '\t', quote = "", comment.char = "")

# Append binarized p-value of each gene in metabolites to metabolites
metabolites$binarized_p_val <- bpv$binarized_p_val[match(metabolites$Gene.Name, bpv$gene_symbol)]  

# Remove rows with no binarized p-value
metabolites <- na.omit(metabolites)


# Count how many significantly mutated genes are associated with each metabolite
metabolites_with_bpv_sum <- aggregate(metabolites$binarized_p_val, by=list(metabolites$KEGG, metabolites$Metabolite), FUN=sum)
colnames(metabolites_with_bpv_sum) <- c("KEGG", "Metabolite", "bpv_sum")

# Count the number of genes associated with each metabolite
metabolites_with_alphas <- aggregate(cbind(num_genes = Gene.Name) ~ Metabolite, 
                                      data = metabolites, 
                                      FUN = function(x){NROW(x)})
# nrow(metabolites_with_alphas[which(metabolites_with_alphas$num_genes == 1),])

# Create a metabolite-per-row data frame associating BPV totals with number of genes associated with each metabolite
expressional_scores <- data.frame(metabolites_with_bpv_sum$Metabolite, metabolites_with_bpv_sum$KEGG, metabolites_with_bpv_sum$bpv_sum, metabolites_with_alphas$num_genes)
colnames(expressional_scores) <- c("Metabolite", "KEGG", "bpv_sum", "num_genes")

# Compute expressional scores
expressional_scores$Expression_Score <- expressional_scores$bpv_sum / expressional_scores$num_genes

# Normalize expressional scores to be between 0 and 1
expressional_scores$Normalized_Score <- ((expressional_scores$Expression_Score - min(expressional_scores$Expression_Score)) / (max(expressional_scores$Expression_Score) - min(expressional_scores$Expression_Score)))

# Remove metabolites with fewer than 4 associated genes
expressional_scores_trimmed <- expressional_scores[expressional_scores$num_genes >= 4, ]

# Re-normalize
expressional_scores_trimmed$Normalized_Score <- ((expressional_scores_trimmed$Expression_Score - min(expressional_scores_trimmed$Expression_Score)) / (max(expressional_scores_trimmed$Expression_Score) - min(expressional_scores_trimmed$Expression_Score)))

# Sort metabolites by expressional scores
expressional_scores_sorted <- expressional_scores_trimmed[with(expressional_scores_trimmed, order(-Expression_Score, -bpv_sum)), ]

# Copy results into new data frames
expressional_results <- expressional_scores_trimmed
expressional_results_sorted <- expressional_scores_sorted

# Modify column names
colnames(expressional_results) <- c("Metabolite", "KEGG_ID", "bpv_sum", "num_genes", "Expression_Score", "Normalized_Score")
colnames(expressional_results_sorted) <- c("Metabolite", "KEGG_ID", "bpv_sum", "num_genes", "Expression_Score", "Normalized_Score")

# Generate a histogram of the mutation scores
hist(expressional_results$Normalized_Score, xlab = 'Normalized Expressional Score', main = 'Distribution of Normalized Expressional Scores')





#----------Retrieve KEGG data----------------------------------------------------------------
#----------Retrieve KEGG data----------------------------------------------------------------
#----------Retrieve KEGG data----------------------------------------------------------------

# Read in raw KEGG data: all metabolic pathways found in humans ("hsa")
paths <- names(keggList("pathway", "hsa"))

# Initialize vectors to hold pathways and their constituent 
# metabolites (paths_with_metabolites) and metabolite names (names)
paths_with_metabolites <- vector(mode="list", length=length(paths))
names <- vector(mode="character", length=length(paths))

# Initialize lists to hold metabolites found in all pathways (met_pool)
# and metablites found only in cancer-linked pathways (cancer_pool)
met_pool <- character()
cancer_pool <- character()

# Iterate over all pathways
for (i in 1:length(paths)){
  print(paste("i:", i))
  # Parse out all metabolites in given pathway
  dat <- keggGet(paths[i])
  names[i] <- dat[[1]]["NAME"][[1]]
  metabolites <- names(dat[[1]]["COMPOUND"][[1]])
  
  # Add all metabolites to met_pool
  met_pool <- c(met_pool, metabolites)
  
  # Add metabolites to cancer_pool if pathway is cancer-linked
  if ((i >= 279) && (i <= 301)) {
    cancer_pool <- c(cancer_pool, metabolites)
  }
  
  # Update paths_with_metabolites with new pathway
  paths_with_metabolites[[i]] <- metabolites
}

# Associate names with the stored pathways
pwp <- paths_with_metabolites
names(pwp) <- names

# Remove repeated metabolites in cancer_pool
cancer_pool <- unique(cancer_pool)





#----------Compute total scores at all alphas from 0 to 1 by 0.05----------------------------------------------------------------
#----------Compute total scores at all alphas from 0 to 1 by 0.05----------------------------------------------------------------
#----------Compute total scores at all alphas from 0 to 1 by 0.05----------------------------------------------------------------
#----------total_score = (mutation_score*alpha) + (expression_score*(1-alpha))---------------------------------------------------

# Remove all metabolites without a KEDD ID
mutation_scores <- mutational_results[mutational_results$KEGG_ID != "None", ]
expression_scores <- expressional_results[expressional_results$KEGG_ID != "None", ]

# Remove all metabolites that are not present in both the mutational and expressional results 
# (we trimmed metabolites with fewer than four associated metabolites in the expressional list)
mutation_scores$expressional_score <- expression_scores$Expression_Score[match(mutation_scores$Metabolite, 
                                                                                expression_scores$Metabolite)]
mutation_scores <- na.omit(mutation_scores)
mutation_scores <- mutation_scores[, -5]

# Remove all metabolites not present in KEGG pathways
mutation_scores <- mutation_scores[mutation_scores$KEGG_ID %in% met_pool, ]
expression_scores <- expression_scores[expression_scores$KEGG_ID %in% met_pool, ]

# Remove last three rows of mutation and expression tables (513 rows -> 510 rows)
mutation_scores <- head(mutation_scores, -3)
expression_scores <- head(expression_scores, -3)

# Re-normalize scores
mutation_scores$Normalized_Score <- ((mutation_scores$Mutation_Score - min(mutation_scores$Mutation_Score)) / (max(mutation_scores$Mutation_Score) - min(mutation_scores$Mutation_Score)))
expression_scores$Normalized_Score <- ((expression_scores$Expression_Score - min(expression_scores$Expression_Score)) / (max(expression_scores$Expression_Score) - min(expression_scores$Expression_Score)))

# Verify that mutation_scores and expression_scores align perfectly
all(mutation_scores$Metabolite == expression_scores$Metabolite)

# Create matrix to store total scores
#alphas <- c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1)

# Create sequence of all alphas from 0 to 1 by 0.01
alphas <- seq(from = 0, to = 1, by = 0.01)

# Create matrix for storing overall scores at all alpha values
total_scores <- matrix(data = NA, nrow = nrow(mutation_scores), ncol = length(alphas))
rownames(total_scores) <- mutation_scores$KEGG_ID
colnames(total_scores) <- alphas

# Calculate overall scores based on (mutation_score*alpha) + (expression_score*(1-alpha)) 
# for all alphas 0-1 by 0.01 using normalized scores
for(column in 1:ncol(total_scores)){
  alpha <- as.numeric(colnames(total_scores)[column])
  total_scores[, column] <- ((mutation_scores$Normalized_Score * alpha) 
                              + (expression_scores$Normalized_Score * (1 - alpha)))
}





#----------Perform machine learning to identify best alpha value----------------------------------------------
#----------Perform machine learning to identify best alpha value----------------------------------------------
#----------Perform machine learning to identify best alpha value----------------------------------------------

# Initialize list to store, for each metabolite, a 1 if the metabolite 
# is present in a KEGG cancer-linked pathway, or a 0 if not
in_cancer_pool <- c()

# Iterate through all metabolites in the overall score matrix
for (i in 1:length(rownames(total_scores))) {
  if (rownames(total_scores)[i] %in% cancer_pool) {
    # Mark a 1 in the list if the metabolite appears in a KEGG cancer-linked pathway
    in_cancer_pool = c(in_cancer_pool, 1)
  } else {
    # Mark a 0 in the list if the metabolite does not appear in a KEGG cancer-linked pathway
    in_cancer_pool = c(in_cancer_pool, 0)
  }
}

# Add the list of 1's and 0's for each metabolite to the overall score matrix
total_scores <- cbind(total_scores, in_cancer_pool)


# Randomly permute data row-wise
total_scores_random <- total_scores[sample(nrow(total_scores)),]

# Split random data into 5 training sets and 5 hold out sets
set_size <- (nrow(total_scores_random) / 5)

set1 <- total_scores_random[(set_size+1):nrow(total_scores_random), ]
set1_holdout <- total_scores_random[1:set_size, ]

set2 <- total_scores_random[1:set_size, ]
set2 <- rbind(set2, total_scores_random[((2*set_size)+1):nrow(total_scores_random), ])
set2_holdout <- total_scores_random[(set_size+1):(2*set_size), ]

set3 <- total_scores_random[1:(2*set_size), ]
set3 <- rbind(set3, total_scores_random[((3*set_size)+1):nrow(total_scores_random), ])
set3_holdout <- total_scores_random[((2*set_size)+1):(3*set_size), ]

set4 <- total_scores_random[1:(3*set_size), ]
set4 <- rbind(set4, total_scores_random[((4*set_size)+1):nrow(total_scores_random), ])
set4_holdout <- total_scores_random[((3*set_size)+1):(4*set_size), ]

set5 <- total_scores_random[1:(4*set_size), ]
set5_holdout <- total_scores_random[((4*set_size)+1):nrow(total_scores_random), ]

sets <- list(set1, set2, set3, set4, set5)
holdout_sets <- list(set1_holdout, set2_holdout, set3_holdout, set4_holdout, set5_holdout)

# Run machine learning simulation on all training set-holdout set pairs, 
# computing a maximum alpha and corresponding maximum holdout AUC for each
alpha_maxes <- c()
final_aucs <- c()
i = 0
for (set in sets) {
  
  i = i + 1
  
  auc_vals <- matrix(data = NA, nrow = 1, ncol = length(alphas))
  rownames(auc_vals) <- c("AUC_Values")
  colnames(auc_vals) <- alphas
  
  for(column in 1:(ncol(total_scores_random)-1)) {
    prediction.obj <- prediction(set[, column], set[, "in_cancer_pool"])
    perf <- performance(prediction.obj, measure = 'auc')
    auc <- perf@y.values[[1]]
    auc_vals[1, column] <- auc
  }
  
  alpha_max <- colnames(auc_vals)[apply(auc_vals,1,which.max)]
  alpha_maxes <- c(alpha_maxes, as.numeric(alpha_max))
  
  holdout_pred <- prediction(holdout_sets[[i]][, alpha_max], holdout_sets[[i]][, "in_cancer_pool"])
  holdout_perf <- performance(holdout_pred, measure = 'auc')
  holdout_auc <- holdout_perf@y.values[[1]]
  final_aucs <- c(final_aucs, holdout_auc)
}

# Set final alphas to be those found in my writeup (due to the random nature 
# of the set divisions, alpha values will vary each time the code is run)
# final_alpha <- mean(c(1.0, 1.0, 0.13, 1.0, 1.0))

# Compute average of alpha maxes (i.e. the optimal alpha)
final_alpha <- mean(alpha_maxes)

# Compute the average of the final AUC values
final_auc <- mean(final_aucs)

# Create a box-and-whiskers plot of the final AUC values
boxplot(final_aucs, main="Final Area Under ROC Curve (AUC) Values")





#----------Compute total scores at best alpha for all metabolites (even those without a KEGG ID)---------------------------------
#----------Compute total scores at best alpha for all metabolites (even those without a KEGG ID)---------------------------------
#----------Compute total scores at best alpha for all metabolites (even those without a KEGG ID)---------------------------------
#----------total_score = (mutation_score*alpha) + (expression_score*(1-alpha))---------------------------------------------------

# Remove all metabolites that are not present in both the mutational and expressional results 
# (we trimmed metabolites with fewer than four associated metabolites in the expressional list)
mutational_results$expressional_score <- expressional_results$Expression_Score[match(mutational_results$Metabolite, expressional_results$Metabolite)]
mutational_results_trimmed <- na.omit(mutational_results)
mutational_results_trimmed <- mutational_results_trimmed[, -5]

# Raw (unweighted) sum of normalized mutation and expression scores
total_results <- data.frame(Metabolite = mutational_results_trimmed$Metabolite, 
                            KEGG_ID = mutational_results_trimmed$KEGG_ID, 
                            Normalized_Mutation_Score = mutational_results_trimmed$Normalized_Score, 
                            Normalized_Expression_Score = expressional_results$Normalized_Score, 
                            Total_Score = (((mutational_results_trimmed$Normalized_Score)*final_alpha) + ((expressional_results$Normalized_Score)*(1 - final_alpha))))

# Compute normalized total score
total_results$Normalized_Total_Score <- ((total_results$Total_Score - min(total_results$Total_Score)) / (max(total_results$Total_Score) - min(total_results$Total_Score)))

# Sort by total scores
total_results_sorted <- total_results[with(total_results, order(-Total_Score)), ]

# Print the 20 metabolites with highest normalized 
# overall scores at the optimal alpha value
head(total_results_sorted, 20)

# Generate a histogram of normalized total scores, computed using the 
# best alpha value as computed by ROC machine learning techniques
hist(total_results$Normalized_Total_Score, xlab = 'Normalized Overall Score', main = 'Distribution of NormalizedOverall Scores at Alpha=0.01')





#----------Compute final p-values using Hypergeometric test----------------------------------------------------------------
#----------Compute final p-values using Hypergeometric test----------------------------------------------------------------
#----------Compute final p-values using Hypergeometric test----------------------------------------------------------------
# This test is optional

'''
# Random sampling function
sample_pval <- function(sample_size, threshold, num_samples, met_pool, cancer_pool) {
  samples <- replicate(num_samples, sample(met_pool, sample_size, replace = TRUE))
  num_ge_threshold = 0
  for (i in 1:num_samples) {
    ge_threshold_count <- length(intersect(cancer_pool, samples[,i]))
    if (ge_threshold_count >= threshold) {
      num_ge_threshold <- num_ge_threshold + 1
    }
  }
  return(num_ge_threshold/num_samples)
}

# Define quantiles to be tested
quantiles <- c(1, 5, 10, 15, 20)
#alphas <- as.numeric(colnames(total_scores))
alphas <- seq(from = 0, to = 1, by = 0.01)

# Initialize data structure to hold all final p-values
final_pvals <- matrix(data = NA, nrow = length(alphas), ncol = length(quantiles))
rownames(final_pvals) <- alphas
colnames(final_pvals) <- quantiles

# Calculate all p-values (one per alpha per quantile)
for (a in 1:length(alphas)) {
  
  alph <- total_scores[,a]
  
  for (q in 1:length(quantiles)) {
    
    # Get metabolites with score in desired quantile
    cuttoff <- quantile(alph, prob=1-quantiles[q]/100)
    in_quantile <- names(alph[alph >= cuttoff])
    sample_size <- length(in_quantile)
    
    # Find out how many metabolites in quantile are in KEGG cancer pool
    threshold <- length(intersect(cancer_pool, in_quantile))
    
    # Calculate p-value via sampling and add it to final table
    num_samples <- 1000
    final_pvals[a,q] <- sample_pval(sample_size, threshold, num_samples, met_pool, cancer_pool)
    
  }
}

final_pvals_df <- data.frame(final_pvals)
final_pvals_df$avg <- rowMeans(final_pvals_df)
final_pvals_df[final_pvals_df$avg == min(final_pvals_df$avg)]
'''
