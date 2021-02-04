# Importing and tidying the data -----------------------------------------------

# path variable where imported counts.txt would be found.
# file.path(getwd() function is way to make code transferable to different
# operating systems.
path <- file.path(getwd(), "data", "counts.txt")

# Actual import 
# read.table() creates dataframe of file
# sep = , and escape sequence, to insert horizontal tab when tab key is pressed. 
# header = is used if the first line of the file should be used for column names. 
counts_pen <- read.table(file = path, sep = '\t', header = FALSE)

# install.packages("tidyverse")
library(tidyverse)

# Subtracted the first four lines/rows of counts
# tail(x,n=) generally returns last n rows of a data frame
counts_pen <- tail(counts_pen, -4)

# Created a sequence of 4-156 with intervals of 4 (kept only every third column)
# [ , is for the rows, and c(1,) keeps the first column
counts_pen <- counts_pen[ , c(1, seq(4, 24, 4))]

# Assigning row_names to be the first element in counts
row_names <- counts_pen[1]
# str_sub: Extract and replace substrings from a character vector 
# (This function starts with the name of V1 at the sixth letter and end in the 
# last one).
# The $ helps extract v1.
row_names <- str_sub(row_names$V1) #tidyverse
# Making the row names of counts to be equivalent to row_names above
rownames(counts_pen) <- row_names

# Took out column v1
counts_pen <- counts_pen[ , 2:7 ]

# Assigning col_names the values of the SRR numbers
col_names<- c("SRR3119155", "SRR3119156",
              "SRR3119162", "SRR3119163",
              "SRR3119169", "SRR3119170")
# Now col_names will be the column names in counts
colnames(counts_pen) <- col_names

# This is to import the metadata.tsv file
metadata <- read.table(file = "./data/metadata.tsv",
                       sep = '\t',
                       header = TRUE)

# A subset of the metadata that extracts only the 6 replicates I'm looking for
subset_meta <- metadata[metadata$name == "pen_S-" |
                          metadata$name == "pen_S+",]

# colnames() assigns the column names as the values below
colnames(counts_pen) <- c("S_pen_before_flower_rep_1",
                          "S_pen_after_flower_rep_1",
                          "S_pen_before_flower_rep_2",
                          "S_pen_after_flower_rep_2",
                          "S_pen_before_flower_rep_3",
                          "S_pen_after_flower_rep_3")

# coldata is a file with the column names as the row names from the counts file
# rep() = Replicates values
coldata_pen <- data.frame(row.names = colnames(counts_pen),
                          condition = factor(rep(c("pre_flower",
                                                   "post_flower"), 3)))
# To make sure counts are in order
str(counts_pen)

# factor() is converting the character coldata_pen$condition into a factor
# R uses reference level for factors based on alphabetical order
# factors save space for storing categorical data by converting characters into 
# numbers
coldata_pen$condition <- factor(coldata_pen$condition, levels = c("pre_flower",
                                                                  "post_flower"))
    
# DESeq2------------------------------------------------------------------------ 

# To download DESeq2 because the package was not found before
# if (!requireNamespace("BiocManager", quietly = TRUE))#R code missing something
# install.packages("BiocManager")

# BiocManager::install("DESeq2")

library("DESeq2")
# browseVignettes("DESeq2")

# Making the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_pen,
                             colData = coldata_pen,
                             design = ~ condition)

# Calculating differential expression-------------------------------------------

# Analysis steps wrapped in DESeq function
# ?DESeq
dds <- DESeq(dds)
res <- results(dds)

# summary() produces results summaries, in this case cooksCutoff threshold and
# independentFiltering from DESeq2 was used
summary(res)
#?results

# To view result names to see correct format of coef
resultsNames(dds)

# To install apeglm because error in code below said to install
# BiocManager::install("apeglm")

# Shrinking LFC estimates for visualization and ranking of genes
resLFC <- lfcShrink(dds,
                    coef = "condition_post_flower_vs_pre_flower",
                    type = "apeglm")

# Vizualization-----------------------------------------------------------------

# MA-plot
# This plot shows log2 fold change attributable to a given variable over the 
# mean of normalized counts for all the samples in the DESeqDataSet.  
plotMA(res, ylim=c(-2,2))
# MA plot of adjusted log fold changes by shrinking (adjusting genes with low 
# numbers to be more accurate)
plotMA(resLFC, ylim=c(-2,2))

# PCA
# First transform data with variance stabilizing transformation to remove
# dependence of variance on the mean. Useful for visualization and downstream
# analysis. (remove returnData = TRUE for basic PCA)
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# ggplot2 is based on the grammar of graphics, ggplot() creates a coordinate 
# system that you can add layers to.
ggplot(pcaData, aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(aes(color = condition, shape = condition), size = 5) +
  labs(title = "Tomato pre/post flower PCA",
       x = paste0("PC1: ",percentVar[1],"% variance"),
       y = paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c("#007282", "#009E73")) +
  theme_bw() +
  theme(axis.title = element_text(size = 26, face = 'bold'),
        axis.text = element_text(size = 22, face = 'bold', color = 'black'),
        axis.text.x = element_text(size = 20, face = 'bold', color = 'black'),
        plot.title = element_text(size = 20, face = 'bold', margin = margin(0,0,10,0)),
        axis.title.x = element_text(margin = margin(10,0,0,0)),
        panel.border = element_blank(),
        axis.line = element_line(size = 1, color = 'black'),
        axis.ticks = element_line(size = 1, color = 'black'),
        axis.ticks.length = unit(8, 'pt'),
        plot.margin = margin(0.5, 0.5, 0.5,0.5, 'cm'),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 20, face = 'bold'))

ggsave(filename = "PCA_S_pen_pre_post_flower.png",
       device = 'png',
       width = 9,
       height = 7.5,
       dpi = 400,
       units = 'in')

# Turning results object into a dataframe
res_df <- as.data.frame(res)

# Filtering to only counts that passed cutoff, took out NA's
res_df <- res_df[complete.cases(res_df), ]

# Ordering results from high Log2Fold values to low
up_expressed_by_lfc <- res_df[order(res_df$log2FoldChange, decreasing = T), ]
# Only values less than adjusted p-value 0.05 kept
up_expressed_by_lfc <- up_expressed_by_lfc[up_expressed_by_lfc$padj < 0.05, ]
# Orders results from low to high Log2Fold values
down_expressed_by_lfc <- res_df[order(res_df$log2FoldChange), ]
# Only values less than adjusted p-value 0.05 kept
down_expressed_by_lfc <- down_expressed_by_lfc[down_expressed_by_lfc$padj < 0.05, ]

# Top 50 counts with greatest Log2Fold values and significant padj
top_pen <- up_expressed_by_lfc[1:50, ]

# Top 50 counts with most negative Log2Fold values and significant padj
bottom_pen <- down_expressed_by_lfc[1:50, ]

# To find top 50 values
#top_pen <- top_n(up_expressed_by_lfc, 50)

# Bottom 50 values
#bottom_pen <- top_n(down_expressed_by_lfc, -50)

# Predicted functions----------------------------------------------------------

# Sol Genomics annotations
path <- file.path(getwd(), "data", "spenn_v2.0_gene_models_annot.gff")
annotations <- read.table(file = path, sep = '\t', header = FALSE)

# To focus on 3rd and ninth columns
annotations <- annotations[ , c(3,9)]

# Extracted all of the rows with mRNA in it
annotations <- annotations %>%
  group_by(V9) %>%
  filter(any(V3 == "mRNA"))

# To check for duplicates
duplicated(annotations)
# To get rid of duplicates if there were any
annotations <- unique(annotations)

# To separate V9 into more columns
annotations <- annotations %>% separate(V9,
                         into = c("ID", "Name","Parent", "Note"),
                         sep = ";")

# Extracting parent and note
annotations <- annotations[ , c(4,5)]

# To get rid of Note=
annotations$Note <- substring(annotations$Note, 6)

# To get rid of Parent=
annotations$Parent <- substring(annotations$Parent, 8)

# Making row names into a column named Parent
res_df <- tibble::rownames_to_column(res_df, "Parent")
top_pen <- tibble::rownames_to_column(top_pen, "Parent")
bottom_pen <- tibble::rownames_to_column(bottom_pen, "Parent")

# To add functional annotations
res_functions <- inner_join(res_df, annotations, by = "Parent")
top_pen <- inner_join(top_pen, annotations, by = "Parent")
bottom_pen <- inner_join(bottom_pen, annotations, by = "Parent")
