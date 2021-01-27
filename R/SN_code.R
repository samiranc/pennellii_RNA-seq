# Importing and tidying the data -----------------------------------------------
# Sets a file path that can work on all operating systems

# path variable where imported counts.txt would be found.
# file.path(getwd() function is way to make code transferable to different
# operating systems.
path <- file.path(getwd(), "data", "counts.txt")

# Actual import 
# ?read.table()
# sep, and escape sequence, to insert horizontal tab when tab key is pressed. 
# header function is used if the first line of the file should be used for
# column names. 
counts_pen <- read.table(file = path, sep = '\t', header = FALSE)

# install.packages("tidyverse")
library(tidyverse)
#####
# Subtracted the first four lines/rows of counts
# tail() returns first of last parts, in this case it returned counts minus 
# first four rows?
counts_pen <- tail(counts_pen, -4)

# From 156 variable to 40? Don't understand what the comma is for
# Created a sequence of 4-156 with intervals of 4 (kept only every third column)
# The c(1,) kept the first column
# Without [ , ] error occurs of undefined columns selected
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

# Assigns subset_meta with unpollinated pistil rows from metadata set
# Does it create a metadata, or what is one?
subset_meta <- metadata[metadata$name == "pen_S-" |
                          metadata$name == "pen_S+",]


# colnames() assigns the column names as the values below
colnames(counts_pen) <- c("S_pen_before_flower_rep_1",
                      "S_pen_before_flower_rep_2",
                      "S_pen_before_flower_rep_3",
                      "S_pen_after_flower_rep_1",
                      "S_pen_after_flower_rep_2",
                      "S_pen_after_flower_rep_3")

# Probably would not need the ones below
# coldata is s file with the column names as the row names from the counts file
# rep() = repetitions? ( Replicates values)
coldata_pen <- data.frame(row.names = colnames(counts_pen),
                          condition = factor(c(rep("pre_flower", 3),
                                              rep("post_flower", 3))))
# factor() is converting the character coldata_pen$condition into a factor
# R uses reference level for factors based on alphabetical order
# factors save space for storing categorical data by converting characters into 
# numbers
coldata_pen$condition <- factor(coldata_pen$condition, levels = c("pre_flower",
                                                                  "post_flower"))
    
# DESeq2------------------------------------------------------------------------                                                                  
# To download DESeq2 because the package was not found before
# if (!requireNamespace("BiocManager", quietly = TRUE))#R code missing something
install.packages("BiocManager")

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

res
summary(res)

# To view result names to see correct format of coef
resultsNames(dds)

# To install apeglm because error in code below said to install
# BiocManager::install("apeglm")

# Shrinking LFC estimates for visualization and ranking of genes
resLFC <- lfcShrink(dds,
                    coef = "condition_pre_flower_vs_post_flower",
                    type = "apeglm")

# Vizualization-----------------------------------------------------------------
# MA-plot
# This plot shows log2 fold change attributable to a given variable over the 
# mean of normalized counts for all the samples in the DESeqDataSet.  
plotMA(res, ylim=c(-2,2))
# Removes low count genes
plotMA(resLFC, ylim=c(-2,2))

# PCA
# First transform data with variance stabilizing transformation to remove
# dependence of variance on the mean. Useful for visualization and downstream
# analysis. (remove returnData = TRUE for basic PCA)?
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# ggplot2 is based on the grammar of graphics
ggplot(pcaData, aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(aes(color = condition, shape = condition), size = 5) +
  labs(title = "Tomato pre/post flower PCA",
       x = paste0("PC1: ",percentVar[1],"% variance"),
       y = paste0("PC1: ",percentVar[2],"% variance")) +
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

# To specify path of the file, not sur if this is what made it work but it did.
path <- file.path(getwd(), "plots", "PCA_S_pen_pre_post_flower.png")

ggsave(filename = "PCA_S_pen_pre_post_flower.png",
       device = 'png',
       width = 9,
       height = 7.5,
       dpi = 400,
       units = 'in')

# Turning results object into a dataframe
res_df <- as.data.frame(res)

# Filtering to only counts that passed cutoff
res_df <- res_df[complete.cases(res_df), ]

# Making some sorted data frames for significant adjusted p-values at 0.05
# Ordering values of res_df from greatest to least.
# ‘$’ refers to a specific column relative to a specific data frame
up_expressed_by_lfc <- res_df[order(res_df$log2FoldChange, decreasing = T), ]
# Only values less than p-value 0.05 kept
up_expressed_by_lfc <- up_expressed_by_lfc[up_expressed_by_lfc$padj < 0.05, ]
# Orders new variable in descending order in log2FoldChange column
down_expressed_by_lfc <- res_df[order(res_df$log2FoldChange), ]
# New variable only has log changes below 0.05
down_expressed_by_lfc <- down_expressed_by_lfc[down_expressed_by_lfc$padj < 0.05, ]

# Predicted functions----------------------------------------------------------
# Sol Genomics annotations
path <- file.path(getwd(), "data", "spenn_v2.0_gene_models_annot.gff")
annotations <- read.table(file = path, sep = '\t', header = FALSE)

# To focus on 3rd and ninth columns
annotations <- annotations[ , c(3,9)]

# Heck yeah, this extracted all of the rows with mRNA in it! 
library(dplyr)

annotations <- annotations %>%
  group_by(V9) %>%
  filter(any(V3 == "mRNA"))

# To separate V9 into more columns
# Did I do this right? Always gives error of discarded additional pieces
library(tidyr)
annotations <- annotations %>% separate(V9,
                         into = c("ID", "Name","Parent", "Note"),
                         sep = ";")

# Extracting parent and note
annotations <- annotations[ , c(4,5)]

# To get rid of duplicates, it does reduce the number of rows but I don't know 
# if its necessary?
#unique(annotations)
#annotations <- annotations %>%
#  distinct(Parent, Note,
#           .keep_all = TRUE)

# To get rid of Note=
annotations$Note <- substring(annotations$Note, 6)

# To make parent names as rownames like res_df
library(tidyverse)
row_names <- annotations[1]
row_names <- str_sub(row_names$Parent, 8, -1)
annotations <- annotations[ , 2]
# Had to use make.name() because (rownames(annotations) <- row_names) kept
# getting an error saying there were duplicates despite the unique() function 
# above.
rownames(annotations) <- make.names(row_names, unique = TRUE)

# Should I make annotations a factor?
# annotations <- factor(annotations) # Don't do this
# class(annotations)


######### Disregard code below. Those were just my failed attempts######### 
############# that I would like to visit later######################
# Want to take out Note=, but line 230 missing something?
# library(stringr)
# Took out the rownames??
# annotations$Note <- substring(annotations$Note, 6)

#annotations <-substring(annotations$Note, 6)

#library(tidyverse)
#Notes <- annotations[1]

#Notes <- str_sub(Notes$Note, 6, -1)

#Notes <- str_sub(annotations[[Notes]])

#annotations[ , Notes]
#str_replace_all(annotations$Note, Notes)
#annotations <- str_sub(annotations, Notes)

######### Disregard code below. Those were just my failed attempts for 
# extracting mRNA rows that I would like to visit later

#annotations %>% slice(mRNA)
#annotations[mRNA.*, ]
#subset(annotations, 'mRNA.*')


#annotations <- data.frame(mRNA = c(NA, "Parent", "Note"))
#annotations %>% separate(mRNA, c("Parent", "Note"))

# str_remove(annotations, "Note=")
#annotations[-grep("CDS", annotations$V3),]

#library(stringr)
#annotations <- annotations<-filter(str_detect(V3, 'mRNA'))

q()

