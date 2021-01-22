# Importing and tidying the data ----------------------
# Sets a file path that should work on all operating systems, [surely]

# path variable where imported counts.txt would be found
# file.path(getwd() function is way to make code transferable to different operating systems

path <- file.path(getwd(), "data", "counts.txt")

# Actual import 
# ?read.table()
# sep, and escape sequence, to insert horizontal tab when tab key is pressed 
# header function is used if the first line of the file should be used for column names 

counts_pen <- read.table(file = path, sep = '\t', header = FALSE)

# install.packages("tidyverse")
library(tidyverse)
#####
# Subtracted the first four lines/rows of counts
# tail() returns first of last parts, in this case it returned counts minus first four rows?
counts_pen <- tail(counts_pen, -4)

# From 156 variable to 40? Don't understand what the comma is for
# Created a sequence of 4-156 with intervals of 4 (kept only every third column)
# The c(1,) kept the first column
# Without [ , ] error occurs of undefined columns selected
counts_pen <- counts_pen[ , c(1, seq(4, 24, 4))]

# Assigning row_names to be the first element in counts
row_names <- counts_pen[1]
# str_sub: Extract and replace substrings from a character vector (This function starts with the name of V1 at the sixth letter and end in the last one)
# The $ helps extract v1
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
# Now there are only 6 columns of two samples, 3 reps, of pennellii
counts_pen <- counts_pen %>% select(c(subset_meta$run)) #tidyverse
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
# factors and strings very similar?
coldata_pen$condition <- factor(coldata_pen$condition, levels = c("pre_flower",
                                                         "post_flower"))

# Making the DESeqDataSet object
dds <-DESeqDataSetFromMatrix(countData = counts,
                             colData = coldata,
                             design = ~ condition)

# Calculating differential expression----------------------

#
dds <- DESeq(dds)
res <- results(dds)

res
summary(res)

#
#
resLFC <- lfcShrink(dds,
                    coef = "condition_post_flower_vs_pre_flower",
                    type = "apeglm")

# Vizualization--------------------------------------------
# MA-plot
# This plot shows log fold chang 
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

#PCA
#
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("condition"), returnDAta = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

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

ggsave(filename = '',
       device = 'png',
       width = 9,
       height = 7.5,
       dpi = 400,
       units = 'in')