# Importing and tidying the data-------------------------

# path variable where imported counts.txt would be found
# file.path(getwd() function is way to make code transferable to different operating systems
path <- file.path(getwd(), "data", "counts.tsv")

# Actual import
# sep, and escape sequence, to insert horizontal tab when tab key is pressed 
# header function is used if the first line of the file should be used for column names 
counts <- read.table(file = path,
                        sep = '\t',
                        header = FALSE)

# Subtracted the first four lines/rows of counts
counts <- tail(counts, -4)

# From 156 variable to 40? The seq() function should've created a sequence of 4-156. Don't understand what the comma is for
# Created a sequence of 4-156 with intervals of 4
# The c(1,) kept the first column
# Without [ , ] error occurs of undefined columns selected
counts <- counts[ , c(1, seq(4, 156, 4))]

# Assigning row_names to be the first element in counts
row_names <- counts[1]
# str_sub: Extract and replace substrings from a character vector (This function starts with the name of V1 at the sixth letter and end in the last one)
# The $ helps extract v1
row_names <- str_sub(row_names$V1, 6, -1) #tidyverse
# Making the row names of counts to be equivalent to row_names above
# This is were it actually changes the counts file, the column before v1 has the names of the genes without "gene:" in front of it
rownames(counts) <- row_names
# Included columns 2 through 40
# Took out column v1
counts <- counts[ , 2:40]

# Assigning col_names the values of the SRR numbers
col_names<- c("SRR3119152", "SRR3119153", "SRR3119154", "SRR3119155",
             "SRR3119156", "SRR3119157", "SRR3119158", "SRR3119159",
             "SRR3119160", "SRR3119161", "SRR3119162", "SRR3119163",
             "SRR3119164", "SRR3119165", "SRR3119166", "SRR3119167",
             "SRR3119168", "SRR3119169", "SRR3119170", "SRR3119171",
             "SRR3119172", "SRR3119173", "SRR3119174", "SRR3119175",
             "SRR3119176", "SRR3119177", "SRR3119178", "SRR3119179",
             "SRR3119180", "SRR3119181", "SRR3119182", "SRR3119183",
             "SRR3119184", "SRR3119185", "SRR3119186", "SRR3119187",
             "SRR3119188", "SRR3119189", "SRR3119190")
# Now col_names will be the column names in counts
colnames(counts) <- col_names

# This is to import the metadata.tsv file
metadata <- read.table(file = "./data/metadata.tsv",
                       sep = '\t',
                       header = TRUE)

# Assigns subset_meta with unpollinated pistil rows from metadata set
# Does it create a metadata, or what is one?
subset_meta <- metadata[metadata$name == "lyc_S-" |
                          metadata$name == "lyc_S+",]
# Now there are only 6 columns of two samples, 3 reps, of lycopersicum
counts <- counts %>% select(c(subset_meta$run)) #tidyverse
# colnames() assigns the column names as the values below
colnames(counts) <- c("S_lyc_before_flower_rep_1",
                      "S_lyc_before_flower_rep_2",
                      "S_lyc_before_flower_rep_3",
                      "S_lyc_after_flower_rep_1",
                      "S_lyc_after_flower_rep_2",
                      "S_lyc_after_flower_rep_3")

# coldata is a data frame with row names as colnames above and one column as the condition of pre and post flowering
# rep() = repetitions?
coldata <- data.frame(row.names = colnames(counts),
                      condition = factor(c(rep("pre_flower", 3),
                                          rep("post_flower", 3))))
# factor()? 
coldata$condition <- factor(coldata$condition, levels = c("pre_flower",
                                                         "post_flower"))
