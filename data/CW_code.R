# Importing and tidying the data-------------------------

# path variable where imported counts.txt would be found
# file.path(getwd() function is way to make code transferable to different operating systems
path <- file.path(getwd(), "data", "counts.tsv")

# Actual import
# sep, and escape sequence, to insert horizontal tab when tab key is pressed 
# header function is used if the first line of the file should be used for column names 
CW_counts <- read.table(file = path,
                        sep = '\t',
                        header = FALSE)

counts <- tail(counts, -4)

counts <- counts[ , c(1, seq(4,156,4))]

row_names <- counts[1]
row_names <- str_sub(row_names$V1, 6, -1) #tidyverse
rownames(counts) <- row_names
counts <- counts[ , 2:40]

colnames<- c("SRR3119152", "SRR3119153", "SRR3119154", "SRR3119155",
             "SRR3119156", "SRR3119157", "SRR3119158", "SRR3119159",
             "SRR3119160", "SRR3119161", "SRR3119162", "SRR3119163",
             "SRR3119164", "SRR3119165", "SRR3119166", "SRR3119167",
             "SRR3119168", "SRR3119169", "SRR3119170", "SRR3119171",
             "SRR3119172", "SRR3119173", "SRR3119174", "SRR3119175",
             "SRR3119176", "SRR3119177", "SRR3119178", "SRR3119179",
             "SRR3119180", "SRR3119181", "SRR3119182", "SRR3119183",
             "SRR3119184", "SRR3119185", "SRR3119186", "SRR3119187",
             "SRR3119188", "SRR3119189", "SRR3119190")
colnames(counts) <- col_names

metadata <- read.table(file = "./data/metadata.tsv",
                       sep = '\t'
                       header = TRUE)

subset_meta <- metadata[metadata$name == "lyc_S-" |
                          metadata$name == "lyc_S+",]
counts <- counts %>% select(c(subset_meta$run)) #tidyverse
colnames(counts) <- c("S_lyc_before_flower_rep_1",
                      "S_lyc_before_flower_rep_2",
                      "S_lyc_before_flower_rep_3",
                      "S_lyc_after_flower_rep_1",
                      "S_lyc_after_flower_rep_2",
                      "S_lyc_after_flower_rep_3")

coldata <- data.frame(row.names = colnames(counts),
                      conditon = factor(c(rep("pre_flower",3),
                                          rep("post_flower", 3))))
coldata$condition <- factor(coldata$conditon, levels = c("pre_flower",
                                                         "post_flower"))