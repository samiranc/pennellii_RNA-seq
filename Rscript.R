# Importing and tidying the data ----------------------
# Sets a file path that should work on all operating systems, [surely]

# path variable where imported counts.txt would be found
# file.path(getwd() function is way to make code transferable to different operating systems

path <- file.path(getwd(C:/Users/RAFAEL/Documents/pennellii_RNA-seq/R/data/counts.txt), "data", "counts.txt")

# Actual import 
?read.table()
# sep, and escape sequence, to insert horizontal tab when tab key is pressed 
# header function is used if the first line of the file should be used for column names 

counts <- read.table(C:\Users\RAFAEL\Documents\R\pennellii_RNA-seq\data\counts.txt = path, sep = '\t', header = FALSE)