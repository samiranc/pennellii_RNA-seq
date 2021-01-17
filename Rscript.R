# Importing and tidying the data ----------------------
# Sets a file path that should work on all operating systems, [surely]

# path variable where imported counts.txt would be found
# file.path(getwd() function is way to make code transferable to different operating systems

path <- file.path(getwd(), "data", "counts.txt")

# Actual import 
?read.table()
# sep, and escape sequence, to insert horizontal tab when tab key is pressed 
# header function is used if the first line of the file should be used for column names 

counts <- read.table(file = path, sep = '\t', header = FALSE)

quit(save = yes)
