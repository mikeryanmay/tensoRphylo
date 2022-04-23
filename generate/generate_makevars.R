# simple script to generate Makevars with all the appropriate sources
# to modify the Makevars, modify the one in generate/Makevars, then run this script

library(tools)

# step 1: find all the .cpp files in src
cpp <- list.files("src/", recursive = TRUE)
cpp <- cpp[file_ext(cpp) == "cpp"]
cpp <- cpp[grepl("RcppExports", cpp) == FALSE]

sources <- paste0(cpp, collapse = " ")

# step 2: read the Makevars file
Makevars <- readLines("generate/Makevars")

# step 3: replace the source placeholder
new_Makevars <- gsub("SOURCEPLACEHOLDER", sources, Makevars)

# step 4: write the Makevars
writeLines(new_Makevars, con = "src/Makevars")
