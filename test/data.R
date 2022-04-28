fn <- "test/data.nex"
data <- readNexusData(fn)

fn <- "test/data.csv"
data_csv <- readDelimitedData(fn, ",", 3)
head(data_csv)

fn <- "test/data.tsv"
data_tsv <- readDelimitedData(fn, "\t", 3)

head(data)
head(data_csv)
head(data_tsv)
