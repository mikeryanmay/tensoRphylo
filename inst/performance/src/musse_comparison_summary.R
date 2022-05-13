library(dplyr)
library(ggplot2)

# read the results
results <- read.table("musse_results.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
results$nstates_f <- factor(results$nstates)
results$ntips_f   <- factor(results$ntips)

# plot runtimes
ggplot(results, aes(x      = ntips_f,
                    y      = mean_time,
                    color  = method,
                    label  = nstates_f,
                    symbol = nstates_f,
                    group  = interaction(method, nstates_f))) +
  scale_y_continuous(trans='log2') +
  geom_line(stat = "summary", fun = "mean") +
  geom_point(aes(shape = nstates_f), stat = "summary", fun = "mean")

# plot errors
ggplot(results, aes(x      = ntips_f,
                    y      = error,
                    color  = method,
                    label  = nstates_f,
                    symbol = nstates_f,
                    group  = interaction(method, nstates_f))) +
  scale_y_continuous(trans='log2') +
  geom_line(stat = "summary", fun = "mean") +
  geom_point(aes(shape = nstates_f), stat = "summary", fun = "mean")

# ggplot(results, aes(x      = ntips_f,
#                     y      = exp(likelihood - true),
#                     color  = method,
#                     label  = nstates_f,
#                     symbol = nstates_f,
#                     group  = interaction(method, nstates_f))) +
#   geom_line(stat = "summary", fun = "median") +
#   geom_point(aes(shape = nstates_f), stat = "summary", fun = "median")






