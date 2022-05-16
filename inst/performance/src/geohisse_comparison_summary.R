library(dplyr)
library(ggplot2)
library(RColorBrewer)

breaks <- 2^(0:30)
# minor_breaks <- rep(1:9, 21)*(2^rep(-10:10, each=9))
minor_breaks <- c()

# read the results
results           <- read.table("geohisse_results.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
results$numTips   <- factor(results$ntips)
results$numHidden <- factor(results$nhidden)
results$method    <- factor(results$method, levels = c("tensorphylo", "hisse"))

# define the colors
colors <- brewer.pal(5, "Set1")[c(5,4)]
names(colors) <- c("tensorphylo", "hisse")

# plot runtimes
ggplot(results, aes(x      = numHidden,
                    y      = mean_time,
                    color  = method,
                    symbol = numTips,
                    group  = interaction(method, numTips))) +
  scale_y_continuous(breaks = breaks, minor_breaks = minor_breaks, trans = "log2") +
  geom_line(stat = "summary", fun = "mean") +
  geom_point(aes(shape = numTips), stat = "summary", fun = "mean") +
  scale_shape_manual(values = 1:6, guide = guide_legend(reverse = TRUE)) +
  scale_color_manual(values = colors) +
  ylab("average time per calculation (milliseconds)") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  guides(color = guide_legend(title="package")) +
  theme_minimal() +
  theme(legend.background = element_rect(fill = "white"),
        legend.box = "horizontal",
        legend.position = c(.8,.05),
        legend.justification = "bottom")

