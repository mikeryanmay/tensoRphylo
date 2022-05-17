library(dplyr)
library(ggplot2)
library(RColorBrewer)

breaks <- 2^(-30:20)
# minor_breaks <- rep(1:9, 21)*(2^rep(-10:10, each=9))
minor_breaks <- c()

# read the results
results <- read.table("results/musse_results.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
results$numStates <- factor(results$nstates)
results$numTips   <- factor(results$ntips)
results$method    <- factor(results$method, levels = c("tensorphylo", "diversitree", "castor"))

# define the colors
colors <- brewer.pal(5, "Set1")[c(5,2,3)]
names(colors) <- c("tensorphylo", "diversitree", "castor")

# plot runtimes
p <- ggplot(results, aes(x      = numStates,
                    y      = mean_time,
                    color  = method,
                    label  = numStates,
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
        legend.box = "vertical",
        legend.justification = "center")


pdf("figures/musse_results.pdf", height = 5, width = 8)
plot(p)
dev.off()
