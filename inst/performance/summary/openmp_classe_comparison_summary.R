library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)

breaks <- 2^(-30:20)
minor_breaks <- c()

# read the results
results <- read.table("results/openmp_classe_results.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
results$numStates <- factor(results$nstates)
results$numTips   <- factor(results$ntips)
results$numCores  <- factor(results$ncores)

# define the colors
# colors <- brewer.pal(6, "Set1")[c(2,5,3,4,6,1)]
# colors <- viridis(length(levels(results$numCores)), begin = 0.2, end = 0.9)
colors <- plasma(length(levels(results$numCores)), begin = 0.2, end = 0.9)

# plot runtimes
p <- ggplot(results, aes(x      = numStates,
                         y      = median_time,
                         color  = numCores,
                         group  = interaction(numTips, numCores))) +
  scale_y_continuous(breaks = breaks, minor_breaks = minor_breaks, trans = "log2") +
  geom_point(data = results %>% filter(numTips %in% c(32, 256, 1024)), aes(shape = numTips), stat = "summary", fun = "median") +
  geom_line(data = results %>% filter(numTips %in% c(32, 256, 1024)), stat = "summary", fun = "median") +
  scale_shape_manual(values = 1:6, guide = guide_legend(reverse = TRUE)) +
  ylab("average time per calculation (milliseconds)") +
  scale_color_manual(values = colors) +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  theme_minimal() +
  theme(legend.background = element_rect(fill = "white"),
        legend.box = "vertical",
        legend.justification = "center")

pdf("figures/openmp_classe_results.pdf", height = 5, width = 8)
plot(p)
dev.off()
