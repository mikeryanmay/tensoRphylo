library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)

breaks <- 2^(-30:20)
# minor_breaks <- rep(1:9, 21)*(2^rep(-10:10, each=9))
minor_breaks <- c()

# read the results
results <- read.table("results/constant_rate_results.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# define the colors
colors <- brewer.pal(5, "Set1")[c(5,2)]
names(colors) <- c("TESS", "TreePar")

# plot runtimes
p <- ggplot(results, aes(x = x_likelihood, y = y_likelihood, color = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2) +
  geom_point(aes(shape = method)) +
  scale_shape_manual(values = 4:3) +
  scale_color_manual(values = colors) +
  xlab("log probability density (tensoRphylo)") + ylab("log probability density (analytical)") +
  theme_minimal() +
  theme(legend.background = element_rect(fill = "white"),
        legend.box = "vertical",
        legend.position = c(0.05, 0.9),
        legend.justification = "left")

pdf("figures/constant_rate_results.pdf", height = 6, width = 6)
plot(p)
dev.off()
