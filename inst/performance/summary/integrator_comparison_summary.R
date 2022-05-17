library(dplyr)
library(ggplot2)
library(RColorBrewer)

breaks <- 2^(-30:20)
minor_breaks <- c()

# read the results
results <- read.table("results/integrator_results.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
results$numStates <- factor(results$nstates)
results$numTips   <- factor(results$ntips)
results$algo      <- factor(results$algo, levels = names(tensoRphylo::integrationScheme))

# define the colors
colors <- brewer.pal(5, "Set1")[c(5,2,3,1)]
names(colors) <- names(tensoRphylo::integrationScheme)

# plot runtimes
p <- ggplot(results, aes(x      = numStates,
                    y      = mean_time,
                    color  = algo,
                    label  = numStates,
                    symbol = numTips,
                    group  = interaction(algo, numTips))) +
  scale_y_continuous(breaks = breaks, minor_breaks = minor_breaks, trans = "log2") +
  geom_line(data = results %>% filter(numTips %in% c(32, 256, 1024)), stat = "summary", fun = "mean") +
  geom_point(data = results %>% filter(numTips %in% c(32, 256, 1024)), aes(shape = numTips), stat = "summary", fun = "mean") +
  scale_shape_manual(values = 1:6, guide = guide_legend(reverse = TRUE)) +
  scale_color_manual(values = colors) +
  ylab("average time per calculation (milliseconds)") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  guides(color = guide_legend(title="integrator")) +
  theme_minimal() +
  theme(legend.background = element_rect(fill = "white"),
        legend.box = "vertical",
        legend.justification = "center")

pdf("figures/integrator_results_time.pdf", height = 5, width = 8)
plot(p)
dev.off()


# plot errors
breaks <- 10^(-30:20)
p <- ggplot(results, aes(x      = numStates,
                    y      = error,
                    color  = algo,
                    label  = numStates,
                    symbol = numTips,
                    group  = interaction(algo, numTips))) +
  scale_y_continuous(breaks = breaks, minor_breaks = minor_breaks, trans = "log10") +
  geom_line(data = results %>% filter(numTips %in% c(32, 256, 1024)), stat = "summary", fun = "mean") +
  geom_point(data = results %>% filter(numTips %in% c(32, 256, 1024)), aes(shape = numTips), stat = "summary", fun = "mean") +
  scale_shape_manual(values = 1:6, guide = guide_legend(reverse = TRUE)) +
  scale_color_manual(values = colors) +
  ylab("error (absolute difference from true likelihood)") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  guides(color = guide_legend(title="integrator")) +
  theme_minimal() +
  theme(legend.background = element_rect(fill = "white"),
        legend.box = "vertical",
        legend.justification = "center")

pdf("figures/integrator_results_accuracy.pdf", height = 5, width = 8)
plot(p)
dev.off()
