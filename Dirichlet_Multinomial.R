library(microbiome)
library(DirichletMultinomial)
library(tidyverse)
library(reshape2)

phyloShit

# To speed up, only consider the core taxa
# that are prevalent at 0.1% relative abundance in 50% of the samples
# (note that this is not strictly correct as information is
# being discarded; one alternative would be to aggregate rare taxa)
pseq.comp <- microbiome::transform(phyloShit, "compositional")
taxa <- core_members(pseq.comp, detection = 0.1/100, prevalence = 50/100)
pseq <- prune_taxa(taxa, phyloShit)

# Pick the OTU count matrix
# and convert it into samples x taxa format
dat <- abundances(phyloShit)
count <- as.matrix(t(dat))

#Fit the DMM model. Let us set the maximum allowed number 
#of community types to 3 to speed up the example.

fit <- mclapply(1:10, dmn, count = count, verbose=TRUE)

#Check model fit with different number of mixture components 
#using standard information criteria

lplc <- sapply(fit, laplace) # AIC / BIC / Laplace
aic  <- sapply(fit, AIC) # AIC / BIC / Laplace
bic  <- sapply(fit, BIC) # AIC / BIC / Laplace
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit", 
     col = "#57068c", lwd = 4, cex.main=2, cex.lab=1.7, cex.sub=1.2, main = "Determining Optimal Number of Dirichlet Components")
lines(aic, type="b", lty = 2, col = "#57068c", lwd = 4)
lines(bic, type="b", lty = 3, col = "#57068c", lwd = 4)

best <- fit[[which.min(lplc)]]
mixturewt(best)

ass <- apply(mixture(best), 1, which.max)

for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("OTU", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))     
  
  p <- ggplot(d, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity", fill="#57068c", colour="#57068c") +
    coord_flip() +
    theme(text = element_text(size=30)) +
    labs(title = paste("Community Type", k))
  print(p)
}


