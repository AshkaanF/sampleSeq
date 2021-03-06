# Sampling model from 
#### Daylight exposure modulates bacterial communities associated with household dust
Fahimipour et al. (2018) doi: [10.1186/s40168-018-0559-4](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0559-4)

This repository contains R scripts that recreate sampling models from our manuscript on the experimental effects of light exposures on the dust microbiome.

#### Instructions
To run:

1. download this repository and set your working directory
2. `source('./R/util_func.R')`

The function `sampleSeq()` takes parameter values as described in the manuscript. Evaluating `sampleSeq(100)` will perform 100 model simulations with the same parameter ranges as in the manuscript. An example usage and visualization is:


```R
## Perform 100 simulations with default parameters
res <- sampleSeq(100)

## color palette
col.dude <- magma(15) %>% rev()

##-
## ggplot
##-
ggplot(data = res, aes(x = log10(1 + nw), y = lf)) +
  theme_classic() +
  xlab(expression('log'['10']*' 1 + RSV Abundance')) +
  ylab(expression('Apparent log'['10']*'-fold change following loss of abundant RSVs')) +
  geom_hex(aes(fill = (..count..)/sum(..count..)), bins = 20, alpha = 0.85, size = 0.5, colour = '#f7f7f7') +
  geom_hline(yintercept = 0, linetype = 3, size = 0.75) +
  scale_fill_gradientn(colours = col.dude, trans = 'identity', name = 'Frequency', breaks = pretty_breaks(n = 5)) +
  theme(axis.title = element_text(size = 10))
  ```