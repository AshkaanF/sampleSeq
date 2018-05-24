# Sampling model from 
## Daylight exposure modulates bacterial communities associated with household dust
### Fahimipour et al. *forthcoming*

This repository contains R scripts that recreate sampling models from our manuscript on the effects of light exposures on the dust microbiome.

#### Instructions
To run simply:

1. download this repository and set your working directory
2. `source(sampleSeq.R)`

The function `sampleSeq()` takes parameter values as described in the manuscript. Evaluating `sampleSeq(100)` will perform 100 model simulations with the same parameter ranges as in the amnuscript. An example usage and visualization is:

```R
## Perform 100 simulations qwith default parameters
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