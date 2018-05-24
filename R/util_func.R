##-
## initialize
##-
library(dplyr)
library(ggplot2)
library(scales)
library(viridis)
set.seed(777)

##-
## logistic viability function
##-
logist <- function(x, L = 0.5, k = -0.01, b = 1e2, floor = 0.01){
  y <- floor + 
    (
      (L - floor) /
        (1 + 
           (exp(1) ^ (-k * (x - b)))))
  
  y
}

##-
## extract top k values from vector
##-
which.part <- function(x, n = 30) {
  nx <- length(x)
  p <- nx-n
  xp <- sort(x, partial=p)[p]
  which(x > xp)
}

##-
## whole number function
##-
is.wholenumber <- function(x, tol = .Machine$double.eps ^ 0.5)  abs(x - round(x)) < tol

##-
## main loop function
##-
sampleSeq <- function(n.sims = 10, ## number of simulations
                      l.mean = c(3, 6),
                      l.sd = c(1, 2),
                      k.in = c(-5e-3, -5e-4),
                      L.in = c(0.075, 0.75),
                      floor.in = c(0.001, 0.05),
                      b.in = c(1e2, 2.5e3),
                      S.in = c(500, 1000),
                      mortality = c(15, 65),
                      seq.depth = 5e4){ 
  
  ## list for abundance ratios
  ab.rat <- list()
  
  ##-
  ## initiate loop n.sim times
  ##-
  for(u in 1:n.sims){
    
    ## suppress errors
    tryCatch({
      
      ## for S species
      S <- sample(S.in[1]:S.in[2], size = 1)
      
      ## draw SAD log mean
      log.mean <- runif(1, min = l.mean[1], max = l.mean[2])
  
      ## draw SAD log sd
      log.sd <- runif(1, min = l.sd[1], max = l.sd[2])
      
      ## draw k parm
      k.parm <- runif(1, min = k.in[1], max = k.in[2])
      
      ## draw L parameter
      L.parm <- runif(1, min = L.in[1], max = L.in[2])
      
      ## draw floor parameter
      floor.parm <- runif(1, min = floor.in[1], max = floor.in[2])
      
      ## draw b parameter
      b.parm <- runif(1, min = b.in[1], max = b.in[2])
      
      ##-
      ## start by establishing the 'ground truth' distributions
      ## of living (viab) and total (SAD) 'microbes'
      ##-
      ## generate a SAD from log normal dist and scale by A
      SAD <- rlnorm(S, log.mean, sdlog = log.sd) %>% round(0)
      
      ## set relationship between abundance class and % viable
      viab <- (SAD * (SAD %>% logist(L = L.parm,
                                     floor = floor.parm,
                                     k = k.parm,
                                     b = b.parm))
               ) %>% round(0)
      
      ##-
      ## generate vectors of individuals
      ##-
      ## generate 'species' vectors to sample SAD
      SAD.indv <- rep(1:S, SAD)
      
      ## generate 'species' vectors to sample viable SAD
      viab.indv <- rep(1:S, viab)
      
      ##---
      ## pool samples at equal depth.
      ## this simulates the lab technician
      ## drawing sequences from DNA aliquots
      ## pooled into equal concentrations
      ##---
      sequencer.depth <- seq.depth
      
      ## draw from SAD (total community)
      draw.total <- sample(SAD.indv, size = sequencer.depth, replace = F) 
      draw.total.tab <- draw.total %>% table()
      
      ## draw from viab (living community)
      draw.living <- sample(viab.indv, size = sequencer.depth, replace = F)
      draw.living.tab <- draw.living %>% table()
      
      ##-
      ## kill bugs here
      ##-
      ## number of tax to remove
      kill.N <- sample(mortality[1]:mortality[2], size = 1)
      
      ## drop k most abundant taxa
      SAD.kill <- SAD[-c(which.part(SAD, n = kill.N))]
      
      ## drop from viab
      viab.kill <- viab[-c(which.part(SAD, n = kill.N))]
      
      ##-
      ## generate vectors of individuals
      ##-
      ## generate 'species' vectors to sample SAD
      SAD.indv.kill <- rep(1:(S - kill.N), SAD.kill)
      
      ## generate 'species' vectors to sample viable SAD
      viab.indv.kill <- rep(1:(S - kill.N), viab.kill)
      
      ## draw from SAD (total community)
      draw.total.2 <- sample(SAD.indv.kill, size = sequencer.depth, replace = F) 
      draw.total.tab.2 <- draw.total.2 %>% table
      
      ## draw from viab (living community)
      draw.living.2 <- sample(viab.indv.kill, size = sequencer.depth, replace = F)
      draw.living.tab.2 <- draw.living.2 %>% table
      
      ##---
      ## calculate abundance ratio
      ##---
      lf.change <- log10(draw.total.tab.2[intersect(names(draw.total.tab), names(draw.total.tab.2))] /
                           draw.total.tab[intersect(names(draw.total.tab), names(draw.total.tab.2))]
                         ) %>% c()
      
      ## extract abundance
      n.w.loss <- draw.total.tab[intersect(names(draw.total.tab), names(draw.total.tab.2))] %>% c()
      
      ## save
      ab.rat[[u]] <- data.frame('lf' = lf.change, 'nw' = n.w.loss, 'sim' = rep(as.character(u)))
      
      ## close error suppression
    }, error = function(e){})
    
    ## drop us a line
    if(is.wholenumber(u / 10)){
      cat(u, 'simulation(s) completed \n')
    }
    
  }
  
  ## output sim results
  ab.rat %>% do.call('rbind', .) %>% as.data.frame()

}




