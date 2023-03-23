library(tidyverse)
library(jtools)
library(haven)
df <- read_sav("Data/SPSS_Chapter2.sav")
df.ch2 <- df %>%
        select(PersonID , cognition , age , grip , sexMW , demgroup , age80 , age85 , age90 , grip6 , grip9 , grip12 , sexWM , demNF , demNC , demFN , demFC , demCN , demCF)

summ(mod <- lm(cognition ~ age*grip9 + sexMW + demNF + demNC, df.ch2))
sim_slopes(mod, pred = grip9, modx = age, jnplot = T)
interact_plot(mod, pred = grip9, modx = age)
