if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')
if (!require(igraph)) install.packages('igraph')
if (!require(bipartite)) install.packages('bipartite')
if (!require(scales)) install.packages('scales')
if (!require(svglite)) install.packages('svglite')
if (!require(geomorph)) install.packages('geomorph')
if (!require(emmeans)) install.packages('emmeans')
if (!require(car)) install.packages('car')
if (!require(glmmTMB)) install.packages('glmmTMB')
if (!require(DHARMa)) install.packages('DHARMa')
if (!require(fitdistrplus)) install.packages('fitdistrplus')

out_wd <- "~/mol_eco_2024/2.5.5_gwas/"
setwd(out_wd)

#
pleio_index <- fread(sep = ",", paste0("pleio_index_sigzone_only.csv")) 

#
size_dt <- fread(sep = "\t", paste0("network_size_treatments_sigzone_only.tsv"))
size_dt$pop_2 <- paste(size_dt$soil, size_dt$herb, size_dt$bee, sep = "")


size_ind <- glmmTMB(1/size_dt$size ~ pleio_index$ind, family = beta_family())
Anova(size_ind)
#pleio_index$ind 13.219  1  0.0002771 ***

