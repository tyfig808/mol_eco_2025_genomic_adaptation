# does morpho share more snps within morpho and volatile
if (!require(tidyverse)) install.packages('tidyverse')
if (!require(vroom)) install.packages('vroom')
if (!require(data.table)) install.packages('data.table')
if (!require(igraph)) install.packages('igraph')
if (!require(bipartite)) install.packages('bipartite')
if (!require(emmeans)) install.packages('emmeans')
if (!require(ggsignif)) install.packages('ggsignif')
if (!require(car)) install.packages('car')
if (!require(effects)) install.packages('effects')
if (!require(ggridges)) install.packages('ggridges')
if (!require(glmmTMB)) install.packages('glmmTMB')
if (!require(DHARMa)) install.packages('DHARMa')
if (!require(UpSetR)) install.packages('UpSetR')

# set wd to place outfile there
file_wd <- "~/mol_eco_2024/2.5.5_gwas/"
out_wd <- "~/mol_eco_2024/2.5.5_gwas/"

# load trait matrix full file
setwd(out_wd)
t_mat <- fread(sep = "\t", paste0("shared_snps_between_traits_sigzone_only.tsv"))

# add merged trait column
t_mat[, pair_trait := paste(trait_1, trait_2, sep = "_")]

# plot to explore, seems like there is a difference
ggplot(t_mat, aes(x = snps_shared, y = shared_class)) + 
  geom_density_ridges_gradient() +
  #scale_fill_viridis_c(name = "F Value", option = "C") +
  labs(x = "# of SNPs shared between trait classes", y = NULL) + 
  theme_minimal(base_size = 19)

# set contrast so that it is within vs between 
t_mat$shared_class <- as.factor(t_mat$shared_class)

# check, make so morph_vol vs the others, then check again
contrasts(t_mat$shared_class)
t_mat$shared_class <- relevel(t_mat$shared_class, ref = "morphology_volatile")
contrasts(t_mat$shared_class)

# check to see what family dist is best
m1 <- glmmTMB(snps_shared ~ factor(shared_class) + (1|treat) + (1|pair_trait/treat), data = t_mat,
		family = "nbinom2") # tweedie seems to have convergence problems so lets go with nbinom2, also looks okay

Anova(m1)
#                      Chisq Df Pr(>Chisq)    
#factor(shared_class) 154.31  2  < 2.2e-16 ***

sim <- simulateResiduals(fittedModel = m1, plot = T, n = 100)

summary(m1)
#                                         Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                                -1.4178     0.1339  -10.59   <2e-16 ***
#factor(shared_class)morphology_morphology   1.3325     0.1506    8.85   <2e-16 ***
#factor(shared_class)volatile_volatile       1.7756     0.1561   11.37   <2e-16 ***

# emmeans 
em <- (emmeans(m1,  ~ shared_class, type = "response"))
em_class <- data.table(as.data.frame(emmeans(m1,  ~ shared_class, type = "response")))
em_class$group <- as.numeric(data.table(as.data.frame(multcomp::cld(em)))$.group) 
em_class$group <-letters[em_class$group]; em_class

#            shared_class  response         SE  df asymp.LCL asymp.UCL group
#1:   morphology_volatile 0.2422523 0.03244679 Inf 0.1863201  0.314975     a
#2: morphology_morphology 0.9182586 0.13906706 Inf 0.6824228  1.235596     b
#3:     volatile_volatile 1.4302074 0.21698006 Inf 1.0623344  1.925470     c

# see if violin plot is better than geom ridges, this removes the 0 though
ggplot(t_mat, aes(x = shared_class, y = snps_shared)) +
  geom_violin(draw_quantiles = c(0.5))

# this includes all the zeros
resolution <- 144
setwd(out_wd)

svglite(paste0("shared_regions_trait_classes_log.svg", sep=""), width = 1440/resolution, height = 960/resolution)
ggplot(t_mat, aes(x = shared_class, y = log(snps_shared+1))) +
  geom_violin(draw_quantiles = c(0.5)) +
  labs(x = "Pairs of Trait Classes", 
    y = "Log(# of Shared Significant Regions)") + 
  theme_bw(base_size = 19)

dev.off()

t_mat[, .N, by = c("shared_class")]
#            shared_class    N
#1: morphology_morphology  728
#2:   morphology_volatile 1456
#3:     volatile_volatile  624

t_mat[, mean(snps_shared), by = "shared_class"]
#            shared_class       V1
#1: morphology_morphology 1.971154
#2:   morphology_volatile 0.345467
#3:     volatile_volatile 2.863782
