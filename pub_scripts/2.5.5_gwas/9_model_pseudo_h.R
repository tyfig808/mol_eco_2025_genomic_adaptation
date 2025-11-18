if (!require(data.table)) install.packages('data.table')
if (!require(car)) install.packages('car')
if (!require(glmmTMB)) install.packages('glmmTMB')
if (!require(DHARMa)) install.packages('DHARMa')
if (!require(fitdistrplus)) install.packages('fitdistrplus')
if (!require(ggplot2)) install.packages('ggplot2')
if (!require(viridis)) install.packages('viridis')
if (!require(ggpubr)) install.packages('ggpubr')
if (!require(scales)) install.packages('scales')
if (!require(ggnewscale)) install.packages('ggnewscale')
if (!require(svglite)) install.packages('svglite')
if (!require(extrafont)) install.packages('extrafont')
if (!require(stringr)) install.packages('stringr')

# model  hertiabiltiy for the trait treatments
setwd("~/mol_eco_2024/2.5.5_gwas")
dt <- fread(sep = "\t", paste0("emmax_estimates_all_treat_trait.tsv"))

#add factor cols
dt$soil <- str_sub(dt$treat, 1, 1)
dt$herb <- str_sub(dt$treat, 2, -2)
dt$bee <- str_sub(dt$treat, -1)

# model
m_beta <- glmmTMB(pseud_h ~ soil*herb*bee - soil:herb:bee + (1|trait), data = dt, family = beta_family())
Anova(m_beta)

#                Chisq Df Pr(>Chisq)    
#soil           3.2760  1  0.0703018 .  
#herb          16.0275  1  6.243e-05 ***
#bee           14.3415  1  0.0001525 ***
#soil:herb      4.2675  1  0.0388491 *  
#soil:bee       0.1753  1  0.6754229    
#herb:bee       2.4631  1  0.1165505    
#soil:herb:bee  0.0525  1  0.8187020    

simulateResiduals(fittedModel = m_beta, plot = T, n = 100)

# check averages, herb makes lower 
dt[, mean(pseud_h), by = c("soil", "herb", "bee")]

# check averages, herb makes lower 
mean_dt <- dt[, mean(pseud_h), by = c("soil", "herb", "bee")]
