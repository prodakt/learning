# ///////////////////////////////////////
# ANOVA 
# ///////////////////////////////////////


# http://www.sthda.com/english/wiki/one-way-anova-test-in-r
install.packages("tidyselect")

library(dplyr)
my_data <- PlantGrowth
head(my_data)
set.seed(1234)
dplyr::sample_n(my_data, 10)
# Show the levels
levels(my_data$group)
my_data$group <- ordered(my_data$group, levels = c("ctrl", "trt1", "trt2")) # jezeli nei sa dobrze posortowane to posortowac

group_by(my_data, group) %>%
  summarise(
    count = n(),
    mean = mean(weight, na.rm = TRUE),
    sd = sd(weight, na.rm = TRUE)
  )

# Compute the analysis of variance
res.aov <- aov(weight ~ group, data = my_data)
# Summary of the analysis
# saov <- as.data.frame(t(as.data.frame(unlist(summary(res.aov)))))
# class(saov)
# saov$`Pr(>F)1`


saov <- summary(res.aov)
saov[[1]]$'Pr(>F)'[1]



# https://www.rdocumentation.org/packages/HybridMTest/versions/1.16.0/topics/row.oneway.anova
####################Three group comparison###################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("HybridMTest", version = "3.8")


install.packages("HybridMTest")
library(HybridMTest)

# load data
data(GroupComp.data)
# Read the expression values   
brain.express.set <- exprs(GroupComp.data)
head(brain.express.set)
# Read the phenotype
brain.pheno.data <- pData(GroupComp.data)
brain.pheno.data[,1] 
head(brain.pheno.data)
brain.pheno.data

# ANOVA test
row.oneway.anova(brain.express.set,brain.pheno.data[,1])



