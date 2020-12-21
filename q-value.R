if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("qvalue")
install.packages("readxl")
install.packages("xlsx")


library("qvalue")
library("readxl")
library("xlsx")

my_data <- read_excel("Downloads/Matrix US45103052_252644010691_normalized-data.xlsx")
head(my_data, 20)
tab <- my_data[13:nrow(my_data),]
head(tab)
tail(tab)

tab <- as.data.frame(tab)
row.names(tab) <- tab[1,]

write.csv2(x = tab, file = "tab.csv", quote = F, row.names = F)
head(read.csv2("tab.csv", header = F))

tab2 <- read.csv2("Downloads/Matrix US45103052_252644010691_normalized-data/Arkusz1-Tabela 1.csv", header = T)
head(tab2)
class(tab2$p)
tab2$p <- as.numeric(as.character(tab2$p))

head(lfdr(tab2$p))
head(qvalue(tab2$p))

qval <- qvalue(tab2$p)

tab2$q_value <- qval$qvalues
tab2$FDR <- qval$lfdr

head(tab2)

write.csv2(tab2, "Downloads/wynik_z_FDR_qvalue.csv", quote = F, row.names = F)
write.table(tab2, "Downloads/wynik_z_FDR_qvalue.tsv", quote = F, sep = "\t", row.names = F)


hist(tab2$q_value)
length(qval$qvalues < 0.05)

nrow(tab2[tab2$q_value < 0.05,])
nrow(tab2[tab2$FDR < 0.25,])
