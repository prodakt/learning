# t-SNE

# https://www.analyticsvidhya.com/blog/2017/01/t-sne-implementation-r-python/
# https://drive.google.com/file/d/0B6E7D59TV2zWYlJLZHdGeUYydlk/view?usp=sharing
install.packages('matrixStats')
library("Rtsne")
library("matrixStats")

## calling the installed package
train <- read.csv(file.choose()) ## Choose the train.csv file downloaded from the link above  
train <- read.csv(file.choose(), sep = "\t")
colnames(train)[1] <- "label"
head(train)
summary(train)
nrow(train)
ncol(train)

## Curating the database for analysis with both t-SNE and PCA
Labels <- train$label
train$label <- as.factor(train$label)

head(Labels)

## for plotting
colors = rainbow(length(unique(train$label)))
names(colors) = unique(train$label)

## Executing the algorithm on curated data
head(train)
tsne <- Rtsne(train[,-1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
exeTimeTsne <- system.time(Rtsne(train[,-1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500))

## Plotting
plot(tsne$Y, t='n', main="tsne")
text(tsne$Y, labels=train$label, col=colors[train$label])

tsne_tab <- as.data.frame(tsne$Y)
head(tsne_tab)

head(tsne)

## plotting the results without clustering
library(ggplot2)
ggplot(tsne_tab, aes(x=V1, y=V2)) +
  geom_point(size=0.25) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("t-SNE") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_brewer(palette = "Set2")


# -----------------
# # load your omic data here as mydata
install.packages('M3C')
library(M3C)
tsne(pollen$data,colvec=c('gold'))
res <- M3C(mydata)
head(res)

tsne(pollen$data,colvec=c('gold'))


# --------------
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)

seqdata <- read.csv2(file = "Downloads/GSE60450_Lactation-GenewiseCounts.txt", sep = "\t")
countdata <- seqdata[,-(1:2)]
rownames(countdata) <- seqdata[,1]
head(countdata)
colnames(countdata) <- substr(colnames(countdata),start=1,stop=7)
myCPM <- cpm(countdata)

fpkm <- read.table("Downloads/genes.fpkm_table", row.names = 1, stringsAsFactors = F)
kolnames <- as.character(unlist(fpkm[1,]))
kolnames
fpkm <- fpkm[-1,]
colnames(fpkm) <- kolnames
head(fpkm)

tab <- apply(fpkm, 2, as.numeric)
rownames(tab) <- rownames(fpkm)
class(tab[1,1])
head(tab)

fpkm <- tab[rowSums(tab) > 1,]
head(fpkm)
fpkm <- t(fpkm)


tsne <- Rtsne(fpkm, dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)

## Plotting
plot(tsne$Y, t='n', main="tsne")
text(tsne$Y, labels=train$label, col=colors[train$label])

tsne_tab <- as.data.frame(tsne$Y)
head(tsne_tab)

head(tsne)

## plotting the results without clustering
library(ggplot2)
ggplot(tsne_tab, aes(x=V1, y=V2)) +
  geom_point(size=0.25) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("t-SNE") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_brewer(palette = "Set2")





install.packages("githubinstall")
library(githubinstall)
githubinstall("M3C")





# -------------------
# https://ajitjohnson.com/tsne-for-biologist-tutorial/


read.csv("https://github.com/ajitjohnson/ajitjohnson.github.io/blob/master/assets/data/tsne_tutorial/exp.csv")






# Rtsne -------
iris_unique <- unique(iris) # Remove duplicates
set.seed(42) # Sets seed for reproducibility
tsne_out <- Rtsne(as.matrix(iris_unique[,1:4])) # Run TSNE
plot(tsne_out$Y,col=iris_unique$Species,asp=1) # Plot the result
class(as.matrix(iris_unique[,1:4])[1,1])

exprTab <- read.csv(file.choose(), sep = "\t", row.names = 1)
extab <- exprTab[rowSums(exprTab) > 1,]
extab2 <- cbind(extab, extab)
extab2 <- unique(extab2)

expr_tsne <- Rtsne(as.matrix(extab2), perplexity = 10)
plot(expr_tsne$Y, asp=1) # Plot the result


head(extab2)

extab <- unique(extab)

exprTab <- t(extab[1:150,])
head(exprTab)
exprTab_unique <- unique(exprTab)
head(exprTab_unique)

expr_tsne <- Rtsne(as.matrix(exprTab_unique), perplexity = 3)
plot(expr_tsne$Y, asp=1) # Plot the result

plot(expr_tsne$Y, col = colnames(exprTab), asp=1) # Plot the result

plot(tsne_out$Y,col=iris_unique$Species,asp=1) # Plot the result


head(extab)
extab <- as.matrix(extab)
sortowanie <- data.frame(Vara = rowVars(extab[,1:6]), Varb = rowVars(extab[,7:12]), 
                         sra = rowMeans2(extab[,1:6]), srb = rowMeans2(extab[,7:12]))
sortowanie$aplusb <- rowSums2(as.matrix(sortowanie[,1:2]))
sortowanie$sraMINsrb <- abs(rowDiffs(as.matrix(sortowanie[,3:4])))
sortowanie$index <- sortowanie[,6]/sortowanie[,5]
head(sortowanie)

rv <- sortowanie[,1:2]
select = order(rv, decreasing=TRUE)[1:500]
head(select)
select
length(select)

seq_len(500)

#select = order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
tSNE <- Rtsne(as.matrix(t(extab[select,])), perplexity = 3, dims = 3) # tSNE analysis

print_tSNE <- function(x, components= c(1,2,3), dims = 3, ntop = 500){
  set.seed(8) # set a random number
  # the function prepares the data, computes the selected principal components
  rv = rowVars(x) # get variance in each row
  select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  tSNE <- Rtsne(as.matrix(t(x[select,])), perplexity = 3, dims = dims) # tSNE analysis
  
  #get the set of the three selected components
  tSNEx <- tSNE$Y[,components[1]]
  tSNEy <- tSNE$Y[,components[2]]
  tSNEz <- tSNE$Y[,components[3]]
  
  tSNE_out <- data.frame(row.names = colnames(x), tSNEx = tSNEx, tSNEy = tSNEy, tSNEz = tSNEz)
  return(tSNE_out)
}

plotData <- print_tSNE(extab)
head(plotData)

plotData$colors <- c(rep("red",6), rep("blue",6))
plotData
plot(plotData$tSNEx, plotData$tSNEy, col = plotData$colors)


rv = rowVars(extab) # get variance in each row
select = order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
head(select)
head(rv)
