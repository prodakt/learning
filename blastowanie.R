##
## info
##

#wczytywanie bibliotek
source("http://bioconductor.org/biocLite.R")

biocLite("stringr")
biocLite("annotate")
biocLite("Biostrings")
biocLite("seqinr")

library("stringr")
library("annotate")
library("Biostrings")
library("seqinr")

getwd()
setwd("G:/NGS/Ciereszko/dane/download_2015-07-28_17-04-00/")

##############################################
## wczytanie danych
NazwaInputu <- "NoGenes_12.csv"
NazwaOutputu <- paste("BLAST_res_",NazwaInputu, sep = "", collapse = "")

contigs <- read.csv(NazwaInputu, sep = ";")                               # lista contigów bez danoacji, ktore wykazuja roznice ekspresji
Sscrofa10.2 <- read.csv("G:/NGS/Ciereszko/dane/Sscrofa10_2.csv", sep=";")    # lista chromosomow z genomu (konkretna wersja)

###########################################################

if (file.exists(NazwaOutputu)){
  warnings("Taki plik istnieje, zostanie nadpisany")
} else {file.create(NazwaOutputu, showWarnings = TRUE)}


for (i in 1:length(contigs$locus)){
  lokus <- contigs[i,]$locus
  # wyciagniecie wartosci: numer chromosomu lub ID, pozycji 'od' i 'do' na chromosomie
  chr<-substr(lokus, 1, str_locate(lokus, ":")[1,]-1)
  od<-substr(lokus,str_locate(lokus, ":")[1,]+1,str_locate(lokus,"-")[1,]-1)
  do<-substr(lokus,str_locate(lokus,"-")[1,]+1,str_length(lokus))
  
  if(str_length(chr)<4){numerID<- Sscrofa10.2[Sscrofa10.2$Name==chr,]$RefSeq} else{numerID<-chr}
  
  adres <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=",numerID,"&strand=2&seq_start=",od,"&seq_stop=",do,"&rettype=fasta&retmode=text", sep="", collapse = "")
  #seq <- readDNAStringSet(filepath = adres, format = "fasta")
  seq <- read.fasta(file = adres, seqtype = "DNA", as.string = TRUE, seqonly = TRUE)
  
  if (str_length(seq)>200){
    dlugiFASTA <- read.fasta(file = adres, seqtype = "DNA", as.string = TRUE)
    write.fasta(sequences = dlugiFASTA, file.out = paste("Long_seq",chr,"=",od,"-",do,".fasta", sep = "", collapse = ""), names = names(dlugiFASTA)) 
  }else{
      blastRes <- blastSequences(x = seq, database = "nr", hitListSize = 10, as = "data.frame", timeout = 300)
      blastRes[,length(colnames(blastRes))+1] <- lokus
      colnames(blastRes)[length(blastRes)] <- "locus"
  
      # write.csv(blastRes, paste("BLAST_res_",chr,"=",od,"-",do,".csv", sep = "", collapse = ""))
      tmp_output <- read.csv(NazwaOutputu)
      tmp_output <- rbind(tmp_output, blastRes)
      write.csv(tmp_output, NazwaOutputu)
  }
  

}

###################################
### czarnym blastem BLAST+
###################################

#test na 1 seq
setwd("C:/Program Files/NCBI/blast-2.4.0+/bin/")
zapyt <- paste("blastn",
               "-db nt", 
               "-query G:/NGS/Ciereszko/dane/Swinie2/1seq.fasta", 
               "-out G:/NGS/Ciereszko/dane/Swinie2/res.out",
               "-task blastn", #predefiniowane ustawienia
               #"-outfmt 7",
               '-outfmt "6 stitle sscinames scomnames staxids qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"',
               "-remote -max_target_seqs 2")
system(zapyt)
test <- read.table("G:/NGS/Ciereszko/dane/Swinie2/res.out", sep = "\t")
test

#test na wiele seq
zapyt <- paste("blastn -db nt", 
               "-query G:/NGS/Ciereszko/dane/Swinie2/312_unknow_seq.fasta", 
               "-out G:/NGS/Ciereszko/dane/Swinie2/wyblastowane.out",
               "-task blastn",
               #"-outfmt 7",
               '-outfmt "6 stitle sscinames scomnames staxids qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"',
               "-remote -max_target_seqs 1")
paste(zapyt)
system(zapyt)
blastn624nc <- read.table("G:/NGS/Ciereszko/dane/Swinie2/res.out", sep = "\t")
colnames(blastn624nc) <- c("title", "orgn", "orgn com", "taxid", "q seq id", "s seq id", "p ident", "length", "mismatch", "gap open","qstart","qend","sstart","send","E-val","score")
head(blastn624nc)


# blast x
zapytx <- paste("blastx -db nr", #refseq_protein", 
                "-query G:/NGS/Ciereszko/dane/Swinie2/312_unknow_seq1.fasta", 
                "-out G:/NGS/Ciereszko/dane/Swinie2/wyblastowane_312_unknown.out",
                '-outfmt "6 stitle sscinames scomnames staxids qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"',
                "-remote",
                "-max_target_seqs 10"#,
                #"-dbsize 1",
                #"-max_hsps 1"
                )
paste(zapytx)
system(zapytx)
blastn312nc <- read.table("G:/NGS/Ciereszko/dane/Swinie2/res.out", sep = "\t")
colnames(blastn312nc) <- c("title", "orgn", "orgn com", "taxid", "q seq id", "s seq id", "p ident", "length", "mismatch", "gap open","qstart","qend","sstart","send","E-val","score")
head(blastn312nc)


# blastowanie z prawie automatycznym nazywaniem kolumn
system("blastn -db nt -query G:/NGS/Ciereszko/dane/Swinie2/624_nc_Sus.fasta -out G:/NGS/Ciereszko/dane/Swinie2/results.out -remote -outfmt 7 -max_target_seqs 1")
setwd("G:/NGS/Ciereszko/dane/Swinie2/")
wynik <- read.table("wyblastowane.out", sep = "\t", header = F)
#nazwy kolumn z pozycji #Fields: - pewnie jakos da sie to wyciagnac
colnames(wynik) <- c('query id', 'subject id', '% identity', 'alignment length', 'mismatches', 'gap opens', 'q. start', 'q. end', 's. start', 's. end', 'evalue', 'bit score')
str(wynik)
dim(wynik)
head(wynik)
system("grep Fields results.out")
wynik$V1
#min(grep("# Fields: ", readLines("results.out"))) #jaki jest numer pierwszej linijki, w ktorej wystepuje "# Fields :"
#readLines("results.out")[5] # wczytaj tylko te linijke z calego pliku - nr 5
nazwykolumn <- readLines("results.out")[min(grep("# Fields: ", readLines("results.out")))]
nazwykolumn <- gsub("# Fields: ", "", nazwykolumn)
nazwykolumn <- strsplit(nazwykolumn, split = ", ")
nazwykolumn
colnames(wynik) <- nazwykolumn[[1]]
head(wynik)
nrow(wynik)




#######################################
# query-blast.R - fukcja z neta

library(dplyr)
library(readr)
library(stringr)
library(reutils)

blastn <- function(fasta, db = "refseq_genomic") 
{
  query_name <- paste0(fasta, ".blast.txt")
  comm <- paste("export BLASTDB=/usr/local/share/blast; blastn -query ", fasta,
                "-db", db, 
                "-task blastn",
                "-max_target_seqs 20000",
                "-outfmt '6  sscinames scomnames staxids qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'",
                "-remote",
                "-entrez_query txid6237[ORGN]",
                ">",
                query_name)
  print(comm)
  system(comm)
  
  r <- read_tsv( query_name, col_names = c("Species",
                                           "Name",
                                           "TaxID",
                                           "QueryID",
                                           "SubjectID",
                                           "Percent_Identity",
                                           "Alignment_Length",
                                           "Mismatches",
                                           "Gap_Openings",
                                           "Q.Start",
                                           "Q.End",
                                           "S.Start",
                                           "S.End",
                                           "E",
                                           "Bits") ) %>%
    separate(SubjectID, into = c("name_drop", "gi", "ref_drop","accession"), sep = "\\|", extra = "drop", convert = T) %>%
    dplyr::select(-name_drop, -ref_drop) %>%
    dplyr::mutate(Name = sapply(unlist(efetch(accession, db="nuccore", "docsum")['//Item[@Name="Title"]/text()']), xmlValue) ) %>%
    dplyr::rename(POS_Start = S.Start, POS_End = S.End) %>%
    dplyr::mutate(CHROM = str_match(Name, "chromosome ([A-Za-z0-9])")[,2]) %>%
    dplyr::select(Species, Name, TaxID, QueryID, CHROM, POS_Start, POS_End, accession,  everything()) 
}

results <- blastn("~/Desktop/test.fa")




#####################
# blastx w petli

mainDir <- "c:/path/to/main/dir"
subDir <- "outputDirectory"

if (file.exists(subDir)){
  setwd(file.path(mainDir, subDir))
} else {
  dir.create(file.path(mainDir, subDir))
  setwd(file.path(mainDir, subDir))
  
}

doplikow <- "G:/NGS/Ciereszko/dane/Swinie2/blastx/"
doblasta <- "C:/Program Files/NCBI/blast-2.4.0+/bin/"

seqfas <- read.fasta(file = "G:/NGS/Ciereszko/dane/Swinie2/312_unknow_seq.fasta", seqtype = "DNA", as.string = TRUE, seqonly = F)
# head(seqfas, n=1)
# length(seqfas)
setwd("C:/Program Files/NCBI/blast-2.4.0+/bin/")

for (i in 1:length(seqfas)){
  write.fasta(seqfas[i], names = names(seqfas[i]),  file.out = paste0(doplikow,names(seqfas[i]),".fasta") )
  
  zapytx <- paste0("blastx", 
                  " -db nr",
                  " -query ", doplikow, names(seqfas[i]),".fasta",
                  " -out ", doplikow, names(seqfas[i]),".out",
                  ' -outfmt "6 stitle sscinames scomnames staxids qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"',
                  " -remote",
                  " -max_target_seqs 10")
  paste(zapytx)
  Sys.sleep(5)
  system(zapytx)
}


# poprawka
doplikow <- "G:/NGS/Ciereszko/dane/Swinie2/blastx/powt/"
plikifasta  <- list.files(path = doplikow, pattern = '\\.fasta')
length(plikifasta)
plikifasta[2]

for (i in 1:length(plikifasta)){
  zapytx <- paste0("blastx", 
                   " -db nr",
                   " -query ", doplikow, plikifasta[i],
                   " -out ", doplikow, plikifasta[i],".out",
                   ' -outfmt "6 stitle sscinames scomnames staxids qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"',
                   " -remote",
                   " -max_target_seqs 10")
  paste(zapytx)
  Sys.sleep(25)
  print(i)
  system(zapytx)
}





# laczenie plikow cos nei dziala

setwd("G:/NGS/Ciereszko/dane/Swinie2/blastx/")
setwd(doplikow)
#Using the answer from here [Importing several files and indexing them ] list files with .csv extension - this assumes that the only .csv files in your working directory are the ones you want to read
files  <- list.files(path = doplikow, pattern = '\\.out')
#read files into a list - are there headers?  tables <- lapply(files, read.csv, header = TRUE)
wyniki <- lapply(files, read.table, sep = "\t")
#rbind files
combined.df <- do.call(rbind , wyniki)
# You can then find the mean - find which columns are numeric
s <- sapply(combined.df, is.numeric)
# find the mean of numeric variables
colMeans(combined.df[s])



  wynik <- read.table(paste0(doplikow, names(seqfas[i]),".out"), sep = "\t")
  write(wynik, file = paste0(doplikow, "_wynik.out"), append = T, sep = "\n")
  
  colnames(blastn312nc) <- c("title", "orgn", "orgn com", "taxid", "q seq id", "s seq id", "p ident", "length", "mismatch", "gap open","qstart","qend","sstart","send","E-val","score")
head(blastn312nc)


  
#################################
### BLASTp Anisakisa
setwd("C:/Program Files/NCBI/blast-2.4.0+/bin/")
setwd("G:/NGS/ElaLopienska-Biernat/")
doblasta <- "\"C:/Program Files/NCBI/blast-2.4.0+/bin/\""
doplikow <- "G:/NGS/ElaLopienska-Biernat/"


zapytp <- paste0("blastp", 
                 " -db nr",
                 #" -query ", choose.files(),
                 " -query G:/NGS/ElaLopienska-Biernat/transport_protein_selected.fasta",
                 " -out G:/NGS/ElaLopienska-Biernat/blastp.out",
                 ' -outfmt "7 qseqid sseqid qstart qend sstart send evalue bitscore length pident qcovs qlen"',
                 " -remote",
                 " -max_target_seqs 1",
                 " -max_hsps 1",
                 " -entrez_query txid6231[ORGN]")

zaphtmlp <- paste0("blastp", 
                 " -db nr",
                 #" -query ", choose.files(),
                 " -query G:/NGS/ElaLopienska-Biernat/transport_protein_selected.fasta",
                 " -out G:/NGS/ElaLopienska-Biernat/blastp_kom.html",
                 " -remote",
                 " -num_descriptions 3",
                 " -num_alignments 3",
                 " -html",
                 " -entrez_query txid6231[ORGN]")
#paste(zapytp)
system(zapytp)
system(zaphtmlp)


wynik <- read.table("G:/NGS/ElaLopienska-Biernat/blastp.out", sep = "\t")
#min(grep("# Fields: ", readLines("results.out"))) #jaki jest numer pierwszej linijki, w ktorej wystepuje "# Fields :"
#readLines("results.out")[5] # wczytaj tylko te linijke z calego pliku - nr 5
nazwykolumn <- readLines("G:/NGS/ElaLopienska-Biernat/blastp.out")[min(grep("# Fields: ", readLines("G:/NGS/ElaLopienska-Biernat/blastp.out")))]
nazwykolumn <- gsub("# Fields: ", "", nazwykolumn)
nazwykolumn <- strsplit(nazwykolumn, split = ", ")
nazwykolumn
colnames(wynik) <- nazwykolumn[[1]]
head(wynik)
nrow(wynik)

write.csv(wynik, "G:/NGS/ElaLopienska-Biernat/blastp_out.csv")
