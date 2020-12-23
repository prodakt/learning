# wczyt. bib ----
BiocManager::install("openPrimeR")
#install.packages("OpenPrimeR")
library("openPrimeR")
library(ggplot2)

setwd("D:/SkryptyR/OpenPrimeR/OpenPrimeR/")


# jpj --------
# wczytywanie
Seq <- read_templates("Exon_Exon.fasta")
head(Seq)

structure <- c( "accession","organism", "gene", "moltype" )
Seq1 <- read_templates("Exon_Exon.fasta", structure, delim = ",")
Seq1$Sequence_Length



# okreslanie miejsc przylaczanie startewor
# Zakres binding regions
B_region <- assign_binding_regions(Seq1, fw = c(20,45), rev= c(20,45))

frs <- 20
fre <- 45
rvs <- 55
rve <- 80
# albo
frs <- 20
fre <- 45
product_min <- 50
product_max <- 300
rvs <- fre + product_min
rve <- fre + product_max + 25

B_region <- assign_binding_regions(Seq1, fw = c(frs,fre), rev= c(Seq1$Sequence_Length - rve,Seq1$Sequence_Length - rvs))

Bin_region <- assign_binding_regions(Seq1, fw =c(60,90), rev= c(15,100))
Bin_region

# wczytywanie ustawien
design_settings <- read_settings("My_settings.xml")
design_settings
str(settings)

# projektowanie primerow
optimal.primers <- design_primers(Bin_region, mode.directionality = "both", settings = design_settings)

# plots

plot_primer_binding_regions(optimal.primers$opti, Bin_region)
plot_constraint_fulfillment(optimal.primers$opti, design_settings)
plot_primer(optimal.primers$opti, Seq1)


pprimer <- plot_primer(optimal.primers$opti, Seq1)
head(pprimer)
pprimer + theme_bw() + scale_x_continuous(breaks = c(-30, 2500))

pprimer + theme_bw() + scale_x_continuous(breaks = c(-30, Seq1$Sequence_Length))
pprimer + theme_bw() + scale_x_continuous(breaks = c(-30, Seq1$Allowed_End_rev))
pprimer + theme_bw() + scale_x_continuous(breaks = c(300, 1500)) + 

# praca na zapisanych i wczytanych danych
# zapisywanie
head(optimal.primers)
write_primers(optimal.primers$opti, "optimal_primers.out")


#wczytywanie
Primers <- read_primers("optimal_primers.out")
# walidacja
constraint.df <- check_constraints(Primers, Bin_region, design_settings, active.constraints = names(constraints(design_settings)))
constraint.df
# wykresy
plot_primer_binding_regions(constraint.df, Bin_region)
plot_constraint_fulfillment(constraint.df, design_settings)
plot_primer(constraint.df, Seq1)





wykresy_z_danych_na_dysku <- function(prajmery, fasta, fwu, rew, struktura){
  Seq1 <- read_templates(fasta, struktura, delim = ",")
  #wczytywanie
  Primers <- read_primers(prajmery)
  Bin_region <- assign_binding_regions(Seq1, fw = fwu, rev= rew)
  # walidacja
  constraint.df <- check_constraints(Primers, Bin_region, design_settings, active.constraints = names(constraints(design_settings)))
  constraint.df
  # wykresy
  plot_primer_binding_regions(constraint.df, Bin_region)
  plot_constraint_fulfillment(constraint.df, design_settings)
  plot_primer(constraint.df, Seq1)
}
wykresy_z_danych_na_dysku(prajmery = "optimal_primers.out", fasta = "Exon_Exon.fasta",fwu = c(60,90), rew = c(15,100), 
                          struktura = structure <- c( "accession","organism", "gene", "moltype" ))






# ////////////////////////////////////////////////////////////////////////
# Pana Mariusza ------
directory <- "c:/Users/mariu/Desktop/R" # Okreslenie lokalizacji pliku
Sekwencja <- file.path(directory, "Exon_Exon.fasta")#okreslanie nazwy pliku w lokalizacji
Seq <- read_templates(Sekwencja) #przypisanie funkcji odczytania zaladowanej sekwencji
structure <- c( "accession","organism", "gene", "moltype" )#Okreslenie struktury naglowka
Seq1 <- read_templates(Sekwencja, structure, delim = ",")
#Zastosowanie stworzonej struktury,
#do zaladowanej matrycy. Scislej mam na mysli to ze kazdy plik fasta z sekwencja ma naglowek
# w ktorym zawarte sa rozne informacje-accession number etc. Tu okreslam w jakiej kolejnosci,
#dane informacje wystepuja i czym sa rodzielone(w tym przypadku przecinkiem)
#Nie jest to w zaden sposob konieczne do przeprowadzenia analizy ale ulatwia polapanie sie

Bin_region <- assign_binding_regions(Seq1, fw =c(60,90), rev= c(15,100))
Bin_region
#
#okreslanie miejsc w ktorych maja przylaczyc sie primery

settings.xml<- file.path(directory, "My_settings.xml")#okreslenie lokalizacji
#i nazwy pliku z ustawieniam i odczytanie pliku
settings <-read_settings(settings.xml)
design.settings<- settings#przypisanie zmiennej do zmieniania ustawien
constraints(design.settings)[["primer_length"]] <- c("min" = 18, "max" = 30)#funkcja 
#okreslajaca mozliwa dlugosc starterow
constraints(design.settings)[["gc_ratio"]] <- c( "min" = 0.4, "max"= 0.6)# funkcja okreslajaca
#mozliwa zawartosc GC w primerach
constraints(design.settings)[["gc_clamp"]] <- c( "min" = 0, "max" = 3)#funkcja okreslajaca
#maksymalna ilosc GC na koncu 3'(ostatnich 5 nukletydach)
constraints(design.settings)[["melting_temp_range"]] <-c(
  "min"= 59, "max" = 61)#funkcja okreslajaca melting temperature primerow

out.file <- tempfile("settings", fileext = ".xml")
write_settings(design.settings, out.file)# ta funkcja i poprzednia to utworzenie pliku z ustawieniami


optimal.primers <- design_primers(Bin_region, mode.directionality = "both",
                                  settings = design.settings)
#funkcja do tworzenia primerow. W argumentach wpisana zmienna z okreslonymi binding regions,
#drugi argument to okreslenie czy ma byc stworzony tylko forward primer "fw" czy tylko reverse "rev", 
#czy jak w tym przypadku oba "both". Ostatni argument to zastosowanie ustawien)

out.file <- tempfile("my_primers", fileext = ".fasta")
write_primers(optimal.primers$opti, out.file)# zapisanie pliku ze stworzonymi primerami.

Primer<- file.path(directory, "Primer1.fasta")# wczytanie pliku ze stworzony primerami
Primers<-read_primers(Primer)#funkcja do oczytania primerow
constraint.df <- check_constraints(Primers, Bin_region, 
                                   design.settings, active.constraints = names(constraints(settings)))
#Nie potrafie do konca wytlumaczyc tej funkcji. Jest to swego rodzaju sprawdzanie wlasciwosc biochemicznych
#stworzonych primerow, konieczne przy pozniejszym tworzeniu wykresow.Po kliknieciu 
# w pakiet openprimer w okienku packages, a nastepnie kliknieciu w user guides... otwiera sie
# swego rodzaju tutorial z funkcjami open primera,i ta funckaj i jej zastosowanie jest niejako tlumaczone w dziale "
#Evaluation of biochemical constraints. Sam do konca nie wiem czy obecnie uzywam jej poprawnie.
#

plot_primer_binding_regions(constraint.df, Bin_region)#wykres wskazujacy binding regions i pozycje primerow

plot_constraint_fulfillment(constraint.df, design.settings)#wykres pokazujacy czy wytworzone primery
# spelniaja poszczegolne ustawienia okreslone w design.settings


plot_primer(constraint.df, Seq1)#wykres wskazujacy kierunki primerow
