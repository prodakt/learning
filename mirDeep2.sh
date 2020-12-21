# miRDeep2
# from
# https://github.com/rajewsky-lab/mirdeep2/blob/master/TUTORIAL.md

# budowanie indeksu ----
bowtie-build cel_cluster.fa cel_cluster
bowtie-build Gallus_gallus.GRCg6a.dna.toplevel.fa Gallus6a_bowtie

perl -plane 's/\s+.+$//' < Gallus_gallus.GRCg6a.dna.toplevel.fa > GRCg6a_NoSpace.fa
# trzeba zostawic PELNA nazwe miRNA, ale np. z lacznikiem "_", a nie usuwac wszsytko po spacji, bo potem
# w tabeli ekspresji nazwy miRNA nie sa unikatowe i jest problem.
# tymczasowo w DESeq zrobilem laczenie nazwy miRNA z prekursorem
sed '/^[^>]/s/[R|Y|W|S|M|K|H|B|V|D]/N/g' GRCg6a_NoSpace.fa > GRCg6a_NoSpace_ATGCN.fa

# mapowanie ----
# adapter sequence from https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences-1000000002694-14.pdf
mapper.pl pliki.txt -d -h -c -i -j -k TGGAATTCTCGGGTGCCAAGG  -l 18 -m -p /home/m1/NGS/db/gallus/Gallus6a_bowtie -s reads_collapsed.fa -t reads_collapsed_vs_genome.arf -v

# mapper.pl reads.fa -c -j -k TCGTATGCCGTCTTCTGCTTGT  -l 18 -m -p cel_cluster  -s reads_collapsed.fa -t reads_collapsed_vs_genome.arf -v
# mapper.pl config.txt -d -c -i -j -l 18 -m -p genome_index -s reads.fa -t reads_vs_genome.arf

# zliczanie
quantifier.pl -p /home/m1/NGS/db/miRNA/hairpin.fa -m /home/m1/NGS/db/miRNA/mature.fa -r reads_collapsed.fa -t gga -y 11_20

# identyfikacja 
# miRDeep2.pl reads_collapsed.fa /home/m1/NGS/db/gallus/Gallus_gallus.GRCg6a.dna.toplevel.fa reads_collapsed_vs_genome.arf mature_ref_this_species.fa mature_ref_other_species.fa precursors_ref_this_species.fa -t Chicken 2> report.log

# po quanifier
# miRDeep2.pl reads_collapsed.fa /home/m1/NGS/db/gallus/GRCg6a_NoSpace.fa reads_collapsed_vs_genome.arf /home/m1/NGS/db/miRNA/gga_mat.fa /home/m1/NGS/db/miRNA/cli_mar.fa  /home/m1/NGS/db/miRNA/gga_mat.fa -t Chicken -q expression_analyses/expression_analyses_11_20/miRBase.mrd 2> report2.log

miRDeep2.pl reads_collapsed.fa /home/m1/NGS/db/gallus/GRCg6a_NoSpace_ATGCN.fa reads_collapsed_vs_genome.arf /home/m1/NGS/db/miRNA/gga_mat_NoSpace.fa /home/m1/NGS/db/miRNA/cli_mar_NoSpace.fa  /home/m1/NGS/db/miRNA/gga_pre_NoSpace.fa -t Chicken -q expression_analyses/expression_analyses_11_20/miRBase.mrd 2> report6.log


# z loga 1
# miRDeep2.pl reads.fa genome.fa reads_vs_genome.arf mautre_ref_miRNAs.fa mature_other_miRNAs.fa  hairpin_ref_miRNAs
# or
# miRDeep2.pl reads.fa genome.fa reads_vs_genome.arf none none none


# ------------------------------------------------------------
# wynik mapper
#desc   total   mapped  unmapped        %mapped %unmapped
total: 142974679        127105345       15869334        88.901  11.099
ct1: 13076209   11738357        1337852 89.769  10.231
ct2: 12866460   11882178        984282  92.350  7.650
ct3: 16698096   15363341        1334755 92.007  7.993
ct4: 15632405   14338729        1293676 91.724  8.276
ct5: 17523621   16469533        1054088 93.985  6.015
vr1: 15958970   14295207        1663763 89.575  10.425
vr2: 10518628   8905758 1612870 84.667  15.333
vr3: 13397325   11201967        2195358 83.613  16.387
vr4: 13572414   11505933        2066481 84.774  15.226
vr5: 13730551   11404342        2326209 83.058  16.942


# ///////////////////
# (base) jasiu@ngs3:/home/m4/NGS/JoannaPIWET/miRNA/rr$ miRDeep2.pl reads_collapsed.fa /home/m1/NGS/db/gallus/GRCg6a_NoSpace_ATGCN.fa reads_collapsed_vs_genome.arf /home/m1/NGS/db/miRNA/gga_mat_NoSpace.fa /home/m1/NGS/db/miRNA/cli_mar_NoSpace.fa  /home/m1/NGS/db/miRNA/gga_pre_NoSpace.fa -t Chicken -q expression_analyses/expression_analyses_11_20/miRBase.mrd 2> report6.log

#####################################
#                                   #
# miRDeep2.0.1.2                    #
#                                   #
# last change: 22/01/2019           #
#                                   #
#####################################

miRDeep2 started at 14:37:27


#Starting miRDeep2
#testing input files
#Quantitation of known miRNAs in data
#parsing genome mappings
#excising precursors
#preparing signature
#folding precursors
#computing randfold p-values
#running miRDeep core algorithm
#running permuted controls
#doing survey of accuracy
#producing graphic results


miRDeep runtime:

started: 14:37:27
ended: 19:36:05
total:4h:58m:38s