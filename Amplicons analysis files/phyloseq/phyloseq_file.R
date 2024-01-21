library(qiime2R)
library("phyloseq")
library("ggplot2")      # graphics
library("readxl") # necessary to import the data from Excel file
library("plyr")
library("dplyr")        # filter and reformat data frames
library("tibble")
library("microbiomeSeq")
library("vegan")
library("ape")
library("dunn.test")

dir()
otu_table <- read.table("rarefied-feature-table.tsv",
                       sep = "\t", header = TRUE, row.names = 1)
class(otu_table)
head(otu_table)
otu_table <- as.matrix(otu_table)
class(otu_table)

taxa <- read.table("taxonomy.tsv",
                  sep = "\t", header = TRUE, row.names = 1)
taxa <- as.matrix(taxa)
head(taxa)
class(taxa)

metadata <- read.table("sample-metadata.tsv",
                      sep = "\t", header = TRUE, row.names = 1)
head(metadata)

trefile = read.tree(file = "tree.nwk")
print(trefile)

OTU <- otu_table(otu_table, taxa_are_rows=TRUE)
TAX = tax_table(taxa)
sample = sample_data(metadata)
phyloseq_obj <- phyloseq(OTU, TAX, sample)
print(phyloseq_obj)



#############################################plot_richness###############################################
?plot_richness


plot_richness(phyloseq_obj, measures=c("Shannon","simpson", "Observed"), 
              x="samples", shape = "Type" , color="Group")


##############################################################################################################


## Kruskal wallis test

alpha_Supraglacial1 <- subset_samples(phyloseq_obj, Group=="Supraglacial-1")
alpha_Supraglacial2 <- subset_samples(phyloseq_obj, Group=="Supraglacial-2")
alpha_Proglacial1 <- subset_samples(phyloseq_obj, Group=="Proglacial site-1")
alpha_Proglacial2 <- subset_samples(phyloseq_obj, Group=="Proglacial site-2")

all_alpha <- merge_phyloseq(alpha_Supraglacial1, alpha_Supraglacial2,
                            alpha_Proglacial1, alpha_Proglacial2 )
all_alpha
alpha_observed <- estimate_richness(all_alpha, measures = "Observed")
alpha_Shannon <- estimate_richness(all_alpha, measures = "Shannon")
alpha_Chao1 <- estimate_richness(all_alpha, measures = "Chao1")
alpha_simpson <- estimate_richness(all_alpha, measures = "simpson")
alpha_ACE <- estimate_richness(all_alpha, measures = "ACE")
alpha_InvSimpson <- estimate_richness(all_alpha, measures = "InvSimpson")

?kruskal.test
alpha.stats <- cbind(alpha_observed, sample_data(all_alpha))
kruskal.test(Observed~Group, alpha.stats)

?dunn.test
dunn.test(alpha.stats$Observed, alpha.stats$Group, 
          method = "bonferroni")

alpha.stats <- cbind(alpha_Chao1, sample_data(all_alpha))
kruskal.test(Chao1~Group, alpha.stats)

dunn.test(alpha.stats$Chao1, alpha.stats$Group, 
          method = "bonferroni")

alpha.stats <- cbind(alpha_Shannon, sample_data(all_alpha))
kruskal.test(Shannon~Group, alpha.stats)

dunn.test(alpha.stats$Shannon, alpha.stats$Group, 
          method = "bonferroni")

alpha.stats <- cbind(alpha_simpson, sample_data(all_alpha))
kruskal.test(Simpson~Group, alpha.stats)

dunn.test(alpha.stats$Simpson, alpha.stats$Group, 
          method = "bonferroni")

alpha.stats <- cbind(alpha_ACE, sample_data(all_alpha))
kruskal.test(ACE~Group, alpha.stats)

dunn.test(alpha.stats$ACE, alpha.stats$Group, 
          method = "bonferroni")

alpha.stats <- cbind(alpha_InvSimpson, sample_data(all_alpha))
kruskal.test(InvSimpson~Group, alpha.stats)

dunn.test(alpha.stats$InvSimpson, alpha.stats$Group, 
          method = "bonferroni")
