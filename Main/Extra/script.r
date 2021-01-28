library(phyloseq)
library(ggplot2)
library(dplyr)

ps1 <- readRDS("D:/GitHub/Seagrass Diversity/seagrass-diversity/Main/Output/Merged/final_phyloseq_object_noOutlier.RDS")

## Relative abundance barchart ----

# Creating a dataframe of the otu table which also includes the Location and Structure variables
combineddf <- as.data.frame(ps1@otu_table)
combineddf$Location <- ps1@sam_data$Location
combineddf$Structure <- ps1@sam_data$Structure.DNA.Extracted.from

# Creating a phyloseq object where samples are grouped and merged according to Location and Structure
location_structure_otu <- combineddf %>%
  group_by(Location, Structure) %>%
  summarise_all(.funs=sum)
location_structure_otu <- as.data.frame(location_structure_otu)
location_structure_otu <- location_structure_otu[,c(-1,-2)]
location_structure_meta <- data.frame(Location = c("Cyrene", "Cyrene", "Cyrene", "Cyrene", "Merambong Shoal", "Merambong Shoal", "Merambong Shoal", "Merambong Shoal", "Perhentian Island", "Perhentian Island", "Perhentian Island", "Perhentian Island", "Port Dickson", "Port Dickson", "Port Dickson", "Port Dickson", "Semakau", "Semakau", "Semakau", "Semakau", "Sentosa", "Sentosa", "Sentosa", "Sentosa"), 
                                      Structure = c("Leaf", "Rhizome", "Root", "Sediment", "Leaf", "Rhizome", "Root", "Sediment", "Leaf", "Rhizome", "Root", "Sediment", "Leaf", "Rhizome", "Root", "Sediment", "Leaf", "Rhizome", "Root", "Sediment", "Leaf", "Rhizome", "Root", "Sediment"))

location_structure_ps <- phyloseq(otu_table(location_structure_otu, taxa_are_rows=FALSE), 
                                  sample_data(location_structure_meta), 
                                  tax_table(ps1@tax_table),
                                  refseq(ps1@refseq))


# Calculating Relative abundances and plotting bar plots according to Location and Structure
ps_Location_Structure_ra <- transform_sample_counts(location_structure_ps, function(otu){otu/sum(otu)})

# Defining function to make bar charts without black lines separating samples. Based on phyloseq function "plot_bar".
simple_plot_bar = function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, 
                            facet_wrap = NULL) {
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack")
  p = p + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p = p + labs(y = "Relative Abundance")
  p = p + guides(guide_legend(ncol = 1))
  if (!is.null(facet_wrap)) {
    p <- p + facet_wrap(facet_wrap, nrow = 1)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}


# Making stacked bar charts for relative abundance of taxa

# According to Phylum
ra_Phylum_barplot_location_grouped_by_structure <- simple_plot_bar(ps_Location_Structure_ra, x="Location", fill="Phylum", facet_wrap = "Structure") + theme(legend.text=element_text(size=rel(0.8)), legend.key.size= unit(0.8, "line")) + guides(fill = guide_legend(nrow = 25))
ggsave(ra_Phylum_barplot_location_grouped_by_structure, filename = "./Relative Abundance of Phylum by Location and Structure - Bar Plot.pdf", dpi=300, width = 12, height = 10)

## Venn Diagram ----
library(VennDiagram)

# Creating Presence Absence Table according to structure
df <- as.data.frame(ps1@otu_table)
df[df>0] <- 1
df$Structure <- ps1@sam_data$Structure.DNA.Extracted.from
otutable <- df %>%
  group_by(Structure) %>%
  summarise_all(.funs=mean)
otutable <- as.data.frame(otutable)
row.names(otutable) <- c("Samples_Le", "Samples_Rh", "Samples_Ro", "Samples_Se")
otutable <- otutable[,-1]

# Creating a phyloseq object where samples are grouped and merged according to Structure
structure_otu <- combineddf[,-1586] %>%
  group_by(Structure) %>%
  summarise_all(.funs=sum)
structure_otu <- as.data.frame(structure_otu)
row.names(structure_otu) <- c("Samples_Le", "Samples_Rh", "Samples_Ro", "Samples_Se")
structure_otu <- structure_otu[,-1]
structure_meta <- data.frame(Structure = c("Leaf", "Rhizome", "Root", "Sediment"))
row.names(structure_meta) <- NULL
row.names(structure_meta) <- c("Samples_Le", "Samples_Rh", "Samples_Ro", "Samples_Se")

structure_ps <- phyloseq(otu_table(structure_otu, taxa_are_rows=FALSE), 
                         sample_data(structure_meta), 
                         tax_table(ps1@tax_table),
                         refseq(ps1@refseq))


otutable <- as.data.frame(otu_table(structure_ps))
otutable[otutable>0] <- 1
filterVector <- c()
for (i in 1:ncol(otutable)){
  x <- sum(otutable[,i])/nrow(otutable) >= 0.95
  x <- sum(otutable[,i])/nrow(otutable)
  filterVector <- c(filterVector, x)
}
patable <- t(otutable)

#Converting to presence/absence table
otutable <- as.data.frame(structure_otu)
otutable[otutable>0] <- 1
patable <- as.data.frame(t(otutable))

# Calculating Size of Sets
size_Le <- sum(patable[,"Samples_Le"])
size_Rh <- sum(patable[,"Samples_Rh"])
size_Ro <- sum(patable[,"Samples_Ro"])
size_Se <- sum(patable[,"Samples_Se"])

# Calculating size of overlap between sets
overlap_Le_Rh <- as.numeric(sum(patable[,"Samples_Le"] == patable[,"Samples_Rh"] & patable[,"Samples_Le"] > 0))
overlap_Le_Ro <- as.numeric(sum(patable[,"Samples_Le"] == patable[,"Samples_Ro"] & patable[,"Samples_Le"] > 0))
overlap_Le_Se <- as.numeric(sum(patable[,"Samples_Le"] == patable[,"Samples_Se"] & patable[,"Samples_Le"] > 0))
overlap_Rh_Ro <- as.numeric(sum(patable[,"Samples_Rh"] == patable[,"Samples_Ro"] & patable[,"Samples_Rh"] > 0))
overlap_Rh_Se <- as.numeric(sum(patable[,"Samples_Rh"] == patable[,"Samples_Se"] & patable[,"Samples_Rh"] > 0))
overlap_Ro_Se <- as.numeric(sum(patable[,"Samples_Ro"] == patable[,"Samples_Se"] & patable[,"Samples_Ro"] > 0))

overlap_Le_Rh_Ro <- as.numeric(sum(patable[,"Samples_Le"] == patable[,"Samples_Rh"] & patable[,"Samples_Le"] == patable[,"Samples_Ro"] & patable[,"Samples_Le"] > 0))
overlap_Le_Rh_Se <- as.numeric(sum(patable[,"Samples_Le"] == patable[,"Samples_Rh"] & patable[,"Samples_Le"] == patable[,"Samples_Se"] & patable[,"Samples_Le"] > 0))
overlap_Le_Ro_Se <- as.numeric(sum(patable[,"Samples_Le"] == patable[,"Samples_Ro"] & patable[,"Samples_Le"] == patable[,"Samples_Se"] & patable[,"Samples_Le"] > 0))
overlap_Rh_Ro_Se <- as.numeric(sum(patable[,"Samples_Rh"] == patable[,"Samples_Ro"] & patable[,"Samples_Rh"] == patable[,"Samples_Se"] & patable[,"Samples_Rh"] > 0))

overlap_Le_Rh_Ro_Se <- as.numeric(sum(patable[,"Samples_Le"] == patable[,"Samples_Rh"] & patable[,"Samples_Rh"] == patable[,"Samples_Ro"] & patable[,"Samples_Ro"] == patable[,"Samples_Se"] & patable[,"Samples_Le"] > 0))

# Creating Venn Diagram with 4 sets
venn_diagram <- draw.quad.venn(size_Le, size_Rh, size_Ro, size_Se, 
                               overlap_Le_Rh, overlap_Le_Ro, overlap_Le_Se, overlap_Rh_Ro, overlap_Rh_Se, overlap_Ro_Se,
                               overlap_Le_Rh_Ro, overlap_Le_Rh_Se, overlap_Le_Ro_Se, overlap_Rh_Ro_Se, 
                               overlap_Le_Rh_Ro_Se,
                               category = c("Leaf", "Rhizome", "Root", "Sediment"), 
                               lwd = 0.5, fill = col_pal_struct, alpha = 0.6, cex = 1, cat.cex = 1.3)

ggsave(venn_diagram, filename = "./Output/Venn_Diagram_by_structure.pdf", dpi = 300, width = 9, height = 5)
