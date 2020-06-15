# Loading Packages
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(broom)
#library(stargazer)
#library(sjPlot)
library(RColorBrewer)
library(RgoogleMaps)
library(metagMisc)
library(shiny)
library(viridis)
library(svglite)
library(gdtools)
library(VennDiagram)
library(indicspecies)
library(geosphere)
library(ecodist)
library(SpiecEasi)
library(Matrix)
library(reshape2)
library(tidygraph)
library(tidyverse)
library(igraph)
library(WGCNA)
library(corncob)
source("https://raw.githubusercontent.com/genomewalker/osd2014_analysis/master/osd2014_16S_asv/lib/graph_lib.R")



# Color palette
pal = c("#6b5456","#ec8d1b","#6abf2a","#8b53b7","#70acbe","#01c95b","#c00014","#31332f","#f7d000","#abba00")
col_structure <- viridis(4)
col_location <- viridis(6)
names(col_structure) <- c("Leaf", "Rhizome", "Root", "Sediment")
names(col_location) <- c("Cyrene", "Merambong Shoal", "Port Dickson","Semakau", "Sentosa", "Perhentian Island")

theme_set(theme_bw())

# load processed data
ps = readRDS(file = "D:/GitHub/Seagrass Diversity/seagrass-diversity/Main/Output/Merged/clean_phyloseq_object.RDS")

# quick peek at ps object
sample_names(ps)
rank_names(ps)
sample_variables(ps)
otu_table(ps)[1:5, 1:5]
tax_table(ps)[1:5, 1:6]

# Change sequences to unique IDs to make viewing easier
seqs_16S = taxa_names(ps)
names(seqs_16S) <- 1:length(seqs_16S) # save seqs and IDs combination
taxa_names(ps) <- names(seqs_16S)
tax_table(ps)[1:5, 1:6]

# quick exploratory look at top families by location
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
ps.islands = prune_samples(ps.top20@sam_data$Location != "",ps.top20)
top20_barplot_location <- plot_bar(ps.top20, x="Location", fill="Family")
ggsave(top20_barplot_location, filename = "../Output/Analysis/Top 20 Bacterial Families by Location - Bar Plot.png", dpi=300, width = 12, height = 10)

# quick exploratory look at top families by structure
ps.structure = prune_samples(ps.top20@sam_data$Structure.DNA.Extracted.from != "",ps.top20)
top20_barplot_structure <- plot_bar(ps.top20, x="Structure.DNA.Extracted.from", fill="Family")
ggsave(top20_barplot_structure, filename = "../Output/Analysis/Top 20 Bacterial Families by Structure - Bar Plot.png", dpi=300, width = 12, height = 10)

# Visualizing Locations
meta <- ps@sam_data
par(pty="s")
meta$lat <- as.numeric(sapply(strsplit(meta$GPS.Coordinates, "N"), `[`, 1))
meta$lon <- sapply(strsplit(meta$GPS.Coordinates, " "), `[`, 2)
meta$lon <- as.numeric(sapply(strsplit(meta$lon, "E"), `[`, 1))
svg(filename = "../Output/Analysis/Full Map.svg")
PlotOnStaticMap(lat = unique(meta$lat), lon = unique(meta$lon), size = c(1000,1000), zoom = 150, 
                               cex = 2, pch = 21, col = "black", bg = col_location, lwd = 2, FUN = points, add = F, GRAYSCALE = TRUE) 
legend(170, 450, legend = names(col_location), col = "black", pt.bg = col_location, pch = 21, cex = 1, pt.cex = 2, box.lty = 1)
dev.off()
svg(filename = "../Output/Analysis/SG Map.svg")
PlotOnStaticMap(lat = unique(meta$lat)[c(-3,-6)], lon = unique(meta$lon)[c(-3,-6)], size = c(1000,1000), zoom = 150, 
                cex = 2, pch = 21, col = "black", bg = col_location[c(-3,-6)], lwd = 2, FUN = points, add = F, GRAYSCALE = TRUE) 
legend(170, 450, legend = names(col_location)[c(-3,-6)], col = "black", pt.bg = col_location[c(-3,-6)], pch = 21, cex = 1, pt.cex = 2, box.lty = 1)
dev.off()


# Removing less prevalent sequences
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

# compute mean and total prevalence
phylum_prev = plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
family_prev = plyr::ddply(prevdf, "Family", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Look at abundance vs prevalence
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps, "Phylum"))
abundance_vs_prevalence <- ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.025, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
ggsave(abundance_vs_prevalence, filename = "../Output/Analysis/abundance_vs_prevalence_Phylum.png", dpi=300, width = 12, height = 10)

#  Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps)

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps1 = prune_taxa(keepTaxa, ps)

# Remove samples not associated with an island
ps1 = prune_samples(ps1@sam_data$Location != "", ps1)

#saveRDS(ps1, file = "D:/GitHub/Seagrass Diversity/seagrass-diversity/Main/Output/Merged/final_phyloseq_object.RDS")
#save(ps1, file = "D:/GitHub/Seagrass Diversity/seagrass-diversity/Main/Output/Merged/final_phyloseq_object.RData")
#shiny::runGitHub("shiny-phyloseq","joey711")


# Removing Outliers: Samples with only one ASV in them found from Shannon Diversity data as they are not to be expected
div_shannon <- data.frame(Location = ps1@sam_data$Location,
                  Structure = ps1@sam_data$Structure.DNA.Extracted.from,
                  Shannon = diversity(otu_table(ps1)))
ps1_noOutlier <- prune_samples(div_shannon$Shannon != 0, ps1)
#saveRDS(ps1_noOutlier, file = "D:/GitHub/Seagrass Diversity/seagrass-diversity/Main/Output/Merged/final_phyloseq_object_noOutlier.RDS")

# Note: All work continued does not include the 3 outliers









#### Analysis ####
ps1 <- readRDS("D:/GitHub/Seagrass Diversity/seagrass-diversity/Main/Output/Merged/final_phyloseq_object_noOutlier.RDS")
# Normalize (relative abundance) ####
ps1ra <- transform_sample_counts(ps1, function(otu){otu/sum(otu)})





### Rarefaction Curves
grp <- factor(ps1@sam_data$Structure.DNA.Extracted.from)
cols <- col_structure[grp]
rarefaction_curve <- rarecurve(ps1@otu_table, step = 20, col = col, label = TRUE)
Nmax <- sapply(rarefaction_curve, function(x) max(attr(x, "Subsample")))
Smax <- sapply(rarefaction_curve, max)
svg(filename = "../Output/Analysis/Rarefaction Curves.svg", )
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Number of Sequences",
     ylab = "Number of ASVs", type = "n",
     main = "Rarefaction Curves")
for (i in seq_along(rarefaction_curve)) {
  N <- attr(rarefaction_curve[[i]], "Subsample")
  lines(N, rarefaction_curve[[i]], col = cols[i])
}
legend(155000, 1000, legend = names(col_structure), col = col_structure, lty = 1, cex = 0.8, box.lty = 1)
dev.off()





# Plot abundance function ####
plot_abundance = function(physeq,title = "", 
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = physeq
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "Structure.DNA.Extracted.from",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}

# abundance before and after normalization
plotBefore = plot_abundance(ps1,"") + ggtitle("Before Normalization") + theme(axis.text.x = element_text(angle = 75, hjust = 1))
ggsave(plotBefore, filename = "../Output/Analysis/Abundance_before_norm.png", dpi= 300, width = 12, height = 10)
plotAfter = plot_abundance(ps1ra,"") + ggtitle("After Relative") + theme(axis.text.x = element_text(angle = 75, hjust = 1))
ggsave(plotAfter, filename = "../Output/Analysis/Abundance_after_norm.png", dpi=300, width = 12, height = 10)

# Make final figure of relative abundance of phyla by Location ####
# Arbitrary subset, based on Phylum, for plotting

mphyseq = psmelt(ps1ra)
mphyseq <- subset(mphyseq, Abundance > 0)

relabundplot = ggplot(data = mphyseq, mapping = aes_string(x = "Structure.DNA.Extracted.from",y = "Abundance",
                                                           color = "Phylum")) +
  geom_violin(fill = NA) +
  geom_point(size = 1, alpha = 0.3,
             position = position_jitter(width = 0.3)) +
  facet_wrap(facets = "Phylum") + scale_y_log10()+
  theme_bw() + theme(legend.position="none") + 
  labs(x="Structure.DNA.Extracted.from",y="Relative Abundance") +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 8)) 

ggsave(relabundplot, filename = "../Output/Analysis/Phylum_relabund_by_Structure.png",dpi=300)





### Relative Abundance Bar Plots by Location and Structure

# Creating a dataframe of the otu table which also includes the Location and Structure variables
combineddf <- as.data.frame(ps1@otu_table)
combineddf$Location <- ps1@sam_data$Location
combineddf$Structure <- ps1@sam_data$Structure.DNA.Extracted.from

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

# Creating a phyloseq object where samples are grouped and merged according to Location
location_otu <- combineddf[,-1587] %>%
  group_by(Location) %>%
  summarise_all(.funs=sum)
location_otu <- as.data.frame(location_otu)
row.names(location_otu) <- c("Samples_CY", "Samples_MShoal", "Samples_Per", "Samples_PD", "Samples_Sem", "Samples_Sent")
location_otu <- location_otu[,-1]
location_meta <- data.frame(Location = c("Cyrene", "Merambong Shoal", "Perhentian Island", "Port Dickson", "Semakau", "Sentosa"))
row.names(location_meta) <- NULL
row.names(location_meta) <- c("Samples_CY", "Samples_MShoal", "Samples_Per", "Samples_PD", "Samples_Sem", "Samples_Sent")

location_ps <- phyloseq(otu_table(location_otu, taxa_are_rows=FALSE), 
                         sample_data(location_meta), 
                         tax_table(ps1@tax_table),
                         refseq(ps1@refseq))

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
ps_Location_ra <- transform_sample_counts(location_ps, function(otu){otu/sum(otu)})
ps_Structure_ra <- transform_sample_counts(structure_ps, function(otu){otu/sum(otu)})
ps_Location_Structure_ra <- transform_sample_counts(location_structure_ps, function(otu){otu/sum(otu)})

# Defining function to make bar charts without black lines separating samples. Based on phyloseq function "plot_bar".
simple_plot_bar = function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, 
                        facet_wrap = NULL) {
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack")
  p = p + theme_bw() + theme(axis.text=element_text(size=15), axis.title=element_text(size=17,face="bold"), 
    axis.text.x = element_text(angle = 70, hjust = 1))
  p = p + labs(y = "Relative Abundance")
  p = p + guides(guide_legend(ncol = 1), fill = guide_legend(ncol = 3))
  if (!is.null(facet_wrap)) {
    p <- p + facet_wrap(facet_wrap, nrow = 1) + theme(strip.text = element_text(size=15))
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

# Making stacked bar charts for relative abundance of taxa

# According to Phylum
#ra_Phylum_barplot_location <- simple_plot_bar(ps_Location_ra, x="Location", fill="Phylum")
#ggsave(ra_Phylum_barplot_location, filename = "../Output/Analysis/Relative Abundance of Phylum by Location - Bar Plot.svg", dpi=300, width = 12, height = 10)

#ra_Phylum_barplot_structure <- simple_plot_bar(ps_Structure_ra, x="Structure", fill="Phylum")
#ggsave(ra_Phylum_barplot_structure, filename = "../Output/Analysis/Relative Abundance of Phylum by Structure - Bar Plot.svg", dpi=300, width = 12, height = 10)

ra_Phylum_barplot_location_grouped_by_structure <- simple_plot_bar(ps_Location_Structure_ra, x="Location", fill="Phylum", facet_wrap = "Structure") + guides(fill = guide_legend(ncol = 2)) + theme(axis.text=element_text(size=14), legend.text=element_text(size=rel(1.2)), legend.title=element_text(size =rel(1.5)))
ggsave(ra_Phylum_barplot_location_grouped_by_structure, filename = "../Output/Analysis/Relative Abundance of Phylum by Location Grouped by Structure- Bar Plot.svg", dpi=300, width = 12, height = 10)

# According to Class
#ra_Class_barplot_location <- simple_plot_bar(ps_Location_ra, x="Location", fill="Class")
#ggsave(ra_Class_barplot_location, filename = "../Output/Analysis/Relative Abundance of Class by Location - Bar Plot.svg", dpi=300, width = 12, height = 10)

#ra_Class_barplot_structure <- simple_plot_bar(ps_Structure_ra, x="Structure", fill="Class")
#ggsave(ra_Class_barplot_structure, filename = "../Output/Analysis/Relative Abundance of Class Structure - Bar Plot.svg", dpi=300, width = 12, height = 10)

ra_Class_barplot_location_grouped_by_structure <- simple_plot_bar(ps_Location_Structure_ra, x="Location", fill="Class", facet_wrap = "Structure") + theme(axis.text=element_text(size=13.5), legend.text=element_text(size=rel(1.0)), legend.title=element_text(size =rel(1.5)))
ggsave(ra_Class_barplot_location_grouped_by_structure, filename = "../Output/Analysis/Relative Abundance of Class by Location Grouped by Structure- Bar Plot.svg", dpi=300, width = 12, height = 10)

# According to Order
ra_Order_barplot_location <- simple_plot_bar(ps_Location_ra, x="Location", fill="Order") + theme(legend.text=element_text(size=rel(0.9)), legend.key.size= unit(0.5, "line"), legend.title=element_text(size =rel(1.5))) + guides(fill = guide_legend(nrow = 50))
ggsave(ra_Order_barplot_location, filename = "../Output/Analysis/Relative Abundance of Order by Location - Bar Plot.svg", dpi=300, width = 12, height = 10)

ra_Order_barplot_structure <- simple_plot_bar(ps_Structure_ra, x="Structure", fill="Order") + theme(legend.text=element_text(size=rel(0.9)), legend.key.size= unit(0.5, "line"), legend.title=element_text(size =rel(1.5))) + guides(fill = guide_legend(nrow = 50))
ggsave(ra_Order_barplot_structure, filename = "../Output/Analysis/Relative Abundance of Order by Structure - Bar Plot.svg", dpi=300, width = 12, height = 10)

#ra_Order_barplot_location_grouped_by_structure <- simple_plot_bar(ps_Location_Structure_ra, x="Location", fill="Order", facet_wrap = "Structure") + theme(legend.text=element_text(size=rel(0.5)), legend.key.size= unit(0.2, "line"))
#ggsave(ra_Order_barplot_location_grouped_by_structure, filename = "../Output/Analysis/Relative Abundance of Order by Location Grouped by Structure- Bar Plot.svg", dpi=300, width = 12, height = 10)

# According to Family
ra_Family_barplot_location <- simple_plot_bar(ps_Location_ra, x="Location", fill="Family") + theme(legend.text=element_text(size=rel(0.9)), legend.key.size= unit(0.5, "line"), legend.title=element_text(size =rel(1.5))) + guides(fill = guide_legend(nrow = 50))
ggsave(ra_Family_barplot_location, filename = "../Output/Analysis/Relative Abundance of Family by Location - Bar Plot.svg", dpi=300, width = 12, height = 10)

ra_Family_barplot_structure <- simple_plot_bar(ps_Structure_ra, x="Structure", fill="Family") + theme(legend.text=element_text(size=rel(0.9)), legend.key.size= unit(0.5, "line"), legend.title=element_text(size =rel(1.5))) + guides(fill = guide_legend(nrow = 50))
ggsave(ra_Family_barplot_structure, filename = "../Output/Analysis/Relative Abundance of Family by Structure - Bar Plot.svg", dpi=300, width = 12, height = 10)

#ra_Family_barplot_location_grouped_by_structure <- simple_plot_bar(ps_Location_Structure_ra, x="Location", fill="Family", facet_wrap = "Structure") + theme(legend.text=element_text(size=rel(0.5)), legend.key.size= unit(0.2, "line"))
#ggsave(ra_Family_barplot_location_grouped_by_structure, filename = "../Output/Analysis/Relative Abundance of Family by Location Grouped by Structure- Bar Plot.svg", dpi=300, width = 12, height = 10)

# According to Genus
#ra_Genus_barplot_location <- simple_plot_bar(ps_Location_ra, x="Location", fill="Genus") + theme(legend.text=element_text(size=rel(0.8)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80))
#ggsave(ra_Genus_barplot_location, filename = "../Output/Analysis/Relative Abundance of Genus by Location - Bar Plot.svg", dpi=300, width = 12, height = 10)

#ra_Genus_barplot_structure <- simple_plot_bar(ps_Structure_ra, x="Structure", fill="Genus") + theme(legend.text=element_text(size=rel(0.8)), legend.key.size= unit(0.3, "line")) + guides(fill = guide_legend(nrow = 80))
#ggsave(ra_Genus_barplot_structure, filename = "../Output/Analysis/Relative Abundance of Genus by Structure - Bar Plot.svg", dpi=300, width = 12, height = 10)





### Venn Diagram according to Structure ###
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
#totalSamples_per_Structure <- c(sum(ps1@sam_data$Structure.DNA.Extracted.from=="Leaf"), 
#                                sum(ps1@sam_data$Structure.DNA.Extracted.from=="Rhizome"), 
#                                sum(ps1@sam_data$Structure.DNA.Extracted.from=="Root"), 
#                                sum(ps1@sam_data$Structure.DNA.Extracted.from=="Sediment"))
#for (i in 1:nrow(otutable)){
#  otutable[i,] <- otutable[i,]/totalSamples_per_Structure[i]
#}
otutable <- as.data.frame(otu_table(structure_ps))
otutable[otutable>0] <- 1
filterVector <- c()
for (i in 1:ncol(otutable)){
  x <- sum(otutable[,i])/nrow(otutable) >= 0.95
  x <- sum(otutable[,i])/nrow(otutable)
  filterVector <- c(filterVector, x)
}
patable <- t(otutable)

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
                               lwd = 2, col = col_structure, fill = col_structure, alpha = 0.6, cex = 2, cat.cex = 2.5)

ggsave(venn_diagram, filename = "../Output/Analysis/Venn Diagram according to Structure.svg", dpi = 300)





### Heatmap ###
melted_ps <- psmelt(ps1ra)

# By Phylum
factor_melted_ps_by_sample <- unique(melted_ps$Sample) #factoring by sample
dataframe <- data.frame()  #creating dataframe to fill in information

for (i in factor_melted_ps_by_sample) {                                    #For each sample,
  sub_ps <- melted_ps[melted_ps$Sample == i,]  #subset dataframe.
  
  per_phylum_abundance <- aggregate(Abundance ~ Phylum, data = sub_ps, sum)  #Aggregate abundance based on each phylum
  
  per_phylum_abundance$Sample <- sub_ps$Sample[1:nrow(per_phylum_abundance)] #Adding sample name/structure to each phylum row
  per_phylum_abundance$Structure <- sub_ps$Structure.DNA.Extracted.from[1:nrow(per_phylum_abundance)]
  per_phylum_abundance$Location <- sub_ps$Location[1:nrow(per_phylum_abundance)]
  
  dataframe <- rbind(dataframe, per_phylum_abundance)                       #store this in a dataframe for each row
}

#Sorting dataframe in order of structure, sub-ordered by location
ordered_df <- data.frame()
structures_list <- c("Leaf", "Rhizome", "Root", "Sediment")

for (i in structures_list) {
  struct_frame <- dataframe[dataframe$Structure == i,]
  struct_frame <- struct_frame[order(struct_frame$Location),]
  ordered_df <- rbind(ordered_df, struct_frame)
}

#Plotting
table(ordered_df$Structure) / sum(unique(ordered_df$Phylum) == unique(ordered_df$Phylum))
heatplot <- ggplot(ordered_df, aes(reorder(Sample, -desc(Structure)), reorder(Phylum, desc(Phylum)))) +
  geom_tile(aes(fill = Abundance)) + 
  labs(y = "Phylum", x = "Samples") + 
  theme(axis.text.x = element_blank(), axis.title = element_text(size = 20), axis.text.y = element_text(size = 13), legend.text = element_text(size = 13), legend.title = element_text(size = 15)) + 
  scale_fill_gradient(low = "#FFFFFF", high = "#680000") +
  geom_vline(xintercept = c(58, 117, 177), alpha = 0.2)

#Adding labels
heatplot <- heatplot +
  annotate("rect", xmin = 0, xmax = 58, ymin = 43, ymax = 46,
           alpha = .6, fill = col_structure[1]) +
  annotate("rect", xmin = 58, xmax = 117, ymin = 43, ymax = 46,
           alpha = .6, fill = col_structure[2]) +
  annotate("rect", xmin = 117, xmax = 177, ymin = 43, ymax = 46,
           alpha = .6, fill = col_structure[3]) +
  annotate("rect", xmin = 177, xmax = 238, ymin = 43, ymax = 46,
           alpha = .6, fill = col_structure[4]) + 
  annotate("text", x = 29, y = 44.5, label = "Leaf", size = 10) +
  annotate("text", x = 87.5, y = 44.5, label = "Rhizome", size = 10) +
  annotate("text", x = 147, y = 44.5, label = "Root", size = 10) +
  annotate("text", x = 207, y = 44.5, label = "Sediment", size = 10)
heatplot
ggsave(heatplot, filename = "../Output/Analysis/Heatmap of Phylum Grouped by Structure.svg", dpi = 300, width = 12, height = 10)


# By Class
factor_melted_ps_by_sample <- unique(melted_ps$Sample) #factoring by sample
dataframe <- data.frame()  #creating dataframe to fill in information

for (i in factor_melted_ps_by_sample) {                                    #For each sample,
  sub_ps <- melted_ps[melted_ps$Sample == i,]  #subset dataframe.
  
  per_class_abundance <- aggregate(Abundance ~ Class, data = sub_ps, sum)  #Aggregate abundance based on each class
  
  per_class_abundance$Sample <- sub_ps$Sample[1:nrow(per_class_abundance)] #Adding sample name/structure to each class row
  per_class_abundance$Structure <- sub_ps$Structure.DNA.Extracted.from[1:nrow(per_class_abundance)]
  per_class_abundance$Location <- sub_ps$Location[1:nrow(per_class_abundance)]
  
  dataframe <- rbind(dataframe, per_class_abundance)                       #store this in a dataframe for each row
}

#Sorting dataframe in order of structure, sub-ordered by location
ordered_df <- data.frame()
structures_list <- c("Leaf", "Rhizome", "Root", "Sediment")

for (i in structures_list) {
  struct_frame <- dataframe[dataframe$Structure == i,]
  struct_frame <- struct_frame[order(struct_frame$Location),]
  ordered_df <- rbind(ordered_df, struct_frame)
}

#Plotting
table(ordered_df$Structure) / sum(unique(ordered_df$Class) == unique(ordered_df$Class))
heatplot <- ggplot(ordered_df, aes(reorder(Sample, -desc(Structure)), reorder(Class, desc(Class)))) +
  geom_tile(aes(fill = Abundance)) + 
  labs(y = "Class", x = "Samples") + 
  theme(axis.text.x = element_blank(), axis.title = element_text(size = 20), axis.text.y = element_text(size = 10), legend.text = element_text(size = 9), legend.title = element_text(size = 10)) + 
  scale_fill_gradient(low = "#FFFFFF", high = "#680000") +
  geom_vline(xintercept = c(58, 117, 177), alpha = 0.2)

#Adding labels
heatplot <- heatplot +
  annotate("rect", xmin = 0, xmax = 58, ymin = 79, ymax = 82,
           alpha = .6, fill = col_structure[1]) +
  annotate("rect", xmin = 58, xmax = 117, ymin = 79, ymax = 82,
           alpha = .6, fill = col_structure[2]) +
  annotate("rect", xmin = 117, xmax = 177, ymin = 79, ymax = 82,
           alpha = .6, fill = col_structure[3]) +
  annotate("rect", xmin = 177, xmax = 238, ymin = 79, ymax = 82,
           alpha = .6, fill = col_structure[4]) + 
  annotate("text", x = 29, y = 80.5, label = "Leaf", size = 10) +
  annotate("text", x = 87.5, y = 80.5, label = "Rhizome", size = 10) +
  annotate("text", x = 147, y = 80.5, label = "Root", size = 10) +
  annotate("text", x = 207, y = 80.5, label = "Sediment", size = 10)
heatplot
ggsave(heatplot, filename = "../Output/Analysis/Heatmap of Class Grouped by Structure.svg", dpi = 300, width = 12, height = 10)





### Shannon Diversity Plots ###
div <- data.frame(Location = ps1ra@sam_data$Location,
                  Structure = ps1ra@sam_data$Structure.DNA.Extracted.from,
                  Shannon = diversity(otu_table(ps1ra)))
Richness = colSums(decostand(otu_table(ps1ra), method = "pa"))
write.csv(div, file = "../Output/Analysis/Diversity_table.csv", quote = FALSE)
# By Location
div %>% group_by(Location) %>%
  summarise(N = n(), Mean = mean(Shannon))
ggplot(div, aes(x=Location, y=Shannon, fill = Location)) + 
  geom_boxplot() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 25), axis.title.y = element_text(vjust = 2), axis.title = element_text(size = 30), title = element_text(size = 17), legend.text = element_text(size = 25), legend.title = element_text(size = 30)) + 
  labs(y="Shannon Diversity") +
  #ggtitle("Shannon Diversity by Location") +
  scale_fill_manual(values = col_location)
ggsave(filename = "../Output/Analysis/Shannon_Diversity_by_Location.svg", dpi = 300, width = 12, height = 10)
# By Structure
div %>% group_by(Structure) %>%
  summarise(N = n(), Mean = mean(Shannon))
ggplot(div, aes(x=Structure, y=Shannon, fill = Structure)) + 
  geom_boxplot() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 25), axis.title.y = element_text(vjust = 2), axis.title = element_text(size = 30), title = element_text(size = 17), legend.text = element_text(size = 25), legend.title = element_text(size = 30)) + 
  labs(y="Shannon Diversity") +
  #ggtitle("Shannon Diversity by Structure") +
  scale_fill_manual(values = col_structure)
ggsave(filename = "../Output/Analysis/Shannon_Diversity_by_Structure.svg", dpi = 300, width = 12, height = 10)
# Breakdown Plots
# 1
ggplot(div, aes(x=Location, y=Shannon, fill = Structure)) + 
  geom_boxplot() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + 
  labs(y="Shannon Diversity") +
  ggtitle("Breakdown Plot 1 of Shannon Diversity")
ggsave(filename = "../Output/Analysis/Shannon_Diversity_Breakdown_1.png", dpi = 300)
# 2
ggplot(div, aes(x=Structure, y=Shannon, fill = Location)) + 
  geom_boxplot() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + 
  labs(y="Shannon Diversity") +
  ggtitle("Breakdown Plot 2 of Shannon Diversity")
ggsave(filename = "../Output/Analysis/Shannon_Diversity_Breakdown_2.png", dpi = 300)





### Ordination(s) ###
NMDS = ordinate(ps1ra, method = "NMDS", distance = "bray", trymax = 100)
PCoA = ordinate(ps1ra, method = "PCoA", distance = "bray", trymax = 100)

svg("../Output/Analysis/Full_NMDS_Stress_Plot.svg")
p_stress_full <- stressplot(NMDS, title("Stress Plot for NMDS"))
dev.off()

NMDS2 = data.frame(NMDS1 = NMDS$points[,1], NMDS2 = NMDS$points[,2],group=ps1ra@sam_data$Structure.DNA.Extracted.from, group2=ps1ra@sam_data$Location)
NMDS2.mean=aggregate(NMDS2[,1:2],list(group=ps1ra@sam_data$Structure.DNA.Extracted.from),mean)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame()
for(g in levels(NMDS2$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS2[NMDS2$group==g,],
                                                   veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                ,group=g))
}
p_NMDS1 <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = group2, shape = group), size = 4) +
  ggtitle(paste("NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + theme_bw() + scale_color_manual(values = col_location) + 
  theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
        legend.text = element_text(size = 20)) + 
  guides(shape = guide_legend(override.aes = list(size = 5)), color = guide_legend(override.aes = list(size = 5))) + 
  labs(color = "Location", shape = "Structure")
p_NMDS2 <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = group, shape = group2), size = 4) +
  ggtitle(paste("NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + theme_bw() + scale_color_manual(values = col_structure) + 
  theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
        legend.text = element_text(size = 20)) + 
  guides(shape = guide_legend(override.aes = list(size = 5)), color = guide_legend(override.aes = list(size = 5))) + 
  labs(color = "Structure", shape = "Location")
p_NMDS1_ell <- ggplot(data = NMDS2, aes(NMDS1, NMDS2)) + geom_point(aes(color = group, shape = group2), size = 4) +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group), size=2, linetype=2) +
  ggtitle(paste("NMDS (Stress Value = ", toString(round(NMDS$stress, digits = 3)), ")", sep = "")) + theme_bw() + scale_color_manual(values = col_structure) + 
  theme(axis.title = element_text(size = 20), title = element_text(size = 20), 
        legend.text = element_text(size = 20)) + 
  guides(shape = guide_legend(override.aes = list(size = 5)), color = guide_legend(override.aes = list(size = 5))) + 
  labs(color = "Structure", shape = "Location")
p_NMDS1_ell

ggsave(p_NMDS1, filename = "../Output/Analysis/Full_NMDS_w_Location_colored.svg", dpi = 300, width = 12, height = 10)
ggsave(p_NMDS2, filename = "../Output/Analysis/Full_NMDS_w_Structure_colored.svg", dpi = 300, width = 12, height = 10)
ggsave(p_NMDS1_ell, filename = "../Output/Analysis/Full_NMDS_w_Structure_colored_and_ellipses.svg", dpi = 300, width = 12, height = 10)

# NMDS by Structure
ps1ra_Le <- prune_samples(ps1ra@sam_data$Structure.DNA.Extracted.from == "Leaf", ps1ra)
ps1ra_Rh <- prune_samples(ps1ra@sam_data$Structure.DNA.Extracted.from == "Rhizome", ps1ra)
ps1ra_Ro <- prune_samples(ps1ra@sam_data$Structure.DNA.Extracted.from == "Root", ps1ra)
ps1ra_So <- prune_samples(ps1ra@sam_data$Structure.DNA.Extracted.from == "Sediment", ps1ra)

NMDS_Le <- ordinate(ps1ra_Le, method = "NMDS", distance = "bray", trymax = 100)
NMDS_Rh <- ordinate(ps1ra_Rh, method = "NMDS", distance = "bray", trymax = 100)
NMDS_Ro <- ordinate(ps1ra_Ro, method = "NMDS", distance = "bray", trymax = 100)
NMDS_So <- ordinate(ps1ra_So, method = "NMDS", distance = "bray", trymax = 100)

PCoA_Le <- ordinate(ps1ra_Le, method = "PCoA", distance = "bray", trymax = 100)
PCoA_Rh <- ordinate(ps1ra_Rh, method = "PCoA", distance = "bray", trymax = 100)
PCoA_Ro <- ordinate(ps1ra_Ro, method = "PCoA", distance = "bray", trymax = 100)
PCoA_So <- ordinate(ps1ra_So, method = "PCoA", distance = "bray", trymax = 100)

p_PCoA_Le <- plot_ordination(ps1ra_Le, PCoA_Le, color = "Location", shape = "Structure.DNA.Extracted.from") + ggtitle("(A) Leaf") + 
  theme_bw() + scale_color_manual(values = col_location) + scale_shape_manual(values = 16) + geom_point(size = 3) +
  theme(axis.title = element_text(size = 10), title = element_text(size = 12), 
        legend.text = element_text(size = 10)) + 
  guides(color = guide_legend(override.aes = list(size = 4)), shape = guide_legend(override.aes = list(size = 4, pch = 16))) +
  labs(shape = "Structure")
p_PCoA_Rh <- plot_ordination(ps1ra_Rh, PCoA_Rh, color = "Location", shape = "Structure.DNA.Extracted.from") + ggtitle("(B) Rhizome") + 
  theme_bw() + scale_color_manual(values = col_location) + scale_shape_manual(values = 17) + geom_point(size = 3) +
  theme(axis.title = element_text(size = 10), title = element_text(size = 12), 
        legend.text = element_text(size = 10)) + 
  guides(color = guide_legend(override.aes = list(size = 4)), shape = guide_legend(override.aes = list(size = 4, pch = 17))) +
  labs(shape = "Structure")
p_PCoA_Ro <- plot_ordination(ps1ra_Ro, PCoA_Ro, color = "Location", shape = "Structure.DNA.Extracted.from") + ggtitle("(C) Root") + 
  theme_bw() + scale_color_manual(values = col_location) + scale_shape_manual(values = 15) + geom_point(size = 3) +
  theme(axis.title = element_text(size = 10), title = element_text(size = 12), 
        legend.text = element_text(size = 10)) + 
  guides(color = guide_legend(override.aes = list(size = 4)), shape = guide_legend(override.aes = list(size = 4, pch = 15))) +
  labs(shape = "Structure")
p_PCoA_So <- plot_ordination(ps1ra_So, PCoA_So, color = "Location", shape = "Structure.DNA.Extracted.from") + ggtitle("(D) Sediment") + 
  theme_bw() + scale_color_manual(values = col_location) + scale_shape_manual(values = 3) + geom_point(size = 3) +
  theme(axis.title = element_text(size = 10), title = element_text(size = 12), 
        legend.text = element_text(size = 10)) + 
  guides(color = guide_legend(override.aes = list(size = 4)), shape = guide_legend(override.aes = list(size = 4, pch = 3))) +
  labs(shape = "Structure")

PCoA_by_structure <- ggarrange(p_PCoA_Le, p_PCoA_Rh, p_PCoA_Ro, p_PCoA_So, 
                               ncol = 2, nrow = 2,
                               legend.grob = get_legend(p_NMDS1), legend = "right",
                               common.legend = TRUE)
ggsave(PCoA_by_structure, filename = "../Output/Analysis/PCoA by Structure.svg", dpi = 300, width = 12, height = 10)

NMDS2_Le = data.frame(NMDS1 = NMDS_Le$points[,1], NMDS2 = NMDS_Le$points[,2],group=ps1ra_Le@sam_data$Location)
NMDS2_Le.mean=aggregate(NMDS2_Le[,1:2],list(group=ps1ra_Le@sam_data$Structure.DNA.Extracted.from),mean)
NMDS2_Rh = data.frame(NMDS1 = NMDS_Rh$points[,1], NMDS2 = NMDS_Rh$points[,2],group=ps1ra_Rh@sam_data$Location)
NMDS2_Rh.mean=aggregate(NMDS2_Rh[,1:2],list(group=ps1ra_Rh@sam_data$Structure.DNA.Extracted.from),mean)
NMDS2_Ro = data.frame(NMDS1 = NMDS_Ro$points[,1], NMDS2 = NMDS_Ro$points[,2],group=ps1ra_Ro@sam_data$Location)
NMDS2_Ro.mean=aggregate(NMDS2_Ro[,1:2],list(group=ps1ra_Ro@sam_data$Structure.DNA.Extracted.from),mean)
NMDS2_So = data.frame(NMDS1 = NMDS_So$points[,1], NMDS2 = NMDS_So$points[,2],group=ps1ra_So@sam_data$Location)
NMDS2_So.mean=aggregate(NMDS2_So[,1:2],list(group=ps1ra_So@sam_data$Structure.DNA.Extracted.from),mean)

df_ell_Le <- data.frame()
for(g in levels(NMDS2_Le$group)){
  df_ell_Le <- rbind(df_ell_Le, cbind(as.data.frame(with(NMDS2_Le[NMDS2_Le$group==g,],
                                                   veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                ,group=g))
}
df_ell_Rh <- data.frame()
for(g in levels(NMDS2_Rh$group)){
  df_ell_Rh <- rbind(df_ell_Rh, cbind(as.data.frame(with(NMDS2_Rh[NMDS2_Rh$group==g,],
                                                         veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                      ,group=g))
}
df_ell_Ro <- data.frame()
for(g in levels(NMDS2_Ro$group)){
  df_ell_Ro <- rbind(df_ell_Ro, cbind(as.data.frame(with(NMDS2_Ro[NMDS2_Ro$group==g,],
                                                         veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                      ,group=g))
}
df_ell_So <- data.frame()
for(g in levels(NMDS2_So$group)){
  df_ell_So <- rbind(df_ell_So, cbind(as.data.frame(with(NMDS2_So[NMDS2_So$group==g,],
                                                         veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                      ,group=g))
}

p_NMDS_Le <- ggplot(data = NMDS2_Le, aes(NMDS1, NMDS2)) + geom_point(aes(color = group), shape = 16, size = 4) +
  ggtitle(paste("(A) Leaf (Stress Value = ", toString(round(NMDS_Le$stress, digits = 3)), ")", sep = "")) + theme_bw() + scale_color_manual(values = col_location) +
  theme(axis.title = element_text(size = 15), title = element_text(size = 15), 
        legend.text = element_text(size = 15)) + 
  guides(shape = guide_legend(override.aes = list(size = 4)), color = guide_legend(override.aes = list(size = 2))) + 
  labs(color = "Location")
p_NMDS_ell_Le
p_NMDS_Rh <- ggplot(data = NMDS2_Rh, aes(NMDS1, NMDS2)) + geom_point(aes(color = group), shape = 17, size = 4) +
  ggtitle(paste("(B) Rhizome (Stress Value = ", toString(round(NMDS_Rh$stress, digits = 3)), ")", sep = "")) + theme_bw() + scale_color_manual(values = col_location) +
  theme(axis.title = element_text(size = 15), title = element_text(size = 15), 
        legend.text = element_text(size = 15)) + 
  guides(shape = guide_legend(override.aes = list(size = 4)), color = guide_legend(override.aes = list(size = 2))) + 
  labs(color = "Location")
p_NMDS_Rh
p_NMDS_Ro <- ggplot(data = NMDS2_Ro, aes(NMDS1, NMDS2)) + geom_point(aes(color = group), shape = 15, size = 4) +
  ggtitle(paste("(C) Root (Stress Value = ", toString(round(NMDS_Ro$stress, digits = 3)), ")", sep = "")) + theme_bw() + scale_color_manual(values = col_location) +
  theme(axis.title = element_text(size = 15), title = element_text(size = 15), 
        legend.text = element_text(size = 15)) + 
  guides(shape = guide_legend(override.aes = list(size = 4)), color = guide_legend(override.aes = list(size = 2))) + 
  labs(color = "Location")
p_NMDS_Ro
p_NMDS_So <- ggplot(data = NMDS2_So, aes(NMDS1, NMDS2)) + geom_point(aes(color = group), shape = 3, size = 4) +
  ggtitle(paste("(D) Sediment (Stress Value = ", toString(round(NMDS_So$stress, digits = 3)), ")", sep = "")) + theme_bw() + scale_color_manual(values = col_location) +
  theme(axis.title = element_text(size = 15), title = element_text(size = 15), 
        legend.text = element_text(size = 15)) + 
  guides(shape = guide_legend(override.aes = list(size = 4)), color = guide_legend(override.aes = list(size = 2))) + 
  labs(color = "Location")
p_NMDS_So

NMDS_by_structure <- ggarrange(p_NMDS_Le, p_NMDS_Rh, p_NMDS_Ro, p_NMDS_So, 
                               ncol = 2, nrow = 2,
                               legend.grob = get_legend(p_NMDS1), legend = "right",
                               common.legend = TRUE)
ggsave(NMDS_by_structure, filename = "../Output/Analysis/NMDS by Structure.svg", dpi = 300, width = 12, height = 10)


p_NMDS_ell_Le <- ggplot(data = NMDS2_Le, aes(NMDS1, NMDS2)) + geom_point(aes(color = group), shape = 16, size = 4) +
  geom_path(data=df_ell_Le, aes(x=NMDS1, y=NMDS2,colour=group), size=2, linetype=2) +
  ggtitle(paste("(A) Leaf (Stress Value = ", toString(round(NMDS_Le$stress, digits = 3)), ")", sep = "")) + theme_bw() + scale_color_manual(values = col_location) +
  theme(axis.title = element_text(size = 15), title = element_text(size = 15), 
        legend.text = element_text(size = 15)) + 
  guides(shape = guide_legend(override.aes = list(size = 4)), color = guide_legend(override.aes = list(size = 2))) + 
  labs(color = "Location")
p_NMDS_ell_Le
p_NMDS_ell_Rh <- ggplot(data = NMDS2_Rh, aes(NMDS1, NMDS2)) + geom_point(aes(color = group), shape = 17, size = 4) +
  geom_path(data=df_ell_Rh, aes(x=NMDS1, y=NMDS2,colour=group), size=2, linetype=2) +
  ggtitle(paste("(B) Rhizome (Stress Value = ", toString(round(NMDS_Rh$stress, digits = 3)), ")", sep = "")) + theme_bw() + scale_color_manual(values = col_location) +
  theme(axis.title = element_text(size = 15), title = element_text(size = 15), 
        legend.text = element_text(size = 15)) + 
  guides(shape = guide_legend(override.aes = list(size = 4)), color = guide_legend(override.aes = list(size = 2))) + 
  labs(color = "Location")
p_NMDS_ell_Rh
p_NMDS_ell_Ro <- ggplot(data = NMDS2_Ro, aes(NMDS1, NMDS2)) + geom_point(aes(color = group), shape = 15, size = 4) +
  geom_path(data=df_ell_Ro, aes(x=NMDS1, y=NMDS2,colour=group), size=2, linetype=2) +
  ggtitle(paste("(C) Root (Stress Value = ", toString(round(NMDS_Ro$stress, digits = 3)), ")", sep = "")) + theme_bw() + scale_color_manual(values = col_location) +
  theme(axis.title = element_text(size = 15), title = element_text(size = 15), 
        legend.text = element_text(size = 15)) + 
  guides(shape = guide_legend(override.aes = list(size = 4)), color = guide_legend(override.aes = list(size = 2))) + 
  labs(color = "Location")
p_NMDS_ell_Ro
p_NMDS_ell_So <- ggplot(data = NMDS2_So, aes(NMDS1, NMDS2)) + geom_point(aes(color = group), shape = 3, size = 4) +
  geom_path(data=df_ell_So, aes(x=NMDS1, y=NMDS2,colour=group), size=2, linetype=2) +
  ggtitle(paste("(D) Sediment (Stress Value = ", toString(round(NMDS_So$stress, digits = 3)), ")", sep = "")) + theme_bw() + scale_color_manual(values = col_location) +
  theme(axis.title = element_text(size = 15), title = element_text(size = 15), 
        legend.text = element_text(size = 15)) + 
  guides(shape = guide_legend(override.aes = list(size = 4)), color = guide_legend(override.aes = list(size = 2))) + 
  labs(color = "Location")
p_NMDS_ell_So

NMDS_by_structure_and_ell <- ggarrange(p_NMDS_ell_Le, p_NMDS_ell_Rh, p_NMDS_ell_Ro, p_NMDS_ell_So, 
                               ncol = 2, nrow = 2,
                               legend.grob = get_legend(p_NMDS1), legend = "right",
                               common.legend = TRUE)
ggsave(NMDS_by_structure_and_ell, filename = "../Output/Analysis/NMDS by Structure with Ellipses.svg", dpi = 300, width = 12, height = 10)


# Stress Plots
par(mfrow = c(2,2))
stressplot(NMDS_Le, title("Leaf"), pch = 16, p.col = col_structure[1], l.col = "red", lwd = 3)
stressplot(NMDS_Rh, title("Rhizome"), pch = 16, p.col = col_structure[2], l.col = "red", lwd = 3)
stressplot(NMDS_Ro, title("Root"), pch = 16, p.col = col_structure[3], l.col = "red", lwd = 3)
stressplot(NMDS_So, title("Sediment"), pch = 16, p.col = col_structure[4], l.col = "red", lwd = 3)
par(mfrow = c(1,1))





### PERMANOVA Test ###
otu = as.data.frame(otu_table(ps1ra))
meta = as.data.frame(sample_data(ps1ra))
df = data.frame(SampleID = meta$SampleID, Location = meta$Location, Species = meta$Species, Structure = meta$Structure.DNA.Extracted.from)
# Location
permanova_Location <- adonis(otu ~ Location * Structure, data = df)

sink("../Output/Analysis/adonis_Location_table.txt")
noquote(print("PermANOVA Table:"))
permanova_Location
sink(NULL)

pairwise.adonis <- function(x,factors, sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  co = combn(unique(factors),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(co[1,elem],co[2,elem]),] ~ factors[factors %in% c(co[1,elem],co[2,elem])] , method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}
padonis_Location <- pairwise.adonis(otu,as.character(meta$Location))
sink("../Output/Analysis/adonis_Location_table.txt", append = TRUE)
noquote(print("Pairwise adonis between locations (Bonferroni corrected Pvalues):"))
padonis_Location
sink(NULL)
# Structure
permanova_Structure <- adonis(otu ~ Structure * Location, data = df)

sink("../Output/Analysis/adonis_Structure_table.txt")
noquote(print("PermANOVA Table:"))
permanova_Structure
sink(NULL)

padonis_Structure <- pairwise.adonis(otu,as.character(meta$Structure.DNA.Extracted.from))
sink("../Output/Analysis/adonis_Structure_table.txt", append = TRUE)
noquote(print("Pairwise adonis between structures (Bonferroni corrected Pvalues):"))
padonis_Structure
sink(NULL)

## Checking assumption of homogeneity of multivariate dispersion
# According to Location for each structure
# Leaf
distances_data_leaf <- vegdist(otu[meta$Structure.DNA.Extracted.from == "Leaf",], method="bray")
sink("../Output/Analysis/Assumption_check_of_homogeneity_of_multivariate_dispersion_Location.txt")
noquote(print("Checking assumption of homogeneity of multivariate dispersion according to Location in Leaf:"))
anova(betadisper(distances_data_leaf, meta[meta$Structure.DNA.Extracted.from == "Leaf",]$Location))
sink(NULL)
# Rhizome
distances_data_rhizome <- vegdist(otu[meta$Structure.DNA.Extracted.from == "Rhizome",], method="bray")
sink("../Output/Analysis/Assumption_check_of_homogeneity_of_multivariate_dispersion_Location.txt", append = TRUE)
noquote(print("Checking assumption of homogeneity of multivariate dispersion according to Location in Rhizome:"))
anova(betadisper(distances_data_rhizome, meta[meta$Structure.DNA.Extracted.from == "Rhizome",]$Location))
sink(NULL)
# Root
distances_data_root <- vegdist(otu[meta$Structure.DNA.Extracted.from == "Root",], method="bray")
sink("../Output/Analysis/Assumption_check_of_homogeneity_of_multivariate_dispersion_Location.txt", append = TRUE)
noquote(print("Checking assumption of homogeneity of multivariate dispersion according to Location in Root:"))
anova(betadisper(distances_data_root, meta[meta$Structure.DNA.Extracted.from == "Root",]$Location))
sink(NULL)
# Sediment
distances_data_sediment <- vegdist(otu[meta$Structure.DNA.Extracted.from == "Sediment",], method="bray")
sink("../Output/Analysis/Assumption_check_of_homogeneity_of_multivariate_dispersion_Location.txt", append = TRUE)
noquote(print("Checking assumption of homogeneity of multivariate dispersion according to Location in Sediment:"))
anova(betadisper(distances_data_sediment, meta[meta$Structure.DNA.Extracted.from == "Sediment",]$Location))
sink(NULL)

# According to Structure for each location
# Cyrene
distances_data_CY <- vegdist(otu[meta$Location == "Cyrene",], method="bray")
sink("../Output/Analysis/Assumption_check_of_homogeneity_of_multivariate_dispersion_Structure.txt")
noquote(print("Checking assumption of homogeneity of multivariate dispersion according to Structure in Cyrene:"))
anova(betadisper(distances_data_CY, meta[meta$Location == "Cyrene",]$Structure.DNA.Extracted.from))
sink(NULL)
# Merambong Shoal
distances_data_MShoal <- vegdist(otu[meta$Location == "Merambong Shoal",], method="bray")
sink("../Output/Analysis/Assumption_check_of_homogeneity_of_multivariate_dispersion_Structure.txt", append = TRUE)
noquote(print("Checking assumption of homogeneity of multivariate dispersion according to Structure in Merambong Shoal:"))
anova(betadisper(distances_data_MShoal, meta[meta$Location == "Merambong Shoal",]$Structure.DNA.Extracted.from))
sink(NULL)
# Port Dickson
distances_data_PD <- vegdist(otu[meta$Location == "Port Dickson",], method="bray")
sink("../Output/Analysis/Assumption_check_of_homogeneity_of_multivariate_dispersion_Structure.txt", append = TRUE)
noquote(print("Checking assumption of homogeneity of multivariate dispersion according to Structure in Port Dickson:"))
anova(betadisper(distances_data_PD, meta[meta$Location == "Port Dickson",]$Structure.DNA.Extracted.from))
sink(NULL)
# Semakau
distances_data_Sem <- vegdist(otu[meta$Location == "Semakau",], method="bray")
sink("../Output/Analysis/Assumption_check_of_homogeneity_of_multivariate_dispersion_Structure.txt", append = TRUE)
noquote(print("Checking assumption of homogeneity of multivariate dispersion according to Structure in Semakau:"))
anova(betadisper(distances_data_Sem, meta[meta$Location == "Semakau",]$Structure.DNA.Extracted.from))
sink(NULL)
# Sentosa
distances_data_Sent <- vegdist(otu[meta$Location == "Sentosa",], method="bray")
sink("../Output/Analysis/Assumption_check_of_homogeneity_of_multivariate_dispersion_Structure.txt", append = TRUE)
noquote(print("Checking assumption of homogeneity of multivariate dispersion according to Structure in Sentosa:"))
anova(betadisper(distances_data_Sent, meta[meta$Location == "Sentosa",]$Structure.DNA.Extracted.from))
sink(NULL)
# Perhentian Island
distances_data_Per <- vegdist(otu[meta$Location == "Perhentian Island",], method="bray")
sink("../Output/Analysis/Assumption_check_of_homogeneity_of_multivariate_dispersion_Structure.txt", append = TRUE)
noquote(print("Checking assumption of homogeneity of multivariate dispersion according to Structure in Perhentian Island:"))
anova(betadisper(distances_data_Per, meta[meta$Location == "Perhentian Island",]$Structure.DNA.Extracted.from))
sink(NULL)





### Indicspecies ###

#Preparing data for multipatt
otutable <- as.data.frame(otu_table(ps1ra))
taxtable <- as.data.frame(tax_table(ps1ra))
meta <- as.data.frame(sam_data(ps1ra))

#Checking NAs
sum(is.na(taxtable$Genus))

#Removing NAs
otutable_no_NA <- otutable[,!is.na(taxtable$Genus)]
taxtable_no_NA <- taxtable[!is.na(taxtable$Genus),]
sum(is.na(taxtable_no_NA$Genus)) #verifying

## Comparison between structures and between location
#Converting Structure column to numeric
structures <- as.factor(meta$Structure.DNA.Extracted.from)
structure_levels <- levels(structures)
locations <- as.factor(meta$Location)
location_levels <- levels(locations)

## Changing ASV ids to actual taxa in otu table column names
# generate a vector containing the full taxonomy path for all OTUs
wholetax <- do.call(paste, c(as.data.frame(taxtable_no_NA)
                             ["Genus"], 
                             sep = "__"))  # to distinguish from "_" within tax ranks

# paste wholetax and ASV ids together
tmp <- names(otutable_no_NA)
for(i in 1:length(tmp)){
  names(tmp)[i] = paste(wholetax[i], tmp[i], sep = "__")
}

# overwrite old names
names(otutable_no_NA) <- names(tmp)


## Creating grp vector
for (i in 1:length(structure_levels)) {  #Replacing structure names with numeric representative
  meta$Structure.DNA.Extracted.from[meta$Structure.DNA.Extracted.from == structure_levels[i]] <- i
}
structureID_vector <- as.numeric(meta$Structure.DNA.Extracted.from)
# Leaf = 1, Rhizome = 2, Root = 3, Sediment = 4

for (i in 1:length(location_levels)) {  #Replacing structure names with numeric representative
  meta$Location[meta$Location == location_levels[i]] <- i
}
locationID_vector <- as.numeric(meta$Location)
# Cyrene = 1, Merambong Shoal = 2, Perhential Island = 3, Port Dickson = 4, Semakau = 5, Sentosa = 6

#Converting to matrix
communityData <- as.matrix(otutable_no_NA)

##Processing
indval_structure <- multipatt(communityData, structureID_vector, control = how(nperm = 999))
indval_location <- multipatt(communityData, locationID_vector, control = how(nperm = 999))

##Output
options(max.print=1000000)

summary(indval_structure)
sink("../Output/Analysis/IndVal_Structure.txt")
noquote(print("Output of IndicSpecies grouped according to Structure where Leaf = 1, Rhizome = 2, Root = 3, and Sediment = 4"))
summary(indval_structure)
sink(NULL)

summary(indval_location)
sink("../Output/Analysis/IndVal_Location.txt")
noquote(print("Output of IndicSpecies grouped according to Location where Cyrene = 1, Merambong Shoal = 2, Perhential Island = 3, Port Dickson = 4, Semakau = 5, and Sentosa = 6"))
summary(indval_location)
sink(NULL)

options(max.print=1000)


## Refining list of indicator species and considering species combinations to create prediction models
#otucomb <- combinespecies(communityData, max.order = 2)$XC
#dim(otucomb)

#indvalotucomb <- multipatt(otucomb, structureID_vector, duleg = TRUE,
#                           control = how(nperm=999))
#summary(indvalotucomb, indvalcomp = TRUE)

## Structure wise
## Group 1 (Leaf)
# indicators()
sc <- indicators(X=communityData, cluster=structureID_vector, group = 1,
                 max.order = 3, verbose = TRUE, 
                 At=0.5, Bt=0.2,
                 #nboot.ci = 1000,
)
print(sc, sqrtIVt = 0.6)
summary(sc)
## pruneindicators
sc2=pruneindicators(sc, At=0.8, Bt=0.2, verbose=TRUE)
print(sc2)
sink("../Output/Analysis/IndicSpecies.txt")
noquote("Leaf Indicator Species")
print(sc2)
sink(NULL)
## predict.indicators
p <- predict(sc2, communityData)
pcv <- predict(sc2, cv=TRUE)
data.frame(Group2 = structureID_vector, Prob=p, Prob_CV = pcv)

## Group 2 (Rhizome)
# indicators()
sc <- indicators(X=communityData, cluster=structureID_vector, group = 2,
                 max.order = 3, verbose = TRUE, 
                 At=0.5, Bt=0.2,
                 #nboot.ci = 1000,
                 )
print(sc, sqrtIVt = 0.6)
summary(sc)
## pruneindicators
sc2=pruneindicators(sc, At=0.8, Bt=0.2, verbose=TRUE)
print(sc2)
sink("../Output/Analysis/IndicSpecies.txt", append=TRUE)
noquote("Rhizome Indicator Species")
print(sc2)
sink(NULL)
## predict.indicators
p <- predict(sc2, communityData)
pcv <- predict(sc2, cv=TRUE)
data.frame(Group2 = structureID_vector, Prob=p, Prob_CV = pcv)

## Group 3 (Root)
# indicators()
sc <- indicators(X=communityData, cluster=structureID_vector, group = 3,
                 max.order = 3, verbose = TRUE, 
                 At=0.5, Bt=0.2,
                 #nboot.ci = 1000,
)
print(sc, sqrtIVt = 0.6)
summary(sc)
## pruneindicators
sc2=pruneindicators(sc, At=0.8, Bt=0.2, verbose=TRUE)
print(sc2)
sink("../Output/Analysis/IndicSpecies.txt", append=TRUE)
noquote("Root Indicator Species")
print(sc2)
sink(NULL)
## predict.indicators
p <- predict(sc2, communityData)
pcv <- predict(sc2, cv=TRUE)
data.frame(Group2 = structureID_vector, Prob=p, Prob_CV = pcv)

## Group 4 (Sediment)
# indicators()
sc <- indicators(X=communityData, cluster=structureID_vector, group = 4,
                 max.order = 3, verbose = TRUE, 
                 At=0.5, Bt=0.2,
                 #nboot.ci = 1000,
)
print(sc, sqrtIVt = 0.6)
summary(sc)
## pruneindicators
sc2=pruneindicators(sc, At=0.8, Bt=0.2, verbose=TRUE)
print(sc2)
sink("../Output/Analysis/IndicSpecies.txt", append=TRUE)
noquote("Sediment Indicator Species")
print(sc2)
sink(NULL)
## predict.indicators
p <- predict(sc2, communityData)
pcv <- predict(sc2, cv=TRUE)
data.frame(Group2 = structureID_vector, Prob=p, Prob_CV = pcv)

## Location wise
## Group 1 (Cyrene)
# indicators()
sc <- indicators(X=communityData, cluster=locationID_vector, group = 1,
                 max.order = 3, verbose = TRUE, 
                 At=0.5, Bt=0.2,
                 #nboot.ci = 1000,
)
print(sc, sqrtIVt = 0.6)
summary(sc)
## pruneindicators
sc2=pruneindicators(sc, At=0.8, Bt=0.2, verbose=TRUE)
print(sc2)
sink("../Output/Analysis/IndicSpecies.txt", append=TRUE)
noquote("Cyrene Indicator Species")
print(sc2)
sink(NULL)
## predict.indicators
p <- predict(sc2, communityData)
pcv <- predict(sc2, cv=TRUE)
data.frame(Group2 = locationID_vector, Prob=p, Prob_CV = pcv)



## Alternate Steps for IndicSpecies

## Indicspecies pre-steps----
#Preparing data for multipatt
otutable <- as.data.frame(otu_table(ps1))
taxtable <- as.data.frame(tax_table(ps1))
meta <- as.data.frame(sam_data(ps1))

#Checking NAs
sum(is.na(taxtable$Genus))

#Removing NAs
otutable_no_NA <- otutable[,!is.na(taxtable$Genus)]
taxtable_no_NA <- taxtable[!is.na(taxtable$Genus),]
sum(is.na(taxtable_no_NA$Genus)) #verifying

#Redefining colname in terms of genus
colnames(otutable_no_NA) <- taxtable_no_NA$Genus


#Converting to matrix
communityData <- as.matrix(otutable_no_NA)
#Appending numbers to taxa, so no repeats
for (i in 1: ncol(communityData)) {
  name_list <- colnames(communityData)
  name_list[i] <- paste (name_list[i], i)
  colnames(communityData) <- name_list
}


## Indicspecies structure----
#Converting Structure column to numeric
structures <- as.factor(meta$Structure.DNA.Extracted.from)
structure_levels <- levels(structures)

for (i in 1:length(structure_levels)) {  #Replacing structure names with numeric representative
  meta$Structure.DNA.Extracted.from[meta$Structure.DNA.Extracted.from == structure_levels[i]] <- i
}
structureID_vector <- as.numeric(meta$Structure.DNA.Extracted.from)

#Processing
indval_struct <- multipatt(communityData, structureID_vector, control = how(nperm = 999))

sink("./Outpuut/Indic_Species_Summary_by_structure.txt")
summary(indval_struct)
sink(NULL)


#Presenting results
library(reactable)

#Arranging in a nice table
indic_results <- indval_struct$sign
#remove those with associations for all
indic_results <- indic_results[!(indic_results$index == 15),] #removing all association
indic_results <- indic_results[indic_results$p.value< 0.05,] #removing insignificant

#Table with species names
sp_names <- sapply(strsplit(rownames(indic_results), " "), `[`, 1)
sp_names <- paste(sp_names, "sp.", sep = " ")
#New df for this task
ordered_indic <- data.frame(Species = sp_names, indic_results)
rownames(ordered_indic) <- 1:nrow(ordered_indic)
#Naming columns after structure
colnames(ordered_indic)[2:5] <- structure_levels
#Ordering and subsetting
ordered_indic <- ordered_indic[order(ordered_indic$index),]
ordered_indic <- ordered_indic[,c(1:5, 7, 8)] #remove index col
ordered_indic$stat <- round(ordered_indic$stat, digits = 3)#rounding to 3dp

#Making table with reactable
tick <- function(value) {
  if (value == 0) "\u2718"
  else "\u2713"
}

reactable(ordered_indic,
          columns = list(
            Leaf = colDef(cell = tick),
            Rhizome = colDef(cell = tick),
            Root = colDef(cell = tick),
            Sediment = colDef(cell = tick),
            stat = colDef(name = "Indicator Value"),
            p.value = colDef(name = "P-value")),
          defaultColDef = colDef(align = "center",
                                 headerStyle = list(background = "#f7f7f8")),
          bordered = TRUE,
          defaultPageSize = 50,
          filterable = TRUE
)


## Indicspecies location ----
#Converting Location column to numeric
locations <- as.factor(meta$Location)
location_levels <- levels(locations)

for (i in 1:length(location_levels)) {  #Replacing structure names with numeric representative
  meta$Location[meta$Location == location_levels[i]] <- i
}
locationID_vector <- as.numeric(meta$Location)

#Processing
indval_loc <- multipatt(communityData, locationID_vector, control = how(nperm = 999))

sink("./Outpuut/Indic_Species_Summary_by_location.txt")
summary(indval_loc)
sink(NULL)


#Presenting results
library(reactable)

#Arranging in a nice table
indic_results <- indval_loc$sign
#remove those with associations for all
indic_results <- indic_results[!(indic_results$index == 127),] #removing all association
indic_results <- indic_results[indic_results$p.value< 0.05,] #removing insignificant

#Table with species names
sp_names <- sapply(strsplit(rownames(indic_results), " "), `[`, 1)
sp_names <- paste(sp_names, "sp.", sep = " ")
#New df for this task
ordered_indic <- data.frame(Species = sp_names, indic_results)
rownames(ordered_indic) <- 1:nrow(ordered_indic)
#Naming columns after location
location_levels[c(1,4)] <- c("Chek_Jawa", "Merambong_Shoal") #temporary
colnames(ordered_indic)[2:8] <- location_levels
#Ordering and subsetting
ordered_indic <- ordered_indic[order(ordered_indic$index),]
ordered_indic <- ordered_indic[,c(1:8, 10, 11)]
ordered_indic$stat <- round(ordered_indic$stat, digits = 3)#rounding to 3dp

#Making table with reactable
tick <- function(value) {
  if (value == 0) "\u2718"
  else "\u2713"
}

reactable(ordered_indic,
          columns = list(
            Chek_Jawa = colDef(name = "Chek Jawa", cell = tick),
            Cyrene = colDef(cell = tick),
            Langkawi = colDef(cell = tick),
            Merambong_Shoal = colDef(name = "Merambong Shoal", cell = tick),
            Perhentian = colDef(cell = tick),
            Semakau = colDef(cell = tick),
            Sentosa = colDef(cell = tick),
            stat = colDef(name = "Indicator Value"),
            p.value = colDef(name = "P-value")),
          defaultColDef = colDef(align = "center",
                                 headerStyle = list(background = "#f7f7f8")),
          bordered = TRUE,
          defaultPageSize = 50,
          filterable = TRUE
)





### Mantel Test ###
## Extracting Longitude and Latitude data
meta <- as.data.frame(ps1ra@sam_data)
meta$lon <- sapply(strsplit(meta$GPS.Coordinates, " "), `[`, 2)
meta$lon <- as.double(sapply(strsplit(meta$lon, "E"), `[`, 1), length=10)
meta$lat <- as.double(sapply(strsplit(meta$GPS.Coordinates, "N"), `[`, 1), length=10)

geo_full <- data.frame(meta$lon, meta$lat)
geo_Le <- geo_full[meta$Structure.DNA.Extracted.from == "Leaf",]
geo_Rh <- geo_full[meta$Structure.DNA.Extracted.from == "Rhizome",]
geo_Ro <- geo_full[meta$Structure.DNA.Extracted.from == "Root",]
geo_So <- geo_full[meta$Structure.DNA.Extracted.from == "Sediment",]

## Preparing asv tables
otu_full <- as.data.frame(ps1ra@otu_table)
otu_Le <- otu_full[meta$Structure.DNA.Extracted.from == "Leaf",]
otu_Rh <- otu_full[meta$Structure.DNA.Extracted.from == "Rhizome",]
otu_Ro <- otu_full[meta$Structure.DNA.Extracted.from == "Root",]
otu_So <- otu_full[meta$Structure.DNA.Extracted.from == "Sediment",]

## Making distance matrices
# Adundance data frames - bray curtis dissimilarity
dist.otu_full <- vegdist(otu_full, method = "bray")
dist.otu_Le <- vegdist(otu_Le, method = "bray")
dist.otu_Rh <- vegdist(otu_Rh, method = "bray")
dist.otu_Ro <- vegdist(otu_Ro, method = "bray")
dist.otu_So <- vegdist(otu_So, method = "bray")

# Geographic data frame - haversie distance
d.geo_full <- distm(geo_full, fun = distHaversine)
d.geo_Le <- distm(geo_Le, fun = distHaversine)
d.geo_Rh <- distm(geo_Rh, fun = distHaversine)
d.geo_Ro <- distm(geo_Ro, fun = distHaversine)
d.geo_So <- distm(geo_So, fun = distHaversine)

dist.geo_full <- as.dist(d.geo_full)
dist.geo_Le <- as.dist(d.geo_Le)
dist.geo_Rh <- as.dist(d.geo_Rh)
dist.geo_Ro <- as.dist(d.geo_Ro)
dist.geo_So <- as.dist(d.geo_So)

## Running Mantel Test
# Abundance vs Geographic
abund_geo_full <- mantel(dist.otu_full, dist.geo_full, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_geo_Le <- mantel(dist.otu_Le, dist.geo_Le, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_geo_Rh <- mantel(dist.otu_Rh, dist.geo_Rh, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_geo_Ro <- mantel(dist.otu_Ro, dist.geo_Ro, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_geo_So <- mantel(dist.otu_So, dist.geo_So, method = "spearman", permutations = 9999, na.rm = TRUE)

## Saving Output
sink("../Output/Analysis/Mantel_Test.txt")
noquote("Mantel Test on All Samples")
abund_geo_full
noquote("Mantel Test on Leaf Samples")
abund_geo_Le
noquote("Mantel Test on Rhizome Samples")
abund_geo_Rh
noquote("Mantel Test on Root Samples")
abund_geo_Ro
noquote("Mantel Test on Sediment Samples")
abund_geo_So
sink(NULL)





### Multiple Regression on distance matrices ###
dist_MRM <- MRM(dist.otu_full ~ dist.geo_full,  nperm = 9999)
dist_MRM_Le <- MRM(dist.otu_Le ~ dist.geo_Le,  nperm = 9999)
dist_MRM_Rh <- MRM(dist.otu_Rh ~ dist.geo_Rh,  nperm = 9999)
dist_MRM_Ro <- MRM(dist.otu_Ro ~ dist.geo_Ro,  nperm = 9999)
dist_MRM_So <- MRM(dist.otu_So ~ dist.geo_So,  nperm = 9999)
dist_MRM

sink("../Output/Analysis/MRM_Table.txt")
print("Bray-Curtis distance regressed against spatial distance (Multiple regression on matrices) (All Samples):")
print(dist_MRM)
print("Bray-Curtis distance regressed against spatial distance (Multiple regression on matrices) (Leaf Samples):")
print(dist_MRM_Le)
print("Bray-Curtis distance regressed against spatial distance (Multiple regression on matrices) (Rhizome Samples):")
print(dist_MRM_Rh)
print("Bray-Curtis distance regressed against spatial distance (Multiple regression on matrices) (Root Samples):")
print(dist_MRM_Ro)
print("Bray-Curtis distance regressed against spatial distance (Multiple regression on matrices) (Sediment Samples):")
print(dist_MRM_So)
sink(NULL)





### Network plot ###
ig=make_network(ps1ra, max.dist = .8)
set.seed(13)
np_Location <- plot_network(ig, physeq = ps1ra, color = "Location", shape = "Structure.DNA.Extracted.from",label = NULL,point_size = 4) + 
  theme_bw() + scale_color_manual(values = col_location) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        legend.text = element_text(size = 25), legend.title = element_text(size = 30)) + 
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  labs(color = "Location", shape = "Structure")
ggsave(np_Location, filename="../Output/Analysis/Network_Jaccard_Location.svg", dpi=300, height = 10, width = 12)

set.seed(13)
np_Structure <- plot_network(ig, physeq = ps1ra, color = "Structure.DNA.Extracted.from", shape = "Location",label = NULL,point_size = 4) +
  theme_bw() + scale_color_manual(values = col_structure) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        legend.text = element_text(size = 25), legend.title = element_text(size = 30)) + 
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  labs(color = "Structure", shape = "Location")
ggsave(np_Structure, filename = "../Output/Analysis/Network_Jaccard_Structure.svg", dpi=300, height = 10, width = 12)

'
## Matrix Filtration
# Phyloseq object
OTU = otu_table(t(ps1), taxa_are_rows = TRUE)
TAX = tax_table(ps1)
sampledata = sample_data(ps1)
Phymatrice = phyloseq(OTU,TAX, sampledata)

# Filter OTUs that have 50 sequences in at least 3 samples
subset_3samples = genefilter_sample(Phymatrice, filterfun_sample(function(x) x >= 50), A=3)
Phymatrice2 = prune_taxa(subset_3samples, Phymatrice)
ntaxa(Phymatrice2)
nsamples(Phymatrice2)

# Print the filtered matrix
otus <- otu_table(Phymatrice2)
otus <- as.data.frame(otus)

write.table(otus, file="../Output/Analysis/Phymatrice_M2BiPAT_50seqs3samples_brut.txt", sep="\t")

## Covariance estimation using SPIEC-EASI
# Perform SPIEC-EASI
M2BiPAT.OSD.ctl.gl <- spiec.easi(Phymatrice2, method="glasso", lambda.min.ratio=0.01, nlambda=5, icov.select.params = list(rep.num=5, ncores=1))

# Export SPIEC-EASI correlation matrix and model 
secor <- cov2cor(forceSymmetric(getOptCov(OSD.ctl.gl), ifelse(sum(Matrix::tril(getOptCov(OSD.ctl.gl)))>sum(Matrix::triu(getOptCov(OSD.ctl.gl))), "L", "U")))
refit <- forceSymmetric(getRefit(OSD.ctl.gl), ifelse(sum(Matrix::tril(getRefit(OSD.ctl.gl)))>sum(Matrix::triu(getRefit(OSD.ctl.gl))), "L", "U"))
OSD.ctl.gl_net <- adj2igraph(as.matrix(secor*refit),  vertex.attr=list(name=taxa_names(PhyMatrice)))
E(OSD.ctl.gl_net)$weight %>% hist

graph.density(OSD.ctl.gl_net)
OSD.ctl.gl$est$sparsity[getOptInd(OSD.ctl.gl)]

g <- OSD.ctl.gl_net %>%
  as_tbl_graph()
'





### Core Microbiome ###






### Supervised Learning ###







### corncob ###

## Preparing data
ps1_new <- ps1
taxa_names(ps1_new) <- paste0("ASV", taxa_names(ps1))

## Fitting a Model
ps1_filtered <- ps1_new %>% 
  tax_glom("Phylum")

ps1_filtered
tax_table(ps1_filtered)[1:5,]

## Example
corncob <- bbdml(formula = ASV32 ~ 1,
                 phi.formula = ~ 1,
                 data = ps1_filtered)

## Interpreting a Model
plot(corncob)

plot(corncob, total = TRUE)

plot(corncob, total = TRUE, color = "Structure.DNA.Extracted.from")
plot(corncob, color = "Structure.DNA.Extracted.from")

## Adding covariates
corncob_structure <- bbdml(formula = ASV32 ~ Structure.DNA.Extracted.from,
                           phi.formula = ~ Structure.DNA.Extracted.from,
                           data = ps1_filtered)
corncob_location <- bbdml(formula = ASV10 ~ Location,
                           phi.formula = ~ Location,
                           data = ps1_filtered)

plot(corncob_structure, color = "Structure.DNA.Extracted.from", total = TRUE)
plot(corncob_structure, color = "Structure.DNA.Extracted.from")

plot(corncob_location, color = "Location")

## Model Selection
lrtest(mod_null = corncob, mod = corncob_structure)

## Parameter Interpretation
summary(corncob_structure)
summary(corncob_location)

## Analysis for Multiple Taxa
set.seed(1)
da_analysis <- differentialTest(formula = ~ Location,
                                phi.formula = ~ Location,
                                formula_null = ~ 1,
                                phi.formula_null = ~ Location,
                                test = "Wald", boot = FALSE,
                                data = ps1_filtered,
                                fdr_cutoff = 0.05)
da_analysis
da_analysis$significant_taxa

dv_analysis <- differentialTest(formula = ~ Structure.DNA.Extracted.from,
                                phi.formula = ~ Structure.DNA.Extracted.from,
                                formula_null = ~ Structure.DNA.Extracted.from,
                                phi.formula_null = ~ 1,
                                data = ps1_filtered,
                                test = "LRT", boot = FALSE,
                                fdr_cutoff = 0.05)
dv_analysis$significant_taxa

otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = ps1_filtered)
otu_to_taxonomy(OTU = dv_analysis$significant_taxa, data = ps1_filtered)

da_analysis$p[1:5]
da_analysis$p_fdr[1:5]

plot(da_analysis)
summary(da_analysis)
###
set.seed(1)
da_analysis_loc <- differentialTest(formula = ~ Location,
                                phi.formula = ~ Location,
                                formula_null = ~ 1,
                                phi.formula_null = ~ Location,
                                test = "Wald", boot = FALSE,
                                data = ps1_filtered,
                                fdr_cutoff = 0.05)
da_analysis_loc
da_analysis_loc$significant_taxa

dv_analysis_loc <- differentialTest(formula = ~ Location,
                                phi.formula = ~ Location,
                                formula_null = ~ Location,
                                phi.formula_null = ~ 1,
                                data = ps1_filtered,
                                test = "LRT", boot = FALSE,
                                fdr_cutoff = 0.05)
dv_analysis_loc$significant_taxa

corncob_location <- bbdml(formula = ASV21 ~ Location,
                          phi.formula = ~ Location,
                          data = ps1_filtered)
summary(corncob_location)
