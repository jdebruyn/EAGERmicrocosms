setwd("D:/Project - EAGER Microcosms I/EAGER Microcosm I 16S/R analyses")
getwd()
library("ggpubr")

###Alpha Diversity Stats
#chao, invsimpson, coverage generated from mothur>summary.single
#used filtered OTU table (OTUs with only 1 or 2 reads removed)

##Import and prepare data
#Import and merge - should have 78 obs in each
alpha<-read.csv(file="eager.alpha.csv")
meta<-read.csv(file="eager.metadata.csv")
alpham<-merge.data.frame(alpha, meta, by="sample")

#Split into treatment and time subsets
day11<-subset.data.frame(alpham, alpham$day=="11")
day39<-subset.data.frame(alpham, alpham$day=="39")
treat1<-subset.data.frame(alpham, alpham$treatment=="1")
treat2<-subset.data.frame(alpham, alpham$treatment=="2")
treat3<-subset.data.frame(alpham, alpham$treatment=="3")
treat4<-subset.data.frame(alpham, alpham$treatment=="4")
treatmentnames<-c("1", "2", "3", "4")
treat1234<-subset.data.frame(alpham, alpham$treatment=="1"|alpham$treatment=="2"|alpham$treatment=="3"|alpham$treatment=="4")
temp10<-subset.data.frame(alpham, alpham$temperature=="10")
temp20<-subset.data.frame(alpham, alpham$temperature=="20")
temp30<-subset.data.frame(alpham, alpham$temperature=="30")

#Boxplots of chao and invsimpson
library (ggpubr)
p<-ggboxplot(treat1234, x="temperature", y="chao", fill="treatment", xlab="Temperature (\u00B0C)", ylab="Chao") + facet_wrap("day")
set_palette(p, palette=c("#dc267f", "#648fff", "#785ef0", "#fe6100"))

q<-ggboxplot(treat1234, x="temperature", y="invsimpson", fill="treatment", xlab="Temperature (\u00B0C)", ylab="Inverse Simpson Index") + facet_wrap("day")
set_palette(q, palette=c("#dc267f", "#648fff", "#785ef0", "#fe6100"))

#Stripcharts of chao and invsimpson (Didn't like these as much as boxplots)
ggplot(treat1234, aes(x=temperature, y=chao, color=treatment)) + facet_wrap("day") +
  geom_jitter(pos=position_jitterdodge(jitter.width=0.2, dodge.width=0.8)) +
  labs(x="Temperature", y="Chao") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#ANOVAs
anova<-aov(chao~Treatment*Temperature*Day, data=treat1234)
summary(anova)
TukeyHSD(anova)

anova11c<-aov(chao~Treatment*Temperature, data=day11)
summary(anova11c)
TukeyHSD(anova11c)
anova39c<-aov(chao~Treatment*Temperature, data=day39)
summary(anova39c)
TukeyHSD(anova39c)

anova11i<-aov(invsimpson~Treatment*Temperature, data=day11)
summary(anova11i)
TukeyHSD(anova11i)
anova39i<-aov(invsimpson~Treatment*Temperature, data=day39)
summary(anova39i)
TukeyHSD(anova39i)


#????Means and standard deviations
chao.means<-aggregate(data=alpham, chao~Treatment+Temperature+Day, mean)
chao.sd<-aggregate(data=alpham, chao~Treatment+Temperature+Day, sd)
invsimpson.means<-aggregate(data=alpham, invsimpson~Treatment+Temperature+Day, mean)
invsimpson.sd<-aggregate(data=alpham, invsimpson~Treatment+Temperature+Day, sd)
write.csv(chao.means, file="chao.means.csv")
write.csv(chao.sd, file="chao.sd.csv")
write.csv(invsimpson.means, file="invsimp.means.csv")
write.csv(invsimpson.sd, file="invsimp.sd.csv")


###COMMUNITY COMPOSITION strip charts
#Following Pat Schloss' minimalR tutorial at https://riffomonas.org/minimalR/09_taxonomic_data.html

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
                
library(tidyverse)


#Get taxonomy file ready: remove quote marks, confidence scores, and final semicolon from taxonomy file, then split into columns
taxonomy <- read_tsv(file="eager.opti_mcc.0.03.cons.taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern='[(]\\d*[)]', replacement='')) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern='["]', replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep=";")
write.csv(taxonomy, file="eager.taxonomy.csv")

#Get OTU data ready, get rel abundance (total seqs in each sample=44340)
otu_data <- read_tsv(file="eager.opti_mcc.0.03.subsample.shared", col_types=cols(Group=col_character())) %>%
  select(-label, -numOtus) %>%
  rename(sample=Group) %>%
  pivot_longer(cols=-sample, names_to="otu", values_to="count") %>%
  mutate(rel_abund=count/44340)
write.csv(otu_data, file="eager.otu_data.csv")  

#How to find number of seqs in each sample: 44340
otu_data %>% group_by(sample) %>% summarize (n=sum(count)) %>% summary()

#Join OTU data, taxonomy, and metadata
metadata <- read_csv(file="eager.metadata.csv")

agg_phylum_data <- inner_join(otu_data, taxonomy) %>%
  group_by(sample, phylum) %>%
  summarize(agg_rel_abund=sum(rel_abund)) %>%
  inner_join(., metadata) %>%
  ungroup()

#Get median rel abund to determine top phyla
phyla_median <- agg_phylum_data %>%
  group_by(phylum) %>%
  summarize(median=median(agg_rel_abund)) %>%
  arrange((desc(median)))

top_phyla <- agg_phylum_data %>%
  group_by(phylum) %>%
  summarize(median=median(agg_rel_abund)) %>%
  arrange((desc(median))) %>%
  top_n(13, median) %>%
  pull(phylum) # use pull to convert the names from a data frame to a vector of names

#stripchart
agg_phylum_data %>%
  filter(phylum %in% top_phyla) %>%
  mutate(phylum=factor(phylum, levels=top_phyla)) %>%
  ggplot(aes(x=phylum, y=agg_rel_abund, color=treatment, shape=temperature)) +
  facet_wrap("day")+
  geom_jitter(pos=position_jitterdodge(jitter.width=0.2, dodge.width=0.8)) +
  labs(x="Phylum", y="Relative abundance") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



###COMMUNITY COMPOSITION with PHYLOSEQ
library(phyloseq)
packageVersion("phyloseq")
library(ggplot2)
packageVersion("ggplot2")
library(readxl)
library(ggplot2)#graphing
library(viridis)
library(ggpubr)#graphing
library(vegan)
library(reshape2)
library(phyloseq)
library(lme4)
library(lmerTest)
library(fBasics)
library(knitr)
library(multcomp)
library(multcompView)
library(tidyverse)
library(dplyr)

#import mothur data
mothur_data<-import_mothur(mothur_shared_file = "eager.opti_mcc.0.03.subsample.shared", mothur_constaxonomy_file = "eager.opti_mcc.0.03.cons.taxonomy")
treatments<-read.csv(file="eager.metadata.csv")
treatments<-sample_data(treatments)
rownames(treatments)<-treatments$sample

#merge files into a phyloseq object "moth_merge":
map=treatments #make a copy of the metadata; this copy will be reformatted for phyloseq
map=sample_data(map) #convert metadata into phyloseq format
rownames(map)=map$sample #assign rownames to be group
moth_merge = merge_phyloseq(mothur_data, map) #merge mothurdata with the metadata
colnames(tax_table(moth_merge)) #print column names of tax file
colnames(tax_table(moth_merge)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus") #rename tax column names

#subset by day
moth1=subset_samples(moth_merge, day == "1") 
moth11=subset_samples(moth_merge, day == "11") 
moth39=subset_samples(moth_merge, day == "39")
moth1139=merge_phyloseq(moth11, moth39)


#Community composition: Stacked bar charts
#separate initial samples (CS) from the microcosms (day 11 and 39)
moth_abund = moth_merge %>% transform_sample_counts(function(x) {x/sum(x)}) #take phyloseq file with read counts and make a new phyloseq object containing relative abundances on each OTU

#phylum1
moth1_phylum = moth1 %>%
  tax_glom(taxrank="Phylum") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() %>%
  filter(Abundance > 0.02) %>%
  arrange(Phylum) #take the master phyloseq object with read counts, group all OTUs by phylum, make a new phyloseq object containing relative abundances of each phylum, melt into a dataframe, remove any phyla with rel. abundance less than 2%, and arrange by phylum

moth1139_phylum = moth1139 %>%
  tax_glom(taxrank="Phylum") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() %>%
  filter(Abundance > 0.02) %>%
  arrange(Phylum) #take the master phyloseq object with read counts, group all OTUs by phylum, make a new phyloseq object containing relative abundances of each phylum, melt into a dataframe, remove any phyla with rel. abundance less than 2%, and arrange by phylum


#class
moth1_class = moth1 %>%
  tax_glom(taxrank="Class") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() %>%
  filter(Abundance > 0.02) %>%
  arrange(Class)

moth1139_class = moth1139 %>%
  tax_glom(taxrank="Class") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() %>%
  filter(Abundance > 0.02) %>%
  arrange(Class)

phylum_colors2 <- c(
  "#332288", "#88CCEE", "#44AA99","#117733", "#999933", "#DDCC77", 
  "#CC6677", "#882255", "#AA4499", "#4477AA", "#BBBBBB", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "firebrick2", "#599861")

#plot relative abundance, facetting by origin 
p1= ggplot(moth1_phylum, aes(x=treatment, y=Abundance, fill = Phylum)) + 
  facet_grid(~temperature, scales="free") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors2) + 
  theme(axis.text.x = element_text(size = 12), panel.grid.major=element_blank(), panel.background = element_blank()) +  
  labs(y="Relative Abundance")+
  theme(legend.position="bottom") +
  guides(fill=guide_legend(ncol=3))

p1139= ggplot(moth1139_phylum, aes(x=treatment, y=Abundance, fill = Phylum)) + 
  facet_grid(day~temperature, scales="free") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors2) + 
  theme(axis.text.x = element_text(size = 12), panel.grid.major=element_blank(), panel.background = element_blank()) +  
  labs(y="Relative Abundance")+
  theme(legend.position="bottom") +
  guides(fill=guide_legend(ncol=3))

c1= ggplot(moth_class, aes(x=treatment, y=Abundance, fill = Class)) + 
  facet_grid(day~temperature, scales="free") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors2) + 
  theme(axis.text.x = element_text(size = 12), panel.grid.major=element_blank(), panel.background = element_blank()) +  
  labs(y="Relative Abundance")+
  theme(legend.position="bottom") +
  guides(fill=guide_legend(ncol=3))

c1139= ggplot(moth1139_class, aes(x=treatment, y=Abundance, fill = Class)) + 
  facet_grid(day~temperature, scales="free") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors2) + 
  theme(axis.text.x = element_text(size = 12), panel.grid.major=element_blank(), panel.background = element_blank()) +  
  labs(y="Relative Abundance")+
  theme(legend.position="bottom") +
  guides(fill=guide_legend(ncol=3))



#Heatmaps to look at individuatl OTUs
firm<-subset_taxa(moth_merge, Phylum=="Firmicutes")
plot_heatmap(firm)


clos<-subset_taxa(moth_merge, Class=="Clostridia")
plot_heatmap(clos, sample.label="DayTreatTemp")

bac<-subset_taxa(moth_merge, Class=="Bacilli")
plot_heatmap(bac, sample.label="treatment")





#Bar charts of individual taxa
top_phyla<-moth1139_phylum%>%
  group_by(Phylum) %>%
  summarize(median=median(Abundance)) %>%
  arrange((desc(median))) %>% # keep this so that the phyla are sorted properly
  top_n(5, median) %>%
  pull(Phylum) # use pull to convert the names from a data frame to a vector of names

moth1139_phylum %>%
  filter(Phylum %in% top_phyla) %>%
  ggplot(aes(x=Phylum, y=Abundance, color=treatment)) +
  geom_boxplot()+
  facet_grid(day~temperature)+
  theme_bw()

#Plot firmicutes
moth1139_class %>%
  filter(Class == "Bacilli") %>%
ggplot(aes(x=treatment, y=Abundance, color=treatment, shape=temperature)) +
  geom_point()+
  facet_grid(day~temperature)+
  theme_bw()



###Network analysis
plot_net(moth_merge, maxdist=0.2, color = "treatment", shape = "temperature", point_label = "day")


##SOURCETRACKER
#Documentaion and example: https://github.com/danknights/sourcetracker/blob/master/example.r
#Another example: https://mgaley-004.github.io/MiSeq-Analysis/Tutorials/SourceSink.html

#Load metadata and sample data
meta<-read.table(file="eager.metadata.txt", sep='\t', header=T, row.names=1, comment='')
otu<-read.table(file="eager.shared2.FEAST.txt", sep='\t', header=T, row.names=1, check=F, comment='')
otus<-as.matrix(otu)

#Extract source environments and source/sink indices
train.ix<-which(meta$SourceSink=='source')
test.ix<-which(meta$SourceSink=='sink')
envs<-meta$Env

source('sourcetracker-1.0.1/src/SourceTracker.r')
#tune.results<-tune.st(otu[train.ix,], envs[train.ix]) ##optional
alpha1<-alpha2<-0.001
st<-sourcetracker(otu[train.ix], envs[train.ix])
###ERROr: Error in if (!is.element(class(x), c("matrix", "data.frame", "array"))) x <- matrix(x,  : 
#the condition has length > 1



#Next steps
results<-predict(st, otu[test.ix], alpha1=alpha1, alpha2=alpha2)






##FEAST: Fast Expectation maximization microbial Source Tracking
#Paper: Shenhav et al., Nature Methods 2019 (https://www.nature.com/articles/s41592-019-0431-x).
#Documentaion and example: https://github.com/cozygene/FEAST

#Load packages

Packages <- c("Rcpp", "RcppArmadillo", "vegan", "dplyr", "reshape2", "gridExtra", "ggplot2", "ggthemes")
install.packages(Packages)
lapply(Packages, library, character.only = TRUE)

devtools::install_github("cozygene/FEAST", force = T)



#OTU table
#load mothur shared file, then remove columns of "label" and "numOTUs"
#documentation is confusing about otu table orientation (samples as rows vs. columns), so tried it both ways! 
#Developer confirmed it should be samples as rows, so transform is needed

#format otu table for FEAST
otu<-read.table(file="eager.opti_mcc.0.03.subsample.shared", header=T, row.names=2)
otu<-subset(otu, select=-c(label, numOtus))
write.table(otu, file="eager.otu.txt", row.names=TRUE, col.names=TRUE, sep="\t", quote=F) 
#transpose
otut<-t(as.matrix(otu))
write.table(otut, file="eager.otu.t.FEAST.txt", row.names=TRUE, col.names=TRUE, sep="\t", quote=F) 


#load otu table in FEAST format
otutable2<-Load_CountMatrix(CountMatrix_path = "eager.shared2.FEAST.txt")
otutable2t<-Load_CountMatrix(CountMatrix_path = "eager.shared2t.FEAST.txt") #manually transposed so columns are samples

#load metadata 
metadata<-Load_metadata(metadata_path="eager.metadata.FEAST.txt")

#Run FEAST
FEAST_Output<-FEAST(C=otutable2t, metadata=metadata, different_sources_flag=1, outfile="out1")

### Kept getting the same error. "Error in sprintf("Error: there are %d sample ids in common ") :too few arguments""



###Code from developer: 
metadata <- Load_metadata(metadata_path = "eager.metadata.FEAST.txt")
otus <- Load_CountMatrix(CountMatrix_path = "eager.otu.txt")


#Subset of sinks for a fast run
#metadata_sink <- metadata[metadata$SourceSink == "Sink",][c(1:3),]
#metadata_subsource <- metadata[metadata$SourceSink == "Source",]

#metadata_sub <- rbind(metadata_sink, metadata_subsource)
#metadata_sub$SourceSink <- as.character(metadata_sub$SourceSink)



FEAST_output <- FEAST(C = t(otus), metadata = metadata, 
                      different_sources_flag = 0, dir_path = my_path, outfile="demo")

FEAST_output


##Alex Emmons ran the FEAST code
##Output (means) is in FEASTout.csv
##Plot pie charts of FEAST output using ggplot2

library(ggplot2)
feast<-read.csv(file="FEASTout2.csv", header=T)
ggplot(feast, aes(x=Day, y=Proportion, fill=Source))+
  geom_col()  + coord_polar("y") +
  scale_fill_manual(values=c("#AA4499", "#4477AA", "#bbbbbb"))+
  facet_grid(Temperature~Treatment, labeller=label_both) +
  theme_void()









###NEXT STEPS
#Differential abundance?
#Bar plots of individual taxa????

