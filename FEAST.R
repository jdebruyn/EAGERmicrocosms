library(tidyverse)
library(phyloseq)

meta<-read.delim("eager.metadata.txt") %>% column_to_rownames("sample")
otu<-read.delim("eager.otu.txt") %>% column_to_rownames("SampleID")
otu1<-as.matrix(otu)


#make a phyloseq object
eager<-phyloseq(sample_data(meta),otu_table(otu1,taxa_are_rows = FALSE))

#look at read counts
samplesums<-data.frame(sample_id=sample_names(eager),sums=sample_sums(eager))
min(samplesums$sums)
median(samplesums$sums)
max(samplesums$sums)


meta<-meta %>% rownames_to_column("#SampleID")
#write_delim(meta,"metadata2.txt",delim="\t")


#reformat otu table
otu2<-data.frame(t(otu)) %>% rownames_to_column("#OTU ID")
#write.table(otu2,"otutable_t.txt",quote=FALSE,sep='\t',row.names=FALSE)


#give FEAST a try
#install
devtools::install_github("cozygene/FEAST")
library(FEAST)

#formatting otu table and metadata
fotu<-t(otu)
fmeta<-meta %>% select(SampleID=`#SampleID`,Env,SourceSink)
fmeta<-fmeta %>% arrange(SourceSink) 
fmeta<-fmeta[match(colnames(fotu), fmeta$SampleID),]
fmeta$SourceSink[fmeta$SourceSink=="source"]<-"Source"
fmeta$SourceSink[fmeta$SourceSink=="sink"]<-"Sink"
rownames(fmeta)<-NULL
#apparently sampleid had to be the rownames
fmeta<-fmeta %>% column_to_rownames("SampleID")

fmeta$id<-c(1:78)
fmeta$id[fmeta$SourceSink=="Source"]<-NA

#run v 0.1.0
FEAST_output <- FEAST(C = otu1, metadata = fmeta, different_sources_flag =0, dir_path = ".",
                      outfile="EAGER")


