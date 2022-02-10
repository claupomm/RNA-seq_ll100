#####################################
## R script for reading in salmon quant.sf
## and analysis via DESeq2
## annotation via ensembl meta-data
#####################################



#####################################
## Load Libraries
#####################################
library(DESeq2)
library(gplots)
library(tximportData)
library(GenomicFeatures)
library(tximport)




#####################################
## READ META-DATA
#####################################
pdata = read.table('file2sample.csv',header=T,sep=',')
pdata = pdata
rownames(pdata) = pdata[,1]
dirElem  = strsplit(getwd(),'/')[[1]];
expname  = paste("_salmon_shiny_",dirElem[length(dirElem)],sep='')

files = dir(path='counts/', pattern=".*\\.sf$", recursive=TRUE)
o     = pmatch(paste(pdata$sampleName,".salmon/",sep=""), files)
files = files[o];
files = paste('counts',files,sep='/')
pdata$fileName = files
pdata$group = factor(pdata$group)
pdata$project = factor(pdata$project)




#####################################
## summarise transcripts to genes
#####################################
TxDb <- makeTxDbFromGFF(file = "/path/to/gtf_file/gencode/37/gencode.v37.annotation.gtf") # 
k <- keys(TxDb, keytype = "TXNAME")
tx2gene <- select(TxDb, k, "GENEID", "TXNAME")
head(tx2gene)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
colnames(txi$counts) = pdata$sampleName
colnames(txi$abundance) = pdata$sampleName




#####################################
## NORMALIZE DATA, DESeq2
## negative binomial GLM test
## alternative test: nbionomLRT (chi-square)
#####################################
dds <- DESeqDataSetFromTximport(
      txi,
	  colData = pdata,
	  design= (~ project + group) # experimental design, multiple factors possible
)
dds <- DESeq(dds) # includes estimation of size factors, dispersion + nbinomWaldTest; betaPrior=F for >2 factors

conds = labels(terms(design(dds)))[1]
ens.str <- substr(rownames(dds), 1, 15)
rownames(dds) = ens.str
resultsNames(dds)
nexprs = counts(dds,normalized=T)




#####################################
## annotation via download from website ensembl biomart
#####################################
anno = read.csv(gzfile("/path/to/ensembl_meta_data/biomart_ens103_210308.txt.gz"), header=T, as.is=T, sep="\t")
colnames(anno) = c("ensembl_gene_id", "description", "chromosome_name", "gene_start_position", "gene_end_position", "strand", "external_gene_name", "entrez_gene_id")
anno = anno[order(anno$external_gene_name),]
anno = anno[grep("^CHR", anno$chromosome_name,invert=T),]
anno = anno[!duplicated(anno$ensembl_gene_id),] # one annotation for one gene
id_type= "ensembl_gene_id"
anno = anno[anno$ensembl_gene_id%in%rownames(nexprs), ]

anexprs = data.frame(nexprs)
colnames(anexprs) = colnames(nexprs)
anexprs = merge(anno,anexprs,by.x='ensembl_gene_id',by.y=0,all.y=T)
tpm <- txi$abundance
ens <- substr(rownames(tpm), 1, 15)
rownames(tpm) = ens
atpm = merge(anno,data.frame(tpm),by.x='ensembl_gene_id',by.y=0,all.y=T)

write.table(anexprs,paste('tables/NormData',expname, '.csv',sep=''), row.names=F,quote=F, sep='\t',na="")
write.table(atpm, file=,paste('tables/TPM',expname, '.csv',sep=''), row.names=F,quote=F, sep='\t',na="")


# save R objects for later/further analysis
save(pdata,txi,dds,nexprs,anexprs,atpm,expname,tx2gene,file="R_salmon.rda")



sessionInfo()

