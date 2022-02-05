if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pd.clariom.d.human")

install.packages("pd.clariom.d.human")

BiocManager::install("clariomdhumantranscriptcluster.db")

BiocManager::install("oligo")
BiocManager::install("affycoretools")
BiocManager::install("clariomdhumanprobeset.db")

BiocManager::install("limma")

BiocManager::install("DBI")

library("pd.clariom.d.human")
library("oligo")
library("affycoretools")
library("clariomdhumanprobeset.db")
library("clariomdhumantranscriptcluster.db")
library("limma")
library("DBI")


Today<-c("M","l","S")
for (word in Today) {
  print(word)
}

###To know what is in the database

KeyTypes<-keytypes(clariomdhumanprobeset.db)

index<-1
for (KeyType in KeyTypes) {
  print(index)
  print(KeyType)
  index<-index+1
}

Columns<-keytypes(clariomdhumantranscriptcluster.db)
for (Column in Columns) {
  print(Column)
  print(select(clariomdhumantranscriptcluster.db,keys = keys(clariomdhumantranscriptcluster.db)[30000],column=Column)[1,2])
  
}

###Working with CEL-files
uArrayRaw <- read.celfiles(list.celfiles())

pbmc.expr=exprs(uArrayRaw)

library(hexbin)

plot(hexplom(log2(pbmc.expr)))



transcript.eset=rma(uArrayRaw,target="core") #default level of summarization
#probeset.eset=rma(uArrayRaw,target="probeset")


# annotating using Affymetrix-provided information (NetAffx)
#transcript.eset.affy <- annotateEset(transcript.eset, pd.clariom.d.human)

# annotating using info assembled by the Bioconductor Core Team
transcript.eset <- annotateEset(transcript.eset,
                                clariomdhumantranscriptcluster.db)
#probeset.eset <- annotateEset(probeset.eset, clariomdhumanprobeset.db)

#Filter?
#transcript.eset <- getMainProbes(transcript.eset)

#fData(transcript.eset.affy)[c(50792,50803,50980),] #pd.clariom.d.human
#fData(transcript.eset)[c(50792,50803,50980),] #Based on
clariomdhumantranscriptcluster.db

#Work with Limma
#create design matrix
DesignMatrix<-matrix(c(0,1,1,0,0,1,1,0),nrow=4, ncol=2,byrow = TRUE)

rownames(DesignMatrix)<-sampleNames(transcript.eset)
colnames(DesignMatrix)<-c("Treatment","Control")
DesignMatrixPair<-matrix(c(0,1,1,1,0,1,0,1,0,1,0,0),nrow=4, ncol=3,byrow = TRUE)
rownames(DesignMatrixPair)<-sampleNames(transcript.eset)
colnames(DesignMatrixPair)<-c("Treatment","Control","Subject")


#create contrast matrix
FittedModel<-lmFit(transcript.eset, DesignMatrixPair)
Contrasts<-"Trt=Treatment-Control"
ContrastsMatrix<-makeContrasts(Contrasts,levels = DesignMatrixPair)

FitContrast=contrasts.fit(FittedModel,ContrastsMatrix)
efitContrast=eBayes(FitContrast)

##
##Plot the p-values.
par(mfrow=c(1,2))##to generate one graph
for (i in 1:ncol(efitContrast$p.value)){
  hist(efitContrast$p.value[,i],main=colnames(efitContrast$p.value)[i])
}



## Why probeset or gene number number is 100000??

#Pick out results
ResultsTablePair<-topTable(efitContrast,p.value = 0.05, number =
                             100000,adjust.method = "BH")

select(clariomdhumantranscriptcluster.db,keys = ResultsTable[,"PROBEID"],columns = "GENETYPE")


##Compare control, treatment and recovery

###Working with CEL-files for PBMC Shirin
uArrayRaw2 <- read.celfiles(list.celfiles())

#pbmc.expr=exprs(uArrayRaw2)

#write.csv(pbmc.expr,file = "pbmc.expr.csv")

library(hexbin)

plot(hexplom(log2(pbmc.expr)))

# to get the log2 
pbmc.log2<-log2(pbmc.expr)

write.csv(pbmc.log2,file = "pbmc.log2.H.clariomD.csv")


transcript.eset1=rma(uArrayRaw2,target="core") #default level of summarization
#probeset.eset=rma(uArrayRaw,target="probeset")

write.csv(transcript.eset1,file = "transcript.eset1.csv")
# annotating using Affymetrix-provided information (NetAffx)
#transcript.eset.affy <- annotateEset(transcript.eset, pd.clariom.d.human)

##Get expression data from the normalized values
pbmc.expr.eset1=exprs(transcript.eset1)

## transform pbmc.expr.eset1 into log2
#log2.pbmc.expr.eset1<-log2(pbmc.expr.eset1)

#write.csv(log2.pbmc.expr.eset1,file = "log2.pbmc.expr.eset1.csv")

#to get the normalized expression table
write.csv(pbmc.expr.eset1,file = "pbmc.expr.eset1.csv")



# annotating using info assembled by the Bioconductor Core Team
transcript.eset1 <- annotateEset(transcript.eset1,
                                 clariomdhumantranscriptcluster.db)
#probeset.eset <- annotateEset(probeset.eset, clariomdhumanprobeset.db)

#Filter?
#transcript.eset <- getMainProbes(transcript.eset)

#fData(transcript.eset.affy)[c(50792,50803,50980),] #pd.clariom.d.human
#fData(transcript.eset)[c(50792,50803,50980),] #Based on
clariomdhumantranscriptcluster.db

#Work with Limma
#create design matrix
DesignMatrix1<-matrix(c(0,1,0,1,0,0,0,0,1,0,1,0,1,0,0,0,0,1),nrow=6, ncol=3,byrow = TRUE)

colnames(DesignMatrix1)<-c("Treatment","Control","Recovery")
rownames(DesignMatrix1)<-sampleNames(transcript.eset1)



#create contrast matrix
FittedModel<-lmFit(transcript.eset1, DesignMatrix1)

##All contrast
ContrastsMatrix.all<-makeContrasts(Treatment-Control,Treatment-Recovery,Recovery-Control,levels = DesignMatrix1)
FitContrast.all=contrasts.fit(FittedModel,ContrastsMatrix.all)
efitContrast.all=eBayes(FitContrast.all)

ResultsTable.all.con<-topTable(efitContrast.all, number =
                                 100000,adjust.method = "BH",p.value=0.05)


##Pair contrast1
ContrastsMatrix1<-makeContrasts(Treatment-Control,levels = DesignMatrix1)
FitContrast1=contrasts.fit(FittedModel,ContrastsMatrix1)
efitContrast1=eBayes(FitContrast1)

#Pick out results
ResultsTable1<-topTable(efitContrast1, number =
                          100000,adjust.method = "BH",p.value=0.05)

ResultsTable1.fdr<-topTable(efitContrast1, number =
                              100000,adjust.method = "fdr",p.value=0.05)


## Pair contras Recovery-treatment
ContrastsMatrix.rec.trt<-makeContrasts(Recovery-Treatment,levels = DesignMatrix1)
FitContrast.rec.trt=contrasts.fit(FittedModel,ContrastsMatrix.rec.trt)
efitContrast.rec.trt=eBayes(FitContrast.rec.trt)

ResultsTable.rec.trt<-topTable(efitContrast.rec.trt, number =
                                 100000,adjust.method = "BH",p.value=0.05)

write.csv(ResultsTable.rec.trt,file = "ResultsTable.rec.trt.csv")


##Pair contrast2
ContrastsMatrix2<-makeContrasts(Treatment-Recovery,levels = DesignMatrix1)
FitContrast2=contrasts.fit(FittedModel,ContrastsMatrix2)
efitContrast2=eBayes(FitContrast2)

ResultsTable2<-topTable(efitContrast2, number =
                          100000,adjust.method = "BH",p.value=0.05)

##Pair contrast3
ContrastsMatrix3<-makeContrasts(Recovery-Control,levels = DesignMatrix1)
FitContrast3=contrasts.fit(FittedModel,ContrastsMatrix3)
efitContrast3=eBayes(FitContrast3)

ResultsTable3<-topTable(efitContrast3, number =
                          100000,adjust.method = "BH",p.value=0.05)

##Mixed contrast matrix
mMxContrastsMatrix<-makeContrasts(contrasts="Treatment-(Control+Recovery)/2",levels=DesignMatrix1)

FitContrast.mix=contrasts.fit(FittedModel,mMxContrastsMatrix,)
efitContrast.mix=eBayes(FitContrast.mix)

ResultsTable.mix.cons<-topTable(efitContrast.mix, number =
                                  100000,adjust.method = "BH",p.value=0.05)

write.csv(ResultsTable.mix.cons,file = "ResultsTable.mix.cons.csv")



# Fold-change thresholding
#FitContrast.mix2 <- treat(FitContrast.mix,lfc=0.1)
#topTreat(FitContrast.mix2)

#efitContrast.mix.lfc=eBayes(FitContrast.mix2)
#ResultsTable.lfc3<-topTable(efitContrast.mix.lfc, number =
#100000,adjust.method = "BH",p.value=0.05)
#volcanoplot(FitContrast.mix2)

#volcanoplot(efitContrast.mix)



##To get all genes expression and fold change

ResultsTable.mix.cons.all<-topTable(efitContrast.mix, number =
                                      100000,adjust.method = "BH")

write.csv(ResultsTable.mix.cons.all,file = "ResultsTable.mix.pbmc.cons.all.csv")


##Select all ncRNA

RNA.genetype<-select(clariomdhumantranscriptcluster.db,keys = ResultsTable.mix.cons.all[,"PROBEID"],columns = c("GENETYPE","ENTREZID","GENENAME"))

write.csv(RNA.genetype,file = "RNA.genetype.csv")

select.pbmc..ncRNA<-subset(RNA.genetype, GENETYPE == "ncRNA")

write.csv(select.pbmc..ncRNA,file = "select.pbmc.ncRNA.csv")

select.pbmc.lncRNA<-subset(RNA.genetype, GENETYPE == "lncRNA")



# Other pairs

C1.C64.ContrastsMatrix<-makeContrasts(Control-Recovery,levels=DesignMatrix1)

C1.C64.FitContrast.mix=contrasts.fit(FittedModel,C1.C64.ContrastsMatrix)
C1.C64.efitContrast.mix=eBayes(C1.C64.FitContrast.mix)

C1.C64.ResultsTable.mix.cons<-topTable(C1.C64.efitContrast.mix, number =
                                         100000,adjust.method = "BH",p.value=0.05)


write.csv(C1.C64.ResultsTable.mix.cons,file = "C1.C64.ResultsTable.mix.cons.csv")

##To get all genes expression and fold change

C1.C64.ResultsTable.mix.cons.all<-topTable(C1.C64.efitContrast.mix, number =
                                             100000,adjust.method = "BH")

write.csv(C1.C64.ResultsTable.mix.cons.all,file = "C1.C64.ResultsTable.mix.cons.all.csv")


## Select ncRNA

C1.C64.RNA.genetype<-select(clariomdhumantranscriptcluster.db,keys = C1.C64.ResultsTable.mix.cons.all[,"PROBEID"],columns = c("GENETYPE","ENTREZID","GENENAME"))


select.pbmc.C1.C64.ncRNA<-subset(C1.C64.RNA.genetype, GENETYPE == "ncRNA")

write.csv(select.pbmc..ncRNA,file = "select.pbmc.C1.C64.ncRNA.csv")




## Treatment vs recovery


C64.TD1.ContrastsMatrix<-makeContrasts(Recovery-Treatment,levels=DesignMatrix1)

C64.TD1.FitContrast.mix=contrasts.fit(FittedModel,C64.TD1.ContrastsMatrix)
C64.TD1.efitContrast.mix=eBayes(C64.TD1.FitContrast.mix)

C64.TD1.ResultsTable.mix.cons<-topTable(C64.TD1.efitContrast.mix, number =
                                          100000,adjust.method = "BH",p.value=0.05)


write.csv(C64.TD1.ResultsTable.mix.cons,file = "C64.TD1.ResultsTable.mix.cons.csv")

##To get all genes expression and fold change

C64.TD1.ResultsTable.mix.cons.all<-topTable(C64.TD1.efitContrast.mix, number =
                                              100000,adjust.method = "BH")

KeyTypes<-keytypes(clariomdhumanprobeset.db)

select(clariomdhumantranscriptcluster.db,keys = ResultsTable.mix.cons[,"PROBEID"],columns = "GENETYPE")


ResultsTable.mix.cons.genetype<-select(clariomdhumantranscriptcluster.db,keys = ResultsTable.mix.cons[,"PROBEID"],columns = c("GENETYPE","ENTREZID","GENENAME"))


write.csv(ResultsTable.mix.cons.genetype,file = "ResultsTable.mix.cons.genetype.csv")


keytypes(clariomdhumantranscriptcluster.db)

## Make volcano plots
library(ggplot2)
library(scales)

library(ggrepel)

library(RColorBrewer)

install.packages("ggallin")
library(ggallin)


install.packages("ggbreak")
library(ggbreak)

##Volcano plots for TD1 vs C1

## with Adjusted P value
# add a column of NAs
ResultsTable.mix.cons.all$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
ResultsTable.mix.cons.all$diffexpressed[ResultsTable.mix.cons.all$logFC > 0.6 & ResultsTable.mix.cons.all$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
ResultsTable.mix.cons.all$diffexpressed[ResultsTable.mix.cons.all$logFC < -0.6 & ResultsTable.mix.cons.all$adj.P.Val < 0.05] <- "DOWN"


p <- ggplot(data=ResultsTable.mix.cons.all, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed)) + geom_point() + theme_minimal()


ggsave("pbmc.vol.pdf")



# 1. by default, it is assigned to the categories in an alphabetical order):
p1 <- p + scale_color_manual(values=c("blue", "grey", "red"))+xlim(-6,6)

ggsave("p1.pbmc.pdf")


##Volcano plots for C1 vs C64

## with Adjusted P value
# add a column of NAs
C1.C64.ResultsTable.mix.cons.all$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
#C1.C64.ResultsTable.mix.cons.all$diffexpressed[C1.C64.ResultsTable.mix.cons.all$logFC > 0.6 & C1.C64.ResultsTable.mix.cons.all$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
#C1.C64.ResultsTable.mix.cons.all$diffexpressed[C1.C64.ResultsTable.mix.cons.all$logFC < -0.6 & C1.C64.ResultsTable.mix.cons.all$adj.P.Val < 0.05] <- "DOWN"


plot1.C1.C64<- ggplot(data=C1.C64.ResultsTable.mix.cons.all, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal()


ggsave("plot1.C1.C64_1.pdf")



# 1. by default, it is assigned to the categories in an alphabetical order):
plot2.C1.C64 <- plot1.C1.C64 + scale_color_manual(values=c("grey", "grey", "grey"))

ggsave("plot2.C1.C64.pdf")


##Volcano plots Treatment.vs recovery.C64vs TD1

## with Adjusted P value
# add a column of NAs
TD1.C64.ResultsTable.mix.cons.all$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
TD1.C64.ResultsTable.mix.cons.all$diffexpressed[TD1.C64.ResultsTable.mix.cons.all$logFC > 0.6 & TD1.C64.ResultsTable.mix.cons.all$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
TD1.C64.ResultsTable.mix.cons.all$diffexpressed[TD1.C64.ResultsTable.mix.cons.all$logFC < -0.6 & TD1.C64.ResultsTable.mix.cons.all$adj.P.Val < 0.05] <- "DOWN"



plot1.TD1.C64<- ggplot(data=TD1.C64.ResultsTable.mix.cons.all, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed)) + geom_point() + theme_minimal()+xlim(-6,6)


ggsave("plot1.TD1.C64.pdf")



# 1. by default, it is assigned to the categories in an alphabetical order):
plot2.TD1.C64 <- plot1.TD1.C64 + scale_color_manual(values=c("blue", "grey", "red"))

ggsave("plot2.TD1.C64.pdf")



## with Adjusted P value
# add a column of NAs
C64.TD1.ResultsTable.mix.cons.all$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
C64.TD1.ResultsTable.mix.cons.all$diffexpressed[C64.TD1.ResultsTable.mix.cons.all$logFC > 0.6 & C64.TD1.ResultsTable.mix.cons.all$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
C64.TD1.ResultsTable.mix.cons.all$diffexpressed[C64.TD1.ResultsTable.mix.cons.all$logFC < -0.6 & C64.TD1.ResultsTable.mix.cons.all$adj.P.Val < 0.05] <- "DOWN"



plot1.C64.TD1<- ggplot(data=C64.TD1.ResultsTable.mix.cons.all, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed)) + geom_point() + theme_minimal()+xlim(-6,6)


ggsave("plot1.C64.TD1.pdf")



# 1. by default, it is assigned to the categories in an alphabetical order):
plot2.C64.TD1 <- plot1.C64.TD1 + scale_color_manual(values=c("blue", "grey", "red"))

ggsave("plot2.C64.TD1.pdf")


##Volcano plot for ncRNA

## Treatment-Control+Recovery/2

df.vol.pbmc.ncRNA<-read.csv("select.pbmc.ncRNA.csv", header=TRUE,sep = ";", quote = "\"'", dec = ".")


## with Adjusted P value
# add a column of NAs
df.vol.pbmc.ncRNA$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df.vol.pbmc.ncRNA$diffexpressed[df.vol.pbmc.ncRNA$logFC > 0.6 & df.vol.pbmc.ncRNA$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df.vol.pbmc.ncRNA$diffexpressed[df.vol.pbmc.ncRNA$logFC < -0.6 & df.vol.pbmc.ncRNA$adj.P.Val < 0.05] <- "DOWN"


df.vol.pbmc.ncRNA$delabel <- NA
df.vol.pbmc.ncRNA$delabel[df.vol.pbmc.ncRNA$diffexpressed != "NO"] <- df.vol.pbmc.ncRNA$GENENAME[df.vol.pbmc.ncRNA$diffexpressed != "NO"]

the_plot_pbmc.ncRNA <- ggplot(data=df.vol.pbmc.ncRNA, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed,label=delabel)) + geom_point() + theme_minimal()+geom_text_repel(size = 2)+xlim(-3,3)


ggsave("the_plot_pbmc.ncRNA.pdf")

##Change the color

the_plot2_pbmc.ncRNA <- the_plot_pbmc.ncRNA + scale_color_manual(values=c("blue", "grey", "red"))


ggsave("the_plot2_pbmc.ncRNA.pdf")


## Control-Recovery

df.vol.pbmc.C1.C64.ncRNA<-read.csv("select.pbmc.C1.C64.ncRNA.csv", header=TRUE,sep = ";", quote = "\"'", dec = ".")

## with Adjusted P value
# add a column of NAs
df.vol.pbmc.C1.C64.ncRNA$diffexpressed <- "NO"


the_plot_PBMC.C1.C64.nc.RNA <- ggplot(data=df.vol.pbmc.C1.C64.ncRNA, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal()+xlim(-2,2)


ggsave("the_plot_PBMC.C1.C64.nc.RNA.pdf")


##Change the color

the_plot2_PBMC.C1.C64.nc.RNA <- the_plot_PBMC.C1.C64.nc.RNA + scale_color_manual(values=c("grey", "grey", "grey"))


ggsave("the_plot2_PBMC.C1.C64.nc.RNA.pdf")





## GO annotation with Cluster profiler
###library
library( ggplot2 )
library( scales )


###Install and load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
install.packages("wordcloud")

##Libraries for Clusterprofiler


library(viridis)

###Annotations.Load the annotation org.Hs.eg.db
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

library("org.Hs.eg.db")


library(clusterProfiler)
library(wordcloud)

library(viridis)



genes.pbmc.TD1.C1.clariom.r<-ResultsTable.mix.cons$ENTREZID


go_enrich.pbmc.TD1.C1.clariom <- enrichGO(gene = genes.pbmc.TD1.C1.clariom.r,
                                          OrgDb = 'org.Hs.eg.db', 
                                          keyType = 'ENTREZID',
                                          readable = T,
                                          ont = "BP",
                                          pvalueCutoff = 0.05, 
                                          qvalueCutoff = 0.10)



write.csv(go_enrich.pbmc.TD1.C1.clariom,file = "GO_genes.pbmc.TD1.C1.clariom.r.12.07.21")


go_enrich.pbmc.TD1.C1.clariom.simplify <- simplify(go_enrich.pbmc.TD1.C1.clariom, cutoff=0.72, by="p.adjust", select_fun=min)


write.csv(go_enrich.pbmc.TD1.C1.clariom.simplify,file = "Ggo_enrich.pbmc.TD1.C1.clariom.simplify")

dotplot(go_enrich.pbmc.TD1.C1.clariom.simplify,font.size=5,showCategory=50)

cnetplot(go_enrich.pbmc.TD1.C1.clariom.simplify, showCategory=c("telomere organization","telomere maintenance","histone phosphorylation","DNA damage response, signal transduction by p53 class mediator resulting in cell cycle arrest","aging",size=20))


##it is better to use vector because then all elements can be saved in the figure
p2<-cnetplot(go_enrich.pbmc.TD1.C1.clariom.simplify, showCategory=c("double-strand break repair","antigen processing and presentation of exogenous peptide antigen","regulation of cyclin-dependent protein serine/threonine kinase activity","regulation of DNA repair","positive regulation of neuron death","response to interferon-beta","type I interferon signaling pathway","response to interferon-alpha","regulation of innate immune response","defense response to symbiont","cellular response to reactive oxygen species","telomere organization","telomere maintenance","histone phosphorylation","DNA damage response, signal transduction by p53 class mediator resulting in cell cycle arrest","aging",size=20))

ggsave("cNetplot_PBMC_C1TD1.selected.GO.ageing_immunity.new2.pdf")




cnetplot(go_enrich.pbmc.TD1.C1.clariom.simplify, showCategory=c("chromosome segregation","nuclear division","mitotic nuclear division","sister chromatid segregation","nuclear chromosome segregation","DNA replication","regulation of mitotic cell cycle phase transition","protein-DNA complex assembly","negative regulation of cell cycle process","protein-DNA complex subunit organization","DNA packaging","microtubule cytoskeleton organization involved in mitosis","nucleosome assembly","regulation of sister chromatid segregation","cell cycle G2/M phase transition","G2/M transition of mitotic cell cycle","regulation of chromosome segregation","mitotic cell cycle checkpoint","nucleosome organization","meiotic cell cycle","chromatin assembly or disassembly","cell cycle G1/S phase transition","G1/S transition of mitotic cell cycle","regulation of chromosome organization","regulation of nuclear division","double-strand break repair","anaphase-promoting complex-dependent catabolic process","chromatin organization involved in regulation of transcription","negative regulation of gene expression, epigenetic","protein localization to chromosome","cytokinesis","establishment of chromosome localization","signal transduction by p53 class mediator","chromosome localization","antigen processing and presentation of exogenous peptide antigen","antigen processing and presentation of exogenous antigen","defense response to virus","defense response to symbiont","type I interferon signaling pathway","regulation of cyclin-dependent protein serine/threonine kinase activity","regulation of DNA replication","response to virus","G0 to G1 transition","kinetochore organization","telomere organization","DNA synthesis involved in DNA repair","histone phosphorylation","DNA replication initiation","response to interferon-alpha","DNA damage response, signal transduction by p53 class mediator resulting in cell cycle arrest","regulation of cell division","mitotic nuclear envelope disassembly","regulation of gene silencing","regulation of attachment of spindle microtubules to kinetochore","female meiotic nuclear division","regulation of microtubule-based process","membrane disassembly","telomere maintenance","regulation of DNA repair","negative regulation of viral genome replication","reciprocal homologous recombination","retrograde vesicle-mediated transport, Golgi to endoplasmic reticulum","viral genome replication","regulation of meiotic cell cycle","regulation of microtubule cytoskeleton organization","regulation of cellular amine metabolic process","regulation of response to DNA damage stimulus","spindle localization","interleukin-7-mediated signaling pathway","regulation of innate immune response","regulation of megakaryocyte differentiation","peptidyl-threonine phosphorylation","establishment of mitotic spindle localization","cellular response to reactive oxygen species","deoxyribonucleotide biosynthetic process","regulation of cellular amino acid metabolic process","female gamete generation","DNA damage response, detection of DNA damage","megakaryocyte differentiation","response to interleukin-7","cellular response to interleukin-7","response to drug","pyrimidine nucleoside triphosphate biosynthetic process","regulation of response to biotic stimulus","microtubule-based movement","microtubule depolymerization","response to UV","regulation of ubiquitin protein ligase activity","response to radiation","pyrimidine nucleoside triphosphate metabolic process","astrocyte activation","positive regulation of protein localization to nucleus","aging","covalent chromatin modification","response to estradiol","regulation of binding","folic acid-containing compound metabolic process","nucleobase-containing small molecule interconversion","positive regulation of neuron death","ciliary basal body-plasma membrane docking","response to reactive oxygen species","response to interferon-beta","liver regeneration",size=20))

p1<-cnetplot(go_enrich.pbmc.TD1.C1.clariom.simplify, showCategory=c("chromosome segregation","nuclear division","mitotic nuclear division","sister chromatid segregation","nuclear chromosome segregation","DNA replication","regulation of mitotic cell cycle phase transition","protein-DNA complex assembly","negative regulation of cell cycle process","protein-DNA complex subunit organization","DNA packaging","microtubule cytoskeleton organization involved in mitosis","nucleosome assembly","regulation of sister chromatid segregation","cell cycle G2/M phase transition","G2/M transition of mitotic cell cycle","regulation of chromosome segregation","mitotic cell cycle checkpoint","nucleosome organization","meiotic cell cycle","chromatin assembly or disassembly","cell cycle G1/S phase transition","G1/S transition of mitotic cell cycle","regulation of chromosome organization","regulation of nuclear division","double-strand break repair","anaphase-promoting complex-dependent catabolic process","chromatin organization involved in regulation of transcription","negative regulation of gene expression, epigenetic","protein localization to chromosome","cytokinesis","establishment of chromosome localization","signal transduction by p53 class mediator","chromosome localization","antigen processing and presentation of exogenous peptide antigen","antigen processing and presentation of exogenous antigen","defense response to virus","defense response to symbiont","type I interferon signaling pathway","regulation of cyclin-dependent protein serine/threonine kinase activity","regulation of DNA replication","response to virus","G0 to G1 transition","kinetochore organization","telomere organization","DNA synthesis involved in DNA repair","histone phosphorylation","DNA replication initiation","response to interferon-alpha","DNA damage response, signal transduction by p53 class mediator resulting in cell cycle arrest","regulation of cell division","mitotic nuclear envelope disassembly","regulation of gene silencing","regulation of attachment of spindle microtubules to kinetochore","female meiotic nuclear division","regulation of microtubule-based process","membrane disassembly","telomere maintenance","regulation of DNA repair","negative regulation of viral genome replication","reciprocal homologous recombination","retrograde vesicle-mediated transport, Golgi to endoplasmic reticulum","viral genome replication","regulation of meiotic cell cycle","regulation of microtubule cytoskeleton organization","regulation of cellular amine metabolic process","regulation of response to DNA damage stimulus","spindle localization","interleukin-7-mediated signaling pathway","regulation of innate immune response","regulation of megakaryocyte differentiation","peptidyl-threonine phosphorylation","establishment of mitotic spindle localization","cellular response to reactive oxygen species","deoxyribonucleotide biosynthetic process","regulation of cellular amino acid metabolic process","female gamete generation","DNA damage response, detection of DNA damage","megakaryocyte differentiation","response to interleukin-7","cellular response to interleukin-7","response to drug","pyrimidine nucleoside triphosphate biosynthetic process","regulation of response to biotic stimulus","microtubule-based movement","microtubule depolymerization","response to UV","regulation of ubiquitin protein ligase activity","response to radiation","pyrimidine nucleoside triphosphate metabolic process","astrocyte activation","positive regulation of protein localization to nucleus","aging","covalent chromatin modification","response to estradiol","regulation of binding","folic acid-containing compound metabolic process","nucleobase-containing small molecule interconversion","positive regulation of neuron death","ciliary basal body-plasma membrane docking","response to reactive oxygen species","response to interferon-beta","liver regeneration",size=20))


ggsave("cNetplot_PBMC_C1TD1.all.GO.new.pdf")



emapplot.pbmc<-emapplot(go_enrich.pbmc.TD1.C1.clariom.simplify,showCategory=c("chromosome segregation","nuclear division","mitotic nuclear division","sister chromatid segregation","nuclear chromosome segregation","DNA replication","regulation of mitotic cell cycle phase transition","protein-DNA complex assembly","negative regulation of cell cycle process","protein-DNA complex subunit organization","DNA packaging","microtubule cytoskeleton organization involved in mitosis","nucleosome assembly","regulation of sister chromatid segregation","cell cycle G2/M phase transition","G2/M transition of mitotic cell cycle","regulation of chromosome segregation","mitotic cell cycle checkpoint","nucleosome organization","meiotic cell cycle","chromatin assembly or disassembly","cell cycle G1/S phase transition","G1/S transition of mitotic cell cycle","regulation of chromosome organization","regulation of nuclear division","double-strand break repair","anaphase-promoting complex-dependent catabolic process","chromatin organization involved in regulation of transcription","negative regulation of gene expression, epigenetic","protein localization to chromosome","cytokinesis","establishment of chromosome localization","signal transduction by p53 class mediator","chromosome localization","antigen processing and presentation of exogenous peptide antigen","antigen processing and presentation of exogenous antigen","defense response to virus","defense response to symbiont","type I interferon signaling pathway","regulation of cyclin-dependent protein serine/threonine kinase activity","regulation of DNA replication","response to virus","G0 to G1 transition","kinetochore organization","telomere organization","DNA synthesis involved in DNA repair","histone phosphorylation","DNA replication initiation","response to interferon-alpha","DNA damage response, signal transduction by p53 class mediator resulting in cell cycle arrest","regulation of cell division","mitotic nuclear envelope disassembly","regulation of gene silencing","regulation of attachment of spindle microtubules to kinetochore","female meiotic nuclear division","regulation of microtubule-based process","membrane disassembly","telomere maintenance","regulation of DNA repair","negative regulation of viral genome replication","reciprocal homologous recombination","retrograde vesicle-mediated transport, Golgi to endoplasmic reticulum","viral genome replication","regulation of meiotic cell cycle","regulation of microtubule cytoskeleton organization","regulation of cellular amine metabolic process","regulation of response to DNA damage stimulus","spindle localization","interleukin-7-mediated signaling pathway","regulation of innate immune response","regulation of megakaryocyte differentiation","peptidyl-threonine phosphorylation","establishment of mitotic spindle localization","cellular response to reactive oxygen species","deoxyribonucleotide biosynthetic process","regulation of cellular amino acid metabolic process","female gamete generation","DNA damage response, detection of DNA damage","megakaryocyte differentiation","response to interleukin-7","cellular response to interleukin-7","response to drug","pyrimidine nucleoside triphosphate biosynthetic process","regulation of response to biotic stimulus","microtubule-based movement","microtubule depolymerization","response to UV","regulation of ubiquitin protein ligase activity","response to radiation","pyrimidine nucleoside triphosphate metabolic process","astrocyte activation","positive regulation of protein localization to nucleus","aging","covalent chromatin modification","response to estradiol","regulation of binding","folic acid-containing compound metabolic process","nucleobase-containing small molecule interconversion","positive regulation of neuron death","ciliary basal body-plasma membrane docking","response to reactive oxygen species","response to interferon-beta","liver regeneration"))


##Making emapplot


install.packages("ggnewscale")
library(ggnewscale)

d <- GOSemSim::godata("org.Hs.eg.db", ont = "BP")    
compare_cluster_GO_emap <- enrichplot::pairwise_termsim(go_enrich.pbmc.TD1.C1.clariom.simplify, semData = d,  method="Wang")
emapplot(compare_cluster_GO_emap)

##Web doc:https://rdrr.io/github/GuangchuangYu/enrichplot/man/emapplot.html
##Layout of the map, e.g. 'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'.
ema1.new.star.pbmc.no.grp<-emapplot(compare_cluster_GO_emap,coords = NULL,
                                    color = "p.adjust",
                                    min_edge = 0.2,
                                    cex_label_category = 1,
                                    cex_category = 1,
                                    cex_line = 1,
                                    shadowtext = TRUE,
                                    label_style = "shadowtext",
                                    repel = FALSE,
                                    node_label = "all",
                                    with_edge = TRUE,
                                    group_category = FALSE,
                                    group_legend = FALSE,
                                    cex_label_group = 1,
                                    nWords = 4,
                                    label_format = 30,
                                    clusterFunction = stats::kmeans,
                                    nCluster = NULL,pie_scale=1.5,layout="star",showCategory=104)


ggsave("ema1.new.star.pbmc.no.grp.cat.pdf")


ema2 <- emapplot(compare_cluster_GO_emap, pie_scale=1.5)
ema3 <- emapplot(compare_cluster_GO_emap,layout="kk")

ggsave("ema3.pdf")
ema4 <- emapplot(compare_cluster_GO_emap, pie_scale=1.5,layout="kk",showCategory=104)


cowplot::plot_grid(ema1, ema2, ema3, ema4, ncol=2, labels=LETTERS[1:4])


ggsave("ema4.pdf")


cnetplot(go_enrich.pbmc.TD1.C1.clariom.simplify, showCategory=c("cellular response to reactive oxygen species","telomere organization","telomere maintenance","histone phosphorylation","DNA damage response, signal transduction by p53 class mediator resulting in cell cycle arrest","aging",size=20))




##cNetplot for all categories from simplified..

cnetplot(go_enrich.pbmc.TD1.C1.clariom.simplify,font.size=3)


## Dot plot with all from simplify
dotplot(go_enrich.pbmc.TD1.C1.clariom.simplify,font.size=3,showCategory=104)+geom_point(size=0.1)

dotplot(go_enrich.pbmc.TD1.C1.clariom.simplify,font.size=3,showCategory=104)


##To change the dot size +scale_size(range = c(0.5, 2.5))
dotplot(go_enrich.pbmc.TD1.C1.clariom.simplify,font.size=3,showCategory=104)+scale_size(range = c(0.5, 2.5))


##GO for non-coding RNA

df.pbmc.miRNA<-read.csv("miRNA.PBMC.csv", header=TRUE,sep = ";", quote = "\"'", dec = ".")


genes.pbmc.miRNA<-df.pbmc.miRNA$ENTREZID


go_enrich.pbmc.miRNA <- enrichGO(gene = genes.pbmc.miRNA,
                                 OrgDb = 'org.Hs.eg.db', 
                                 keyType = 'ENTREZID',
                                 readable = T,
                                 ont = "BP",
                                 pvalueCutoff = 0.05, 
                                 qvalueCutoff = 0.10)



write.csv(go_enrich.pbmc.miRNA,file = "go_enrich.pbmc.miRNA.csv")

write.csv(go_enrich.pbmc.miRNA,file = "go_enrich.pbmc.miRNA")

go_enrich.pbmc.miRNA.simplify <- simplify(go_enrich.pbmc.miRNA, cutoff=0.70, by="p.adjust", select_fun=min)


write.csv(go_enrich.pbmc.miRNA.simplify,file = "go_enrich.pbmc.miRNA.simplify.csv")
write.csv(go_enrich.pbmc.miRNA.simplify,file = "go_enrich.pbmc.miRNA.simplify")


dotplot(go_enrich.pbmc.miRNA.simplify,font.size=3,showCategory=67)+scale_size(range = c(0.5, 2.5))


##Coding and non coding selected

df.pbmc.coding.ncRNA<-read.csv("PBMC_coding_ncRNA.csv", header=TRUE,sep = ";", quote = "\"'", dec = ".")



genes.pbmc.coding.ncRNA<-df.pbmc.coding.ncRNA$ENTREZID


go_enrich.pbmc.coding.ncRNA <- enrichGO(gene = genes.pbmc.coding.ncRNA,
                                        OrgDb = 'org.Hs.eg.db', 
                                        keyType = 'ENTREZID',
                                        readable = T,
                                        ont = "BP",
                                        pvalueCutoff = 0.05, 
                                        qvalueCutoff = 0.10)

write.csv(go_enrich.pbmc.coding.ncRNA,file = "go_enrich.pbmc.coding.ncRNA")


go_enrich.pbmc.coding.ncRNA.simplify <- simplify(go_enrich.pbmc.coding.ncRNA, cutoff=0.70, by="p.adjust", select_fun=min)



write.csv(go_enrich.pbmc.coding.ncRNA.simplify,file = "go_enrich.pbmc.coding.ncRNA.simplify")


cnetplot(go_enrich.pbmc.coding.ncRNA.simplify, showCategory=c("cell cycle G1/S phase transition","G1/S transition of mitotic cell cycle","negative regulation of G1/S transition of mitotic cell cycle","negative regulation of cell cycle G1/S phase transition","regulation of cyclin-dependent protein serine/threonine kinase activity","negative regulation of blood vessel endothelial cell proliferation involved in sprouting angiogenesis","sprouting angiogenesis","negative regulation of translation","negative regulation of cellular amide metabolic process",size=20))


##Network analysis with WGCNA

BiocManager::install("WGCNA")

library(WGCNA)

#Install GO database for the species.The organism-speci_c packages have names of the form org.Xx.eg.db, where Xx stands for organism code, for example, Mm for mouse, Hs for human, etc

BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)

##
GOenr = GOenrichmentAnalysis(genes.pbmc.TD1.C1.clariom.r, organism = "human", nBestP = 10)


GOenrichmentAnalysis(genes.pbmc.TD1.C1.clariom.r,
                     yeastORFs = NULL,
                     organism = "human", 
                     ontologies = c("BP", "CC", "MF"), 
                     evidence = "all",
                     includeOffspring = TRUE, 
                     backgroundType = "givenInGO",
                     removeDuplicates = TRUE,
                     leaveOutLabel = NULL,
                     nBestP = 10, pCut = NULL, 
                     nBiggest = 0, 
                     getTermDetails = TRUE,
                     verbose = 2, indent = 0)

GOenrichmentAnalysis(labels, 
                     entrezCodes, 
                     yeastORFs = NULL,
                     organism = "human", 
                     ontologies = c("BP", "CC", "MF"), 
                     evidence = "all",
                     includeOffspring = TRUE, 
                     backgroundType = "givenInGO",
                     removeDuplicates = TRUE,
                     leaveOutLabel = NULL,
                     nBestP = 10, pCut = NULL, 
                     nBiggest = 0, 
                     getTermDetails = TRUE,
                     verbose = 2, indent = 0)

###pd.clariom.d.human


transcript.eset.affy <- annotateEset(transcript.eset1, pd.clariom.d.human)


#Work with Limma
#create design matrix
DesignMatrix1<-matrix(c(0,1,0,1,0,0,0,0,1,0,1,0,1,0,0,0,0,1),nrow=6, ncol=3,byrow = TRUE)

colnames(DesignMatrix1)<-c("Treatment","Control","Recovery")
rownames(DesignMatrix1)<-sampleNames(transcript.eset.affy)



#create contrast matrix
FittedModel<-lmFit(transcript.eset.affy, DesignMatrix1)

##Mixed contrast matrix
mMxContrastsMatrix1<-makeContrasts(contrasts="Treatment-(Control+Recovery)/2",levels=DesignMatrix1)

FitContrast.mix1=contrasts.fit(FittedModel,mMxContrastsMatrix1)
efitContrast.mix1=eBayes(FitContrast.mix1)

ResultsTable.mix.cons1<-topTable(efitContrast.mix1, number =
                                   100000,adjust.method = "BH",p.value=0.05)




## Make general heatmap


install.packages("pheatmap")
library(pheatmap)

df.PBMC.htmap<-read.csv("ResultsTable.mix.cons.csv", header=TRUE,sep = ";", quote = "\"'", dec = ".")


df.PBMC.htmap.matrix = as.matrix(df.PBMC.htmap[,18:20])

pheatmap(df.PBMC.htmap.matrix,color=colorRampPalette(c("navy", "white", "red"))(50),fontsize=8,border_color=NA,cluster_cols = T,main = "PBMC")


## Make annotation heatmap 30.07.21

df.pbmc.annot.htmap<-read.csv("df.pbmc.annot.heatmap.new.csv",row.names = 1, header=TRUE,sep = ";", quote = "\"'", dec = ".")

df.pbmc.annot.htmap.matrix = as.matrix(df.pbmc.annot.htmap[,2:4])

GO.pbmc.df = data.frame("GO" = df.pbmc.annot.htmap$GO)
rownames(GO.pbmc.df) = rownames(df.pbmc.annot.htmap.matrix)


newCols <- colorRampPalette(grDevices::rainbow(length(unique(GO.pbmc.df$GO))))
mycolors.pbmc <- newCols(length(unique(GO.pbmc.df$GO)))
names(mycolors.pbmc) <- unique(GO.pbmc.df$GO)
mycolors.pbmc <- list(GO = mycolors.pbmc)



pheatmap(df.pbmc.annot.htmap.matrix,
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_row=c(10,13,14,20,21,25,28,30,36,37,45,56,59,61,62,63,66,67,68,69,70,72,73,74,76,77,79,80,81,82,85,89,92,93,94),
         gaps_col=c(1,2),
         cellheight = 6,
         cellwidth = 50,
         border_color=NA,
         fontsize_row = 6,
         main="PBMC/genes grouped by GO",
         filename = "Annotation heatmap.PBMC.row.col.gap.30.07.21.pdf",
         annotation_row = GO.pbmc.df,
         annotation_colors = mycolors.pbmc)











##old one


df.pbmc.annot.htmap<-read.csv("df.pbmc.annot.heatmap.csv",row.names = 1, header=TRUE,sep = ";", quote = "\"'", dec = ".")

df.pbmc.annot.htmap.matrix = as.matrix(df.pbmc.annot.htmap[,2:4])

GO.pbmc.df = data.frame("GO" = df.pbmc.annot.htmap$GO)
rownames(GO.pbmc.df) = rownames(df.pbmc.annot.htmap.matrix)


newCols <- colorRampPalette(grDevices::rainbow(length(unique(GO.pbmc.df$GO))))
mycolors.pbmc <- newCols(length(unique(GO.pbmc.df$GO)))
names(mycolors.pbmc) <- unique(GO.pbmc.df$GO)
mycolors.pbmc <- list(GO = mycolors.pbmc)


pheatmap(df.pbmc.annot.htmap.matrix,
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_row=c(3,11,12,13,16,19,21,22,30),
         gaps_col=1,
         cellheight = 6,
         cellwidth = 50,
         border_color=NA,
         fontsize_row = 6,
         main="PBMC/genes grouped by GO",
         filename = "Annotation heatmap.PBMC.pdf",
         annotation_row = GO.pbmc.df,
         annotation_colors = mycolors.pbmc)


pheatmap(df.pbmc.annot.htmap.matrix,
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col=1,
         cellheight = 6,
         cellwidth = 50,
         border_color=NA,
         fontsize_row = 6,
         main="PBMC/genes grouped by GO",
         filename = "Annotation heatmap.PBMC.no.break.pdf",
         annotation_row = GO.pbmc.df,
         annotation_colors = mycolors.pbmc)



## Heatmap for miRNA

df.miRNA.PBMC<-read.csv("miRNA.PBMC.csv", header=TRUE,sep = ";", quote = "\"'", dec = ".",row.names = 1)


df.miRNA.PBMC.matrix = as.matrix(df.miRNA.PBMC[,5:7])

pheatmap(df.miRNA.PBMC.matrix,color=colorRampPalette(c("navy", "white", "red"))(50),fontsize=8,border_color=NA,cluster_cols = T,main = "PBMC/miRNA")


##To check the session info and all libraries and versions
print(sessionInfo())



#How to extract info from pd.clariom.d.human database? It is a
#sql-lite it seems.
#Try to print all info in all tables in whole database but NO
#Genesymbol/annotation/entrez info present. Unexpected.

DBTables<-dbListTables(pd.clariom.d.human@getdb())
for(Table in DBTables)
{
  print(Table)
  print(dbListFields(pd.clariom.d.human@getdb(), Table))
}
DBTables<-dbListTables(pd.clariom.d.human@getdb())
for(Table in DBTables)
{
  print(Table)
  print(dbReadTable(pd.clariom.d.human@getdb(), Table)[153713:153733,])
}


pd.clariom.d.human@tableInfo

featureSet.table<- dbReadTable(pd.clariom.d.human@getdb(), "featureSet")

type_dict.table<- dbReadTable(pd.clariom.d.human@getdb(), "type_dict")

nonNas<-which(!is.na(type_dict.table[,"type"]))

nonNas<-which(!is.na(featureSet.table[,"type"]))


featureSet.table[nonNas,]


featureSet.table[nonNas,"transcript_cluster_id"]

##Check ncRNA from NCBI

install.packages("rentrez")
library('rentrez')
result<-entrez_summary(db = "gene", id = 102465452)
result[["description"]]
Link <- entrez_link(dbfrom = "gene", id = 102465452, db = "nucleotide")
genlinkid <- Link[["links"]][["gene_nuccore_refseqrna"]]
result<-entrez_summary(db = "nucleotide", id = genlinkid)
result$slen


##Read ID
All_ncRNA_NCBI_ID_PBMC <- read.csv("All_ncRNA_NCBI_ID_PBMC.csv", header=TRUE,sep = ";", quote = "\"'", dec = ".")

result<-entrez_summary(db = "gene", id = (All_ncRNA_NCBI_ID_PBMC$ids))
result[["description"]]
Link <- entrez_link(dbfrom = "gene", id = All_ncRNA_NCBI_ID_PBMC, db = "nucleotide")
genlinkid <- Link[["links"]][["gene_nuccore_refseqrna"]]
result<-entrez_summary(db = "nucleotide", id = genlinkid)
result$slen

















