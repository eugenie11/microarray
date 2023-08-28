# Loading libraries
library(limma)
library(affy)
library(ggplot2)
#library(promor)
library(affyio)
library(arrayQualityMetrics)
library(gcrma)
library(hgu133plus2.db)
library(gplots)
#library(aheatmap)
library(RColorBrewer)

#-----------Looking at the data
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31524

# specify the path on your computer where the folder that contains the CEL-files is located
celpath = "C:/Users/12har/OneDrive/Desktop/H9N2/8H_10H"
#celpath = "C:/Users/12har/Desktop/H9N2/8H_10H"

# import CEL files containing raw probe-level data into an R AffyBatch object
data = ReadAffy(celfile.path=celpath)

# obtain intensity data
expr = exprs(data)
expr[1:5,]
pm(data)[1:5,]

# retrieve the sample annotation of the data
ph = data@phenoData # pData(data)
ph

# retrieve the probe annotation of the data
feat = data@featureData
feat

# retrieve the experiment annotation
exp = data@experimentData
exp
featureNames(data)

#---------Quality control of the data
ph@data[ ,1] = c("8_1","8_2","8_3","10_1","10_2","10_3")
ph

# Create plots to assess the quality of the data
arrayQualityMetrics(expressionset = data,
                    outdir = "Report_for_raw_data",
                    force = TRUE,
                    do.logtransform = TRUE)

# normalizing the data using rma
arrayQualityMetrics(expressionset = rma(data),
                    outdir = "Report_for_rma_data",
                    force = TRUE)
data.matrix = exprs(rma(data))

# normalizing the data using gcrma
arrayQualityMetrics(expressionset = gcrma(data),
                    outdir = "Report_for_gcrma_data",
                    force = TRUE)

#-------adding gene annotation to normalized expression output
# Strategy is to create data frame objects and merge them together - put expression info into a data frame
my_frame <- data.frame(exprs(rma(data)))

Annot <- data.frame(ACCNUM=sapply(contents(hgu133plus2ACCNUM), paste, collapse=", "), 
                    SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "), 
                    DESC=sapply(contents(hgu133plus2GENENAME), paste, collapse=", "))

# Merge data frames together (like a database table join)
data_with_gene_ann <- merge(Annot, my_frame, by.x=0, by.y=0, all=T)
temp <- data_with_gene_ann[,-1]
rownames(temp) <- data_with_gene_ann[,1]

# Write out to a file:
#write.table(all,file="data.ann.txt",sep="\t")


#-------Identification of DE genes
ph@data[ ,2] = c("8H","8H","8H","10H","10H","10H")
colnames(ph@data)[2]="Hour"
ph@data[ ,3] = c("1","2","3","1","2","3")
colnames(ph@data)[3]="Rep"
ph@data

groupsR = ph@data$Rep
groupsH = ph@data$Hour
fr = factor(groupsR,levels=c("1","2","3"))
fh = factor(groupsH,levels=c("8H","10H"))

paired.design = model.matrix(~ fr + fh)
colnames(paired.design)=c("Intercept","2vs1","3vs1","8Hvs10H")
data.fit = lmFit(rma(data),paired.design)
data.fit$coefficients

data.fit.eb = eBayes(data.fit)
#data.fit.eb$p.value[,1:20]

# decide on the number of DE genes to select using volcano plot
name = "Volcano.jpg"
jpeg(name)
volcanoplot(data.fit.eb,coef=4,highlight=39)
dev.off()

options(digits=2)
tab = topTable(data.fit.eb,coef=4,number=200,adjust.method='BH')

topgenes = tab[tab[, "adj.P.Val"] < 0.005, ]
dim(topgenes) # 39 genes
topups = topgenes[topgenes[, "logFC"] > 1, ]
dim(topups)
topdowns = topgenes[topgenes[, "logFC"] < -1, ]
dim(topdowns)

data.matrix.de = data.matrix[(rownames(topgenes)),]

colnames(data.matrix.de) <- c('8H.1','8H.2','8H.3','10H.1','10H.2','10H.3')
heatmap.2(data.matrix.de, col=rev(brewer.pal(9,"RdBu")))
dev.off()

temp["DESC"][(rownames(data.matrix.de)),]


#--------funtional annotation--------GO

##Bimap interface: 
x<-hgu133plus2GO 
#Get the manufacturer identifiers that are mapped to a GO ID 
mapped_genes<-mappedkeys(x) 
#Convert to a list 
xx<-as.list(x[mapped_genes][1:300]) 
if(length(xx)>0){ 
  #Try the first one 
  got<-xx[[1]] 
  got[[1]][["GOID"]] 
  got[[1]][["Ontology"]] 
  got[[1]][["Evidence"]] }


# volcano plot
total_tab = topTable(data.fit.eb,coef=4,number=54675,adjust.method='BH')
de=total_tab
# add a column of NAs
de$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$logFC > 0.85 & de$P.Value < 0.005] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$logFC < -0.85 & de$P.Value < 0.005] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=de, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.85, 0.85), col="red") +
  geom_hline(yintercept=-log10(0.005), col="red")

## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]
de$delabel[de$diffexpressed == "UP"] <- rownames(de)
de$delabel[de$diffexpressed == "DOWN"] <- rownames(de)
library(ggrepel)
ggplot(data=de, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.85, 0.85), col="red")


dim(de$delabel[de$diffexpressed != "NO"])
de$adj.P.Val[de$diffexpressed == "DOWN"]  # 0.0093














