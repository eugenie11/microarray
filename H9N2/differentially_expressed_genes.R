library(limma)
library(affy)
# specify the path on your computer where the folder that contains the CEL-files is located
celpath = "C:/Users/12har/OneDrive/Desktop/H9N2/8H_10H"

# import CEL files containing raw probe-level data into an R AffyBatch object
data = ReadAffy(celfile.path=celpath)

#Normalizing data with RMA
data.rma = rma(data)
data.matrix = exprs(data.rma)

#Retrieve sample annotation using affy
ph = data@phenoData
ph$sample
feat = data@featureData
exp = data@experimentData

#Retrieving other types of information
cdfName(data)
featureNames(data)
length(featureNames(data))
length(probeNames(data))

### for two groups of paired data
ph@data[ ,2] = c("8H","8H","8H","10H","10H","10H")
colnames(ph@data)[2]="Time"
ph@data[ ,3] = c("1","2","3","1","2","3")
colnames(ph@data)[3]="Rep"
ph@data

groupsT = ph@data$Time
groupsR = ph@data$Rep 

ft = factor(groupsT,levels=c("8H","10H"))
fr = factor(groupsR,levels=c("1","2","3"))

paired.design = model.matrix(~ ft + fr)
colnames(paired.design)=c("Intercept","2vs1","3vs1","8Hvs10H")
data.fit = lmFit(data.rma,paired.design)
data.fit$coefficients

data.fit.eb = eBayes(data.fit)
data.fit.eb$p.value[,1:23]
############

options(digits=2)
tab = topTable(data.fit.eb,coef=2,number=200,adjust.method="BH")

#Select set of genes with adjusted p-values below threshold
topgenes = tab[tab[, "adj.P.Val"] < 0.01, ]
dim(topgenes)

topups = topgenes[topgenes[, "logFC"] > 1, ]
dim(topups)
topdowns = topgenes[topgenes[, "logFC"] < -0.8, ] # Adjusted from -1 to -0.8
dim(topdowns)

#Volcano plot
name = "Volcano.jpg"
jpeg(name)
volcanoplot(data.fit.eb,coef=2,highlight=10)
#dev.off()

#Getting up regulated, down regulated genes
IDs.up = rownames(topups)
IDs.down = rownames(topdowns)

data.matrix.up = data.matrix[(rownames(topups)),]
sampleNames = vector()
featureNames = vector()
heatlogs = vector()
for (i in 1:6)
{
  sampleNames = c(sampleNames,rep(ph@data[i,1],dim(topdowns)[1]))
  featureNames = c(featureNames,rownames(data.matrix.up[1:dim(topdowns)[1],]))
  heatlogs = c(heatlogs,data.matrix.up[1:dim(topdowns)[1],i])
}
heatData = data.frame(norm_logInt=heatlogs,sampleName=sampleNames,featureName=featureNames)
dataHeat = ggplot(heatData, aes(sampleName,featureName))
dataHeat + geom_tile(aes(fill=norm_logInt)) + scale_fill_gradient(low="green", high="red")
