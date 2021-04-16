library(ArrayExpress)
library("affy")
library("limma")
library(DataExplorer)


AEset1 = ArrayExpress("E-MEXP-3936")
AEsetnorm1 = rma(AEset1)
fac1 = grep("Factor.Value",colnames(pData(AEsetnorm1)), value=T)
if (suppressWarnings(require("arrayQualityMetrics", quietly=TRUE))) {
qanorm = arrayQualityMetrics(AEsetnorm1,outdir = "QAnorm",intgroup = fac)}

express1 = exprs(AEsetnorm1)

qx <- as.numeric(quantile(express1,c(0.0,0.25,0.5,0.75,0.99,1.0),na.rm=T))
LogC <- (qx[5] > 100 || (qx[6] - qx[1] > 50 && qx[2] > 0))
if(LogC) { express1[which (express1<= 0)] <- NaN
exprs(AEsetnorm1) <- log2(express1)}    #Log transformation of exprs(gset) for DE gene analysis


# group membership for all samples
gsms <- paste0("0X10X10X10X10X10X10X10X10X10X10X10X10X1")
#gsms <- paste0("X01X01X01X01X01X01X01X01X01X01X01X01X01")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
AEsetnorm1 <- AEsetnorm1[ ,sel]


# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("BaseLine","SleepRestriction"))
levels(gs) <- groups
AEsetnorm1$group <- gs
design <- model.matrix(~group + 0, AEsetnorm1)
colnames(design) <- levels(gs)

fit <- lmFit(AEsetnorm1, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

#tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))

de_genes <-tT[abs(tT$logFC) < 1.2,] #Extract genes which have logFC value < 0.2
#de_genes <- tT[tT$P.Value<0.05,]
de_genes <-de_genes[de_genes$P.Val < 0.05,] #Extract genes which have adjPval < 0.05
#Here we got 792 differentially expressed genes


