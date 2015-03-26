########################################################################
# Transcriptome. Differential expression. limma                        #
# Diana Coman Schmid                                                   #
# Eawag 2015                                                           #
# diana.comanschmid@eawag.ch                                           #
########################################################################

library(limma)

# read in the normalized expression data (limma)
eset <- read.table("/.../zebrafishNorm.txt",header=TRUE, row.names=1,sep = '\t')

# read in the sample description
targets = readTargets("/.../targets.txt",sep="")

# limma differential expression
f = factor(targets$Group, levels = unique(targets$Group))
design = model.matrix(~0 + f)
colnames(design) = levels(f)
rownames(design) <- colnames(eset)
fit <- lmFit(eset, design)

# make contrast matrix for all chemicals vs. control
contrast.matrix <- makeContrasts(DEP.h.SC=(DEPHighchem - DEPPermSolvC)-(DEPPermSolvC - DEPPermCtrl),
                                 DEP.l.SC=(DEPLowchem - DEPPermSolvC)-(DEPPermSolvC - DEPPermCtrl),
                                 DEP.h=DEPHighchem - DEPPermSolvC,
                                 DEP.l=DEPLowchem - DEPPermSolvC,
                                 ctrl.DEP=DEPPermSolvC - DEPPermCtrl,
                                 
                                 DMP.h.SC=(DMPHighchem - DMPDinoSolvC)-(DMPDinoSolvC - DMPDinoCtrl),
                                 DMP.l.SC=(DMPLowchem - DMPDinoSolvC)-(DMPDinoSolvC - DMPDinoCtrl),
                                 DMP.h=DMPHighchem - DMPDinoSolvC,
                                 DMP.l=DMPLowchem - DMPDinoSolvC,
                                 ctrl.DMP=DMPDinoSolvC - DMPDinoCtrl,
                                 
                                 SDS.h.SC=(SDSHighchem - SDSPCPSolvC)-(SDSPCPSolvC - SDSPCPCtrl),
                                 SDS.l.SC=(SDSLowchem - SDSPCPSolvC)-(SDSPCPSolvC - SDSPCPCtrl),
                                 SDS.h=SDSHighchem - SDSPCPSolvC,
                                 SDS.l=SDSLowchem - SDSPCPSolvC,
                                 ctrl.SDS=SDSPCPSolvC - SDSPCPCtrl,
				 
				 DCP.h.SC=(DCPHighchem - DCPFlucSolvC)-(DCPFlucSolvC - DCPFlucCtrl),
                                 DCP.l.SC=(DCPLowchem - DCPFlucSolvC)-(DCPFlucSolvC - DCPFlucCtrl),
                                 DCP.h=DCPHighchem - DCPFlucSolvC,
                                 DCP.l=DCPLowchem - DCPFlucSolvC,
                                 ctrl.DCP=DCPFlucSolvC - DCPFlucCtrl,
                                 
                                 PCP.h.SC=(PCPHighchem - SDSPCPSolvC)-(SDSPCPSolvC - SDSPCPCtrl),
                                 PCP.l.SC=(PCPLowchem - SDSPCPSolvC)-(SDSPCPSolvC- SDSPCPCtrl),
                                 PCP.h=PCPHighchem - SDSPCPSolvC,
                                 PCP.l=PCPLowchem - SDSPCPSolvC,
                                 ctrl.PCP=SDSPCPSolvC- SDSPCPCtrl,
                                 
				 Dino.h.SC=(DinoHighchem - DMPDinoSolvC)-(DMPDinoSolvC - DMPDinoCtrl),
                                 Dino.l.SC=(DinoLowchem - DMPDinoSolvC)-(DMPDinoSolvC - DMPDinoCtrl),
                                 Dino.h=DinoHighchem - DMPDinoSolvC,
                                 Dino.l=DinoLowchem - DMPDinoSolvC,
                                 ctrl.Dino=DMPDinoSolvC - DMPDinoCtrl,
                                 
				 Acro.h.SC=(AcroHighchem - AcroEsfenSolvC)-(AcroEsfenSolvC - AcroEsfenCtrl),
                                 Acro.l.SC=(AcroLowchem - AcroEsfenSolvC)-(AcroEsfenSolvC - AcroEsfenCtrl),
                                 Acro.h=AcroHighchem - AcroEsfenSolvC,
                                 Acro.l=AcroLowchem - AcroEsfenSolvC,
                                 ctrl.Acro=AcroEsfenSolvC - AcroEsfenCtrl,
                                 
				 AA.h.SC=(AAHighchem - AALindSolvC)-(AALindSolvC - AALindCtrl),
                                 AA.l.SC=(AALowchem - AALindSolvC)-(AALindSolvC - AALindCtrl),
                                 AA.h=AAHighchem - AALindSolvC,
                                 AA.l=AALowchem - AALindSolvC,
                                 ctrl.AA=AALindSolvC - AALindCtrl,
                                 
				 Fluc.h.SC=(FlucHighchem - DCPFlucSolvC)-(DCPFlucSolvC - DCPFlucCtrl),
                                 Fluc.l.SC=(FlucLowchem - DCPFlucSolvC)-(DCPFlucSolvC - DCPFlucCtrl),
                                 Fluc.h=FlucHighchem - DCPFlucSolvC,
                                 Fluc.l=FlucLowchem - DCPFlucSolvC,
                                 ctrl.Fluc=DCPFlucSolvC - DCPFlucCtrl,
                                 
                                 Perm.h.SC=(PermHighchem - DEPPermSolvC)-(DEPPermSolvC - DEPPermCtrl),
                                 Perm.l.SC=(PermLowchem - DEPPermSolvC)-(DEPPermSolvC - DEPPermCtrl),
                                 Perm.h=PermHighchem - DEPPermSolvC,
                                 Perm.l=PermLowchem - DEPPermSolvC,
                                 ctrl.Perm=DEPPermSolvC - DEPPermCtrl,
                                 
				 Esfen.h.SC=(EsfenHighchem - AcroEsfenSolvC)-(AcroEsfenSolvC - AcroEsfenCtrl),
                                 Esfen.l.SC=(EsfenLowchem - AcroEsfenSolvC)-(AcroEsfenSolvC - AcroEsfenCtrl),
                                 Esfen.h=EsfenHighchem - AcroEsfenSolvC,
                                 Esfen.l=EsfenLowchem - AcroEsfenSolvC,
                                 ctrl.Esfen=AcroEsfenSolvC - AcroEsfenCtrl,
                                 
                                 Lind.h.SC=(LindHighchem - AALindSolvC)-(AALindSolvC - AALindCtrl),
                                 Lind.l.SC=(LindLowchem - AALindSolvC)-(AALindSolvC - AALindCtrl),
                                 Lind.h=LindHighchem - AALindSolvC,
                                 Lind.l=LindLowchem - AALindSolvC,
                                 ctrl.Lind=AALindSolvC - AALindCtrl,
				 
                                 levels=design)                                 

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# store limma differential expression results as list
diffE <- list()
for (d in 1:ncol(fit2$contrasts)){
  diffE[[d]] <- topTable(fit2,coef=colnames(fit2$contrasts)[d],number=1000000,adjust="fdr")
  diffE[[d]]$locusID <- rownames(diffE[[d]])
  names(diffE)[d] <- colnames(fit2$contrasts)[d]
}

# select DMSO controls (only a subset of the chemicals were disolved in DMSO)
diffE.ctrl.u <- diffE[c("ctrl.Fluc","ctrl.Perm","ctrl.Dino","ctrl.Esfen","ctrl.Lind","ctrl.PCP")]

# filter differentially expressed loci: FDR <= 0.1
# c1<- list()
# for (i in 1:length(names(diffE.ctrl.u))){
#   c1[[i]] <- diffE.ctrl.u[[i]][which(diffE.ctrl.u[[i]][,"adj.P.Val"] <= 0.1),"logFC"]
# }
# names(c1) <- names(diffE.ctrl.u)
# summary(c1)

# filter differentially expressed loci: FDR <= 0.1 and no fold change threshold
# store the locusID and the fold change (log2)
diffEsig.FC.c<- list()
for (i in 1:length(names(diffE.ctrl.u))){
  diffEsig.FC.c[[i]] <- diffE.ctrl.u[[i]][which(diffE.ctrl.u[[i]][,"adj.P.Val"] <= 0.1 & (diffE.ctrl.u[[i]][,"logFC"] > 0 | diffE.ctrl.u[[i]][,'logFC'] < 0)),c("locusID","logFC")]
}
names(diffEsig.FC.c) <- names(diffE.ctrl.u)
summary(diffEsig.FC.c)

# merge diff. expressed for all controls in one matrix
for (n in names(diffEsig.FC.c)){
   colnames(diffEsig.FC.c[[n]])[2]<-paste(n,'logFC',sep='_')
 }

summary.diffE.c <- Reduce(function(x,y) merge(x,y, all=T,by.x='locusID',by.y='locusID'),diffEsig.FC.c, accumulate=F)
summary.diffE.c[is.na(summary.diffE.c)] <- 0
rownames(summary.diffE.c) <- summary.diffE.c$locusID
summary.diffE.c$locusID <- NULL

# colnames(summary.diffE.c) <- gsub("_logFC","",colnames(summary.diffE.c))

# plot heatmap of expression levels for the diff. expressed loci in all controls
breaks.all=c(-3,-2.5,-1.5,-1.35,-0.1,0.1,1.35,1.5,2.5,3)
col1=colorRampPalette(c("green4","green","lightgray","violet","purple"))

# dev.off()
par(oma=c(8,4,4,2))

library(gplots)
heatmap.2(as.matrix(summary.diffE.c),labR=NA,cexCol=1,col=col1,tracecol="lightgoldenrod",breaks=breaks.all)


