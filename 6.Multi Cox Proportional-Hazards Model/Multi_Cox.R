library(survival)
library(survminer)
library(glmnet)
args = commandArgs(trailingOnly=TRUE)
coxPfilter=0.05  

rt=read.table(args[1],header=T,sep="\t",check.names=F,row.names=1)
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")    
multiCoxSum=summary(multiCox)


outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
write.table(outTab,file="multi.Cox.txt",sep="\t",row.names=F,quote=F)


riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
riskOut=cbind(rt[,outCol],riskScore,risk)
riskOut=cbind(id=rownames(riskOut),riskOut)
write.table(riskOut,file="geneRisk.txt",sep="\t",quote=F,row.names=F)


pdf(file="multi.forest.pdf",width = 10,height = 6,onefile = FALSE)
ggforest(multiCox,
         main = "Hazard ratio",
         cpositions = c(0.02,0.22, 0.4), 
         fontsize = 0.7, 
         refLabel = "reference", 
         noDigits = 2)
dev.off()
