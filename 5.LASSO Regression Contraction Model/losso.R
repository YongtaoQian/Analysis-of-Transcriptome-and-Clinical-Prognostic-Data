install.packages('survminer')
install.packages("glmnet")
args = commandArgs(trailingOnly=TRUE)


library(survival)
library(survminer)
library(glmnet)
coxPfilter=0.05          

rt=read.table(args[1],header=T,sep="\t",check.names=F,row.names=1)     
rt$futime=rt$futime/365


outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
	if(sd(rt[,i])<0.001){next}

	cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]

	if(coxP<coxPfilter){
	        sigGenes=c(sigGenes,i)
			outTab=rbind(outTab,
			             cbind(id=i,
			             HR=coxSummary$conf.int[,"exp(coef)"],
			             HR.95L=coxSummary$conf.int[,"lower .95"],
			             HR.95H=coxSummary$conf.int[,"upper .95"],
			             pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
			             )
	}
}
write.table(outTab,file="uni.Cox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uni.SigExp.txt",sep="\t",row.names=F,quote=F)


rt=read.table("uni.SigExp.txt",header=T,sep="\t",row.names=1,check.names=F)       
rt$futime[rt$futime<=0]=0.003
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit <- glmnet(x, y, family = "cox", maxit = 1000)
pdf("lasso.lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()
cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("lasso.cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene=c("futime","fustat",lassoGene)
lassoSigExp=rt[,lassoGene]
lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
write.table(lassoSigExp,file="lasso.SigExp.txt",sep="\t",row.names=F,quote=F)





