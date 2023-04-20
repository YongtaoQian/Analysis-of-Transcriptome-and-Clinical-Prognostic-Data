
library(survival)
args = commandArgs(trailingOnly=TRUE)


library(survival)
pFilter=0.05                                                              
rt=read.table(args[1],header=T,sep="\t",check.names=F,row.names=1)    
rt$futime=rt$futime                                                   

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
	   if(sd(rt[,gene])<0.01){
	      next}

	   a=rt[,gene]<=median(rt[,gene])
	   diff=survdiff(Surv(futime, fustat) ~a,data = rt)
	   pValue=1-pchisq(diff$chisq,df=1)
	   fit=survfit(Surv(futime, fustat) ~ a, data = rt)

	   cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
	   coxSummary = summary(cox)
	   coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	   if((pValue<pFilter) & (coxP<pFilter)){
	         sigGenes=c(sigGenes,gene)
	       	 outTab=rbind(outTab,
	                      cbind(gene=gene,
	                            KM=pValue,
	                            B=coxSummary$coefficients[,"coef"],
	                            SE=coxSummary$coefficients[,"se(coef)"],
	                            HR=coxSummary$conf.int[,"exp(coef)"],
	                            HR.95L=coxSummary$conf.int[,"lower .95"],
	                            HR.95H=coxSummary$conf.int[,"upper .95"],
			                    pvalue=coxP) )
	  }
}
write.table(outTab,file="uniCox.xls",sep="\t",row.names=F,quote=F)   
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)
