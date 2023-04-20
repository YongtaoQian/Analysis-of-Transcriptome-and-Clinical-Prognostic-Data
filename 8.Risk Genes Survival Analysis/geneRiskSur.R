install.packages('survminer')
library(survival)
library(survminer)

args = commandArgs(trailingOnly=TRUE)



rt=read.table(args[1],header=T,sep="\t")


diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}

surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=T,
		           pval=pValue,
		           pval.size=5,
		           risk.table=TRUE,
		           legend.labs=c("High risk", "Low risk"),
		           legend.title="Risk",
		           xlab="Time(years)",
		           break.time.by = 1,
		           risk.table.title="",
		           palette=c("red", "blue"),
		           risk.table.height=.25)
pdf(file="survival.pdf",onefile = FALSE,width = 12.5,height =9)
print(surPlot)
dev.off()




