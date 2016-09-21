#!/usr/bin/Rscript

library(ggplot2)
library(reshape)
library(getopt)

spec = matrix(c('verbose','v',0,"logical",'help','h',0,"logical",'llfile','l',1,"character",'ofile','o',1,"character"),byrow=TRUE,ncol=4)

opt=getopt(spec)

# if help was asked for print a friendly message 
# and exit with a non-zero error code 
if( !is.null(opt$help)) {
    cat(getopt(spec, usage=TRUE)); 
    q(status=1);
}

llfile <- opt$llfile

LL <- read.csv(llfile,header=TRUE)
LLR <- LL[LL$G >= 3,]
p <- ggplot(LLR, aes(x = G, y = Dev)) + geom_point() + geom_smooth() + theme_bw() 

p <- p + ylab("Mean posterior deviance") + xlab("Number of strains - G") + theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

pdf(opt$ofile)
plot(p)
dev.off()
