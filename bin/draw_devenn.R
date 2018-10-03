library(VennDiagram)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE) 

if (length(args)!=9) {
  stop("usage : RScript --vanilla graphs.R <size 1> <size 2> <common size> <name1> <name2> <color 1> <color 2> <out> <title>", call.=FALSE)
}

a1=args[1]
a1=strtoi(a1)
a2=args[2] 
a2=strtoi(a2)
a3=args[3]
a3=strtoi(a3)

p=draw.pairwise.venn(a1,a2,a3,fill=c(args[6],args[7]),alpha=c(0.3,0.3),category=c(args[4],args[5]),cat.pos=c(190,170),cat.just=list(c(0,0),c(0,0)),cat.dist=c(0.025,0.025),cat.cex=c(1,1),ext.dist=c(-0.06,0.1),ext.length=c(0.9,0.95),ext.pos=c(0,20),margin=c(0.2))
title=textGrob(args[9], gp=gpar(fontsize=18))

pdf(args[8],onefile=FALSE) 
grid.arrange(grid.arrange(gTree(children=p), top=title)) 
dev.off() 
