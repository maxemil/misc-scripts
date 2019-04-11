tab = read.csv('report.0.tab', sep='\t', header=TRUE)

ident.97 = data.frame(tab$X1.97[which(tab$X1.97 > 0)])
rownames(ident.97) = tab$Environment_source[which(tab$X1.97 > 0)]
colnames(ident.97) = 'ident97'

ident.97.o = data.frame(ident.97[with(ident.97, order(ident97)), ])
rownames(ident.97.o) = rownames(ident.97)[order(ident.97$ident97)]
colnames(ident.97.o) = 'ident97'

ident.99 = data.frame(tab$X1.99[which(tab$X1.99 > 0)])
rownames(ident.99) = tab$Environment_source[which(tab$X1.99 > 0)]
colnames(ident.99) = 'ident99'

ident.99.o = data.frame(ident.99[with(ident.99, order(ident99)), ])
rownames(ident.99.o) = rownames(ident.99)[order(ident.99$ident99)]
colnames(ident.99.o) = 'ident99'

pdf('all_abundances.pdf', width=7, height=9.8994)
  par(mfrow = c(2, 1), mar=c(5,15,4,2))
  par(las=2)
  barplot(ident.97.o$ident97, names.arg=rownames(ident.97.o), horiz=TRUE, main="97%")
  barplot(ident.99.o$ident99, names.arg=rownames(ident.99.o), horiz=TRUE, main="99%")
dev.off()
