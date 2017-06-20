#
# random<=>H0<=>P(AAA/ANA|edited)=P(AAA/ANA|random)
# heterozygote advantage<=>H1<=>P(AAA/ANA|edited)>P(AAA/ANA|random)
#
M<-as.table(rbind(c(17810, 16888), c(232, 299)))
dimnames(M) <- list(types = c("AAA", "non-AAA"),
                    samples = c("edited","random"))
x<-chisq.test(M)
x$p.value
f<-fisher.test(M)
f$p.value

# nonsyn
M<-as.table(rbind(c(15921, 14708), c(135, 180)))
dimnames(M) <- list(types = c("AAA", "non-AAA"),
                    samples = c("edited","random"))
x<-chisq.test(M)
x$p.value
f<-fisher.test(M)
f$p.value

# syn
M<-as.table(rbind(c(1889, 2180), c(97, 119)))
dimnames(M) <- list(types = c("AAA", "non-AAA"),
                    samples = c("edited","random"))
x<-chisq.test(M)
x$p.value
f<-fisher.test(M)
f$p.value

# plot

nsc <- c("#FD6C9E", "#A6E7FF")

non<-(1-c(0.991591928, 0.987909726))*100
syn<-(1-c(0.951158107, 0.948238365))*100
names(non)<-c("Edited", "Random")
names(syn)<-c("Edited", "Random")

pdf("pdf/ANA_ratio.pdf", height = 3, width = 4)
layout(matrix(c(1,2),1,2,byrow = TRUE))
par(mar=c(2.25,4.5,2.8,0),xpd=T, cex=0.8)

x<-barplot(non, ylim=c(0,1.2), ylab="non-AAA %", col=nsc[1])
y <- 1.31
offset <- 0.08
lines(x[1:2],c(y, y))
lines(x[c(1,1)],c(y, y-offset))
lines(x[c(2,2)],c(y, y-offset))
text(x[1]+((x[2]-x[1])/2),y+offset, expression(italic(P) == 1.5 %*% 10^-3))

x<-barplot(syn, ylim=c(0,5), ylab="non-AAA %", col=nsc[2])
y <- 5.5
offset <- 0.3
lines(x[1:2],c(y, y))
lines(x[c(1,1)],c(y, y-offset))
lines(x[c(2,2)],c(y, y-offset))
text(x[1]+((x[2]-x[1])/2),y+offset, expression(italic(P) == "0.71"))

dev.off()