error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length/1.5, ...)
}

red.colors<-c("red", "#FF5700", "#4077C2", "#A7DD28", "gray")
pdf("pdf/red_sam.a2a.nonsyn.pdf", width = 6, height = 3.5, encoding = "MacRoman")
par(mar=c(4,4,0,0), xpd=T)
a2a<-read.table("sample_a2a/a2a.red.txt", nrows = 15)
red_lvl<-red.colors[a2a$Color+1]
barx<-barplot(as.matrix(t(a2a$Freq)),
              beside = T, col = red_lvl, ylim = c(0, max(a2a)+800),
              xlab = "Amino acid substitutions", axes = F, space = c(0.8,0.3), cex.axis =0.8,
              ylab = "Number of events")
a2a_sam_mean<-read.table("sample_a2a/a2a.sam.mean.txt", nrows = 15)
points(barx, t(a2a_sam_mean), type = "p")
a2a_sam_sd<-read.table("sample_a2a/a2a.sam.sd.txt", nrows = 15)
error.bar(barx, t(a2a_sam_mean), t(a2a_sam_sd))
axis(2, cex.axis = 0.8)
text(barx[1]-0.3, -500, expression(K %->% E), xpd=T, cex=0.8, srt=60)
text(barx[2]-0.3, -500, expression(S %->% G), xpd=T, cex=0.8, srt=60)
text(barx[3]-0.3, -500, expression(N %->% D), xpd=T, cex=0.8, srt=60)
text(barx[4]-0.3, -500, expression(R %->% G), xpd=T, cex=0.8, srt=60)
text(barx[5]-0.3, -500, expression(Y %->% C), xpd=T, cex=0.8, srt=60)
text(barx[6]-0.3, -500, expression(I %->% V), xpd=T, cex=0.8, srt=60)
text(barx[7]-0.3, -500, expression(T %->% A), xpd=T, cex=0.8, srt=60)
text(barx[8]-0.3, -500, expression(M %->% V), xpd=T, cex=0.8, srt=60)
text(barx[9]-0.3, -500, expression(I %->% M), xpd=T, cex=0.8, srt=60)
text(barx[10]-0.3, -500, expression(K %->% R), xpd=T, cex=0.8, srt=60)
text(barx[11]-0.3, -500, expression(Q %->% R), xpd=T, cex=0.8, srt=60)
text(barx[12]-0.3, -500, expression(E %->% G), xpd=T, cex=0.8, srt=60)
text(barx[13]-0.3, -500, expression(H %->% R), xpd=T, cex=0.8, srt=60)
text(barx[14]-0.3, -500, expression(N %->% S), xpd=T, cex=0.8, srt=60)
text(barx[15]-0.3, -500, expression(D %->% G), xpd=T, cex=0.8, srt=60)

# stat test
e1<-read.table("sample_a2a/a2a.red.stat.txt", header = T, check.names=FALSE)
r1<-read.table("sample_a2a/a2a.sam.stat.txt", header = T, check.names=FALSE)
p<-e1
t<-e1
f<-e1
a2a<-names(e1)
for (i in 1:26) {
  n<-a2a[i]
  p[i]<-t.test(r1[[n]], mu = e1[[n]])$p.value
  t[i]<-t.test(r1[[n]], mu = e1[[n]])$statistic
}
fdr<-p.adjust(p, method = "BH")[1:15]

pvs<-c(" ", "****")[as.numeric(fdr<0.00001)+1]
pvs_col = c("red", "blue")[as.numeric(t>0)+1]
text(barx, -100, pvs, col = pvs_col)

legend("topright", pch = c(16,16,16),
       c("Very different", "Somewhat different", "Simialr"),
       cex = 0.8,
       box.lwd = NA,
       col = red.colors[2:4]
       )

dev.off()
