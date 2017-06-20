library(RColorBrewer)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length/2, ...)
}
colSd <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=sd, na.rm=na.rm)

x<-read.table("level_group_ns/ns_ratio.txt", header = T)
t0<-t.test(x[2:11,1], mu = x[1,1])$p.value
t1<-t.test(x[2:11,2], mu = x[1,2])$p.value
t2<-t.test(x[2:11,3], mu = x[1,3])$p.value
t3<-t.test(x[2:11,4], mu = x[1,4])$p.value
t4<-t.test(x[2:11,5], mu = x[1,5])$p.value
t5<-t.test(x[2:11,6], mu = x[1,6])$p.value

tp<-c(t0, t1, t2, t3, t4, t5)
tpf<-p.adjust(tp, method = "BH")

pdf("pdf/NC.nsyn.eq_interval.pdf", width = 4.5, height = 4.5, encoding = "MacRoman")
par(mar=c(4.2,4,0,0))
red_lvl<-brewer.pal(6,"Blues")
red_lvl[1]<-"red"
barx<-barplot(as.vector(t(x[1,])), col = red_lvl,
              xpd=T,
              ylim=c(0,8.5),
              xlab="Editing levels (%)",
              ylab="Nonsynonymous editing / Synonymous editing")

xn<-c("0-100","0-20", "20-40", "40-60", "60-80", "80-100")
axis(1, barx, xn, cex.axis=0.85, tick = F)

cm<-colMeans(x[2:11,])
cs<-colSd(x[2:11,])
points(barx, cm, pch=1)
error.bar(barx, cm, cs)
text(barx, x[1,]+0.25, "****")
dev.off()

