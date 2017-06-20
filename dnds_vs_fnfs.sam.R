require(plotrix)
nsc <- c("#FF5300", "#2678EB")
cor_meth<-"p"
mpc <- c(16,17)

pdf("pdf/dnds_vs_fnfs.sam.pdf", 7, 7, encoding="MacRoman")

dnds<-read.table("dnds_vs_fnfs/sam.dnds.txt")$V1
fn<-read.table("dnds_vs_fnfs/sam.fn.txt")
fs<-read.table("dnds_vs_fnfs/sam.fs.txt")

fn_median<-rep(NA,20)
fn_05ci<-rep(NA,20)
fn_95ci<-rep(NA,20)
for (i in 1:20) {
  fn_median[i]<-median(fn[[i]])
  fn_05ci[i]<-quantile(fn[[i]], c(0.05))
  fn_95ci[i]<-quantile(fn[[i]], c(0.95))
}

fs_median<-rep(NA,20)
fs_05ci<-rep(NA,20)
fs_95ci<-rep(NA,20)
for (i in 1:20) {
  fs_median[i]<-median(fs[[i]])
  fs_05ci[i]<-quantile(fs[[i]], c(0.05))
  fs_95ci[i]<-quantile(fs[[i]], c(0.95))
}

plotCI(dnds,fn_median, li=fn_05ci, ui=fn_95ci, slty = 3, 
       xlim=c(0,0.8),
       ylim=c(10,30),
       xlab = ("dN/dS"),
       ylab = expression(paste("Density of random editable sites ", "(\u2030)")),
       col = nsc[1],
       pch = mpc[1]
       )
plotCI(dnds,fs_median, li=fs_05ci, ui=fs_95ci, add=T, slty = 3,
       col = nsc[2],
       pch = mpc[2]
       )

cor.test(dnds,fn_median, method = cor_meth, alternative = "t")
cor.test(dnds,fs_median, method = cor_meth, alternative = "t")
abline(lm(fn_median~dnds), col=nsc[1], lwd=2, lty=1)
abline(lm(fs_median~dnds), col=nsc[2], lwd=2, lty=1)

text(0.4, 25, expression(italic(r) == 0.85), col = nsc[1], cex = 0.8)
text(0.4, 24, expression(italic(P) == "2.5" %*% 10^-6), col = nsc[1], cex = 0.8)
text(0.6, 17, expression(italic(r) == -0.02), col = nsc[2], cex = 0.8)
text(0.6, 16, expression(italic(P) == "0.93"), col = nsc[2], cex = 0.8)

legend("topright",
       c("Nonsyn", "Syn"),
       pch = mpc,
       cex = 0.8,
       col = nsc
)

dev.off()