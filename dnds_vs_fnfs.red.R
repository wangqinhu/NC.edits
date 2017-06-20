nsc <- c("#FF5300", "#2678EB")
cor_meth<-"p"
mpc <- c(16,17)

pdf("pdf/dnds_vs_fnfs.red.pdf", 7, 7, encoding="MacRoman")
x<-read.table("dnds_vs_fnfs/red.txt", header = T)
plot(x$dnds,x$fn, type = "p", pch=mpc[1], col=nsc[1],
     xlab = "dN/dS",
     xlim = c(0, 0.8),
     ylim = c(10,30),
     ylab = expression(paste("Density of editing sites ", "(\u2030)"))
)
points(x$dnds,x$fs, type = "p", pch=mpc[2], col=nsc[2])

cor.test(x$dnds,x$fn, method = cor_meth, alternative = "t")
cor.test(x$dnds,x$fs, method = cor_meth, alternative = "t")
abline(lm(x$fn~x$dnds), col=nsc[1], lwd=2, lty=1)
abline(lm(x$fs~x$dnds), col=nsc[2], lwd=2, lty=1)

text(0.4, 21, expression(italic(r) == -0.76), col = nsc[1], cex = 0.8)
text(0.4, 20, expression(italic(P) == "9.0" %*% 10^-5), col = nsc[1], cex = 0.8)
text(0.1, 14, expression(italic(r) == -0.58), col = nsc[2], cex = 0.8)
text(0.1, 13, expression(italic(P) == "8.0" %*% 10^-3), col = nsc[2], cex = 0.8)

legend("topright",
       c("Nonsyn", "Syn"),
       pch = mpc,
       cex = 0.8,
       col = nsc
)
dev.off()
