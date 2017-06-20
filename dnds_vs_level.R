nsc <- c("#FF5300", "#2678EB")
cor_meth<-"p"
mpc <- c(16,17)


pdf("pdf/dnds_lvl.sum.pdf", 7, 7, encoding="MacRoman")
x<-read.table("dnds_vs_level/red.dnds_vs_lvl.txt", header = T)
plot(x$dnds,x$nlvl, type = "p", pch=mpc[1], col=nsc[1],
     xlab = "dN/dS",
     xlim = c(0, 0.8),
     ylim = c(0, 30000),
     ylab = expression(paste("Cumulative editing level (%)"))
)
points(x$dnds,x$slvl, type = "p", pch=mpc[2], col=nsc[2])

cor.test(x$dnds,x$nlvl, method = cor_meth, alternative = "t")
cor.test(x$dnds,x$slvl, method = cor_meth, alternative = "t")
abline(lm(x$nlvl~x$dnds), col=nsc[1], lwd=2, lty=1)
abline(lm(x$slvl~x$dnds), col=nsc[2], lwd=2, lty=1)

text(0.4, 20000, expression(italic(r) == -0.66), col = nsc[1], cex = 0.8)
text(0.4, 18000, expression(italic(P) == "1.5" %*% 10^-3), col = nsc[1], cex = 0.8)
text(0.2, 10000, expression(italic(r) == 0.33), col = nsc[2], cex = 0.8)
text(0.2,  8000, expression(italic(P) == "0.16"), col = nsc[2], cex = 0.8)

legend("topright",
       c("Nonsyn", "Syn"),
       pch = mpc,
       cex = 0.8,
       col = nsc
)
dev.off()
