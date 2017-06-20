nsc <- c("#FF5300", "#2678EB")
cor_meth<-"p"
mpc <- c(16,17)

pdf("pdf/expr_fnfs.pdf",7,7, encoding="MacRoman")
x<-read.table("expr_vs_fnfs/expr_fnfs.txt", header = T)
plot(x$expr,x$fn, type = "p", pch=mpc[1], col=nsc[1],
     xlab = "Expression (FPKM)",
     xlim = c(0, 60),
     ylim = c(0, 30),
     ylab = expression(paste("Density of editing sites ", "(\u2030)"))
)
points(x$expr,x$fs, type = "p", pch=mpc[2], col=nsc[2])

cor.test(x$expr,x$fn, method = cor_meth, alternative = "t")
cor.test(x$expr,x$fs, method = cor_meth, alternative = "t")
abline(lm(x$fn~x$expr), col=nsc[1], lwd=2, lty=1)
abline(lm(x$fs~x$expr), col=nsc[2], lwd=2, lty=1)

text(10, 25, expression(italic(r) == 0.40), col = nsc[1], cex = 0.8)
text(10, 23.5, expression(italic(P) == "0.08"), col = nsc[1], cex = 0.8)
text(40, 10, expression(italic(r) == 0.29), col = nsc[2], cex = 0.8)
text(40, 8.5, expression(italic(P) == "0.21"), col = nsc[2], cex = 0.8)

legend("topright",
       c("Nonsyn", "Syn"),
       pch = mpc,
       cex = 0.8,
       col = nsc
)
dev.off()
