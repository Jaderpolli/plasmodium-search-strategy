library(adehabitatHR)

mu1 = 1.5

lev <- simm.levy(date = 1:1000, mu = mu1, l0 = 0.00001, x0 = c(0, 0), typeII = TRUE)
pdf(file= paste("levy",num2str(mu1),".pdf"), height = 4.5, width = 4.5)
plot(lev, xlab = "x", ylab = "y", pch = 20, font.main = 1, main = paste("μ = ", num2str(mu1)), lmin = 0.00001", final = FALSE)
dev.off()