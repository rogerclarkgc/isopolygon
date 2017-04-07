library(animation)
for(i in 1:ani.options("nmax")){
polydebug <- iso.polydebug(s_poly, mix, c_poly, its = 1)
windows()
V <- polydebug$V
hull_a <- polydebug$hull_a
m <- polydebug$m
plot(0, 0, type = "n", xlim = c(-8, 30), ylim = c(-2, 24))
points(mix[,1], mix[,2], col = "black", pch = 16)
points(V[hull_a, 1], V[hull_a, 2], col = "blue", pch = 4, cex = 1.5)
lines(V[hull_a, 1], V[hull_a, 2], col = "darkgreen", lwd = 1.5)
m_sample <- sample(62500, 1000, rep = F)
points(m$x[m_sample], m$y_f[m_sample], col = "orangered", cex = .5)
ani.pause()
savePlot(filename = i, type = "png")
}

