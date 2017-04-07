##USER INPUT
sources <- read.table("Sources_example.csv",header=T,sep=",") #always put 13C(x) before 15N(y)
mixture <- read.table("Mixture_example.csv",header=T,sep=",") #some error is required for every value
TEF <-  read.table("TEF_example.csv", header=T,sep=",")
its <- 1500  #specify the number of iterations ("its")
min_C <- -50  #specifiy the dimensions and resolution for the mixing region figure #must specifiy according to your sources and mixtures data
max_C <- -20    #choose values outside the 95% mixing region
min_N <- -2  
max_N <- 10
res <- 250 #resolution of the mixing region figure; reducing this improves performance

##Now RUN the simulation 
step_C <- (max_C - min_C)/(res - 1)
step_N <- (max_N - min_N)/(res - 1)   
C_g <- seq(min_C,max_C,by=step_C) #values must be in ascending order
N_g <- seq(min_N,max_N,by=step_N) #values must be in ascending order
mgrid <- function(a,b) {  #create a grid of values to test for P-I-P
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}
m <- mgrid(C_g,N_g)
Par_values <- array(0, c(its,(nrow(sources)*4+3)))  #create files to store data
p <- array(0, c(its,(nrow(mixture))))
mix_reg <- array(0, c(res,res))
for (i in 1:its) {    #run loops to generate source isotopic signatures, for iterations = 'its'
  v <- array(0, c(nrow(sources),2))
  f <- array(0, c(nrow(TEF),2))
  for (j in 1:nrow(sources)) {
    v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  #generate values from norm. dist. for d13C
    v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  #generate values from norm. dist. for d15N
    f[j,1] <- rnorm(1, mean=TEF[j,1], sd=TEF[j,2])  #generate values from norm. dist. for d13C enrichment
    f[j,2] <- rnorm(1, mean=TEF[j,3], sd=TEF[j,4])  #generate values from norm. dist. for d15N enrichment
  }
  V <- v+f
  hull <- chull(V)  #create a 2D convex hull from the enriched sources, find the boundary of point
  hull_a <- append(hull,hull[1])  #closes the polygon
  P <- point.in.polygon(mixture[,1], mixture[,2], V[hull_a,1], V[hull_a,2]) #calculate P_I_P 
  P_n <- as.numeric(P)                        
  p[i,] <- P_n
  poly_a <- areapl(V[hull_a,])  #calculate polygon area, for evaluating the quantity of iterations
  m$y_f <- m$y[res:1,]  #flip y grid data to resemble axes (d13C=x, d15N=y)  why??
  m_r <- point.in.polygon(m$x, m$y_f, V[hull_a,1], V[hull_a,2]) #calculate P-I-P for the mixing region
  m_r_s <- matrix(m_r,nrow=res,byrow=F)  #convert vector into square matrix
  m_r_s[m_r_s > 1] <- 1  #point.in.polygon can return '2' or '3'
  mix_reg <- mix_reg + m_r_s
  vals <- c(v[,1],v[,2],f[,1],f[,2],0,0,0)  #concatenate values for this iteration
  Par_values[i,] <- vals  #store values
  Par_values[i,ncol(Par_values)-2] <- poly_a
  Par_values[i,ncol(Par_values)-1] <- i
  Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2]) #caculate the variance of ploygon area
  if (i %% 10 ==0) cat(paste("iteration", i, "\n")) 
}
##FIGURE 1: variance of polygon area during simulation
Iterations <- Par_values[,ncol(Par_values)-1]
Variance <- Par_values[,ncol(Par_values)]
plot(Iterations, Variance, type="n")  #plots 
lines(Iterations, Variance, lty=1, lwd=1.5, col="blue")
##FIGURE 2: proportion of iterations that each consumer was inside mixing polygon
p[p > 1] <- 1   #point.in.polygon can return '2' or '3'
Probabilities <- colSums(p)/its
print(Probabilities)
windows() 
barplot(Probabilities, xlab="Consumer", ylab="Probability consumer is within mixing polygon", 
        ylim=c(0,1), names.arg=seq(1,nrow(mixture),by=1))
##FIGURE 3: mixing region, consumers, average enriched source signatures
mix_reg <- mix_reg/its  #convert to 0-1 scale
mix_reg[mix_reg==0] <- NA #make the zeros white
mix_regt <- t(mix_reg[ncol(mix_reg):1,])  #transpose matrix for use with 'image'
windows()
image(C_g, N_g, mix_regt, col=colorRampPalette(c("blue", "light blue", "green", "light green", 
                                                 "yellow", "red"))(100), xlab="d13C", ylab="d15N", useRaster=TRUE)
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TEF
points(sources_TEF[,1], sources_TEF[,3], col="white", pch=4, lwd=2, cex=1.5)
points(mixture, pch=19, cex=1.3)
dev.copy2pdf(file="Mix_Region.pdf")
windows()  #create colour bar for figure 3
cust_color <- colorRampPalette(c("blue", "light blue", "green", "light green", "yellow", "red"))
z <- matrix(1:100, nrow=1)
x <- 1
y <- seq(0,1,len=100)
image(x,y,z,col=colorRampPalette(c("blue", "light blue", "green", "light green", "yellow", "red"))(100), 
      xaxt="n", xlab="", ylab="", useRaster=TRUE, bty="n", las=1)
##FIGURE 4: Same as Fig. 3, but black and white
windows()
plot(C_g, N_g, type = "n", xlab="d13C", ylab="d15N")
cont <- c(0.05, seq(0.1, 1, by=0.1))
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TEF
points(sources_TEF[,1], sources_TEF[,3], col="black", pch=4, lwd=2, cex=1.5)
points(mixture, pch=19, cex=1.3)
##FIGURE 5: Bi-plot with single 95% contour line
windows()
plot(C_g, N_g, type = "n", xlab="d13C", ylab="d15N")
cont <- c(0.05)
contour(C_g, N_g, mix_regt, levels=cont, add=TRUE, drawlabels=FALSE, lwd=1.9)
sources_TEF <- sources + TEF
points(sources_TEF[,1], sources_TEF[,3], col="black", pch=15, lwd=2, cex=1.5)
arrows(sources_TEF[,1]-sources_TEF[,2], sources_TEF[,3], #horizontal error bars (SD)
       sources_TEF[,1]+sources_TEF[,2], sources_TEF[,3], length=0, angle=90, code=3)
arrows(sources_TEF[,1], sources_TEF[,3]-sources_TEF[,4], #vertical error bars (SD)
       sources_TEF[,1], sources_TEF[,3]+sources_TEF[,4], length=0, angle=90, code=3)
points(mixture, pch=1, cex=1.3)
labels <- c("s1", "s2", "s3", "s4", "s5") #CHANE THIS FOR YOUR OWN DATA
text(sources_TEF[,1], sources_TEF[,3], labels=labels, pos=3)

##Write data to files
p_a <- rbind(p,Probabilities)
write.table(p_a, file = "Consumer_Probabilities.csv", sep = ",", row.names=FALSE)
col_names <- c(rep("d13C",nrow(sources)),rep("d15N",nrow(sources)),
               rep("13C_TEF",nrow(sources)),rep("15N_TEF",nrow(sources)),
               "Poly_Area","Iteration","Variance") 
col_nums <- c(rep(1:nrow(sources),4),0,0,0)
col_n <- paste(col_names, col_nums)
write.table(Par_values, file = "Parameter_Values.csv", sep = ",", row.names=FALSE, col.names=col_n)
