mat1 <- matrix(sample(1:30,10),nrow=5,ncol=2)
ses <- sample(1:3,10,replace=T)
vect <- seq(20,100,20)


# I added the vect vector as 'x' values.
# No more need to define the axes
# Just supressed the box around the plot with bty='n'
# to fit the look of your original plot
matplot(x=vect, y=mat1,
        pch=c('x','o'), type = "b", lwd = 2, lty = c(1,2),
        col = c("green","black"),
        main = "Graph 1", xlab = "Numbers 1", ylab = "Numbers 2",
        cex.main = 1.8, cex=2, cex.lab=1.5, cex.axis = 1.6, bty='n')

# I changed ses2 to ses and matrices to vectors
plotCI(x=rep(vect,2), y= as.vector(mat1),
       col=rep(c("green","black"),each=nrow(mat1)),
       add=T) 
