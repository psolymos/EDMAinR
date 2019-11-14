#remotes::install_github("psolymos/EDMAinR")
library(EDMAinR)
B <- 9

file1 <- system.file("extdata/crouzon/Crouzon_P0_Global_MUT.xyz",
    package="EDMAinR")
x1 <- read_xyz(file1)
x1

file2 <- system.file("extdata/crouzon/Crouzon_P0_Global_NON-MUT.xyz",
    package="EDMAinR")
x2 <- read_xyz(file2)
x2

dim(x1)
dimnames(x1)
landmark_names(x1)

x1[1:10, 2:3, 1:5]
subset(x1, 1:10)

str(as.matrix(x1))
str(as.data.frame(x1))
str(stack(x1))
str(as.array(x1))

fit <- edma_fit(x1[1:5,,])
fit
Meanform(fit)
SigmaKstar(fit)

as.dist(fit)
stack(as.dist(fit))

get_fm(fit)
get_fm(fit, sort=TRUE, decreasing=TRUE)
get_fm(fit, sort=TRUE, decreasing=FALSE)

pc <- get_pca(edma_fit(x1))
plot(pc, pch=3)
polygon(pc[chull(pc),], col="#ff000022", lty=2, border=2)

numerator <- edma_fit(x1, B=B)
denominator <- edma_fit(x2, B=B)
#numerator <- edma_fit(x1[1:25,,])
#denominator <- edma_fit(x2[1:25,,])
fd <- formdiff(numerator, denominator)
str(fd)

fdm <- edma_fdm(numerator, denominator, B)

head(get_fdm(fdm))
head(get_fdm(fdm, sort=TRUE, decreasing=TRUE))
head(get_fdm(fdm, sort=TRUE, decreasing=FALSE))

test(fdm)
plot(fdm, "global")
plot(fdm, "local_p")
plot(fdm, "local_ci")

pc2 <- get_pca(fdm)
plot(pc2)

if (interactive()) {

#library(plot3D)
#library(plot3Drgl)

#fit <- edma_fit(x1)
#mf <- data.frame(Meanform(fit))

#scatter3D(mf$X, mf$Y, mf$Z)

#plotrgl()

}




