#remotes::install_github("psolymos/EDMAinR")

set.seed(123)
library(EDMAinR)
library(geomorph)

## this all return a K*D*n (from single or multiple files)
read.morphologika
readland.nts
readland.tps
readmulti.nts
readmulti.tps

## K*D
readland.fcsv


## K*D*n type with semilandmarks
readland.shapes

## this returns mesh3d
read.ply

## we need a function which puts a 3D array into an edma_data object

data(plethodon)
Y.gpa<-gpagen(plethodon$land)
#GPA-alignmentplot
plotAllSpecimens(Y.gpa$coords,links=plethodon$links)

plot(as.edma_data(plethodon$land), which=NULL, ask=NA, hull=TRUE,
     xlim=c(-15,10), ylim=c(-3,3))
