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


library(Morpho)

f <- "~/Dropbox/consulting/2020/edma/data/HlavUSNM143565LAND.txt"
Hlav <- ply2mesh(f, adnormals = TRUE, readnormals = FALSE, readcol = FALSE, silent = FALSE)

Hlav.mesh <- tps3d(Hlav, ref, lsmeans[,,"Hlar"])

shade3d(Hlav.mesh, col="green")

ref <- mshape(Y$coords)

avg.mesh1 <- tps3d(Hlav, F[,,"HlavM6_USNM143565"], ref)

shade3d(avg.mesh1, col="navy")

shade3d(Hlav.mesh, col="darkturquoise")


