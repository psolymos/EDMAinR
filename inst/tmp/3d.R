library(EDMAinR)

file1 <- system.file("extdata/crouzon/Crouzon_P0_Global_MUT.xyz",
    package="EDMAinR")
x1 <- read_xyz(file1)

file2 <- system.file("extdata/crouzon/Crouzon_P0_Global_NON-MUT.xyz",
    package="EDMAinR")
x2 <- read_xyz(file2)

numerator <- edma_fit(x1)
denominator <- edma_fit(x2)
fdm <- edma_fdm(numerator, denominator, B=0)

.formdiff=EDMAinR:::.formdiff
.get_data=EDMAinR:::.get_data
.compare_data=EDMAinR:::.compare_data

system.time(EDMAinR:::.Ttest_fit(numerator, denominator, B=100, ncores=1))
system.time(EDMAinR:::.Ttest_fit(numerator, denominator, B=100, ncores=2))

plot_3d(fdm)
.plot_d_dm(fdm, all=TRUE, midpoints=TRUE, cex=0.5)
.plot_d_dm(fdm, all=TRUE, midpoints=TRUE, cex=0.5, breaks=c(0, 0.9, 1))

.plot_d_dm(fdm, all=TRUE, midpoints=TRUE, d3=FALSE)
.plot_d_dm(fdm, all=FALSE, midpoints=TRUE, d3=FALSE, alpha=0.1)

.plot_d_dm(fdm, all=TRUE, midpoints=F, d3=FALSE, alpha=0.5)


