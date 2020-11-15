#remotes::install_github("psolymos/EDMAinR")
library(EDMAinR)
B <- 3

## --- data ---

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
landmarks(x1)
specimens(x1)
dimensions(x1)

x1[1:10, 2:3, 1:5]
subset(x1, 1:10)

str(as.matrix(x1))
str(as.data.frame(x1))
str(stack(x1))
str(as.array(x1))

plot(x1, which=1)
plot_2d(x1, which=2)
plot_ord(x1)
plot_clust(x1)

## simulate data
K <- 3 # number of landmarks
D <- 2 # dimension, 2 or 3
sig <- 0.75
rho <- 0
SigmaK <- sig^2*diag(1, K, K) + sig^2*rho*(1-diag(1, K, K))
M <- matrix(c(0,1,0,0,0,1), 3, 2)
M[,1] <- M[,1] - mean(M[,1])
M[,2] <- M[,2] - mean(M[,2])
M <- 10*M
edma_simulate_data(10, M, SigmaK)

## --- nonparametric fit ---

fit <- edma_fit(x1, B=B)
fit
str(Meanform(fit))
str(SigmaKstar(fit))

plot_ord(fit)
plot_clust(fit)
plot(fit)
plot_2d(fit)
if (interactive()) plot_3d(fit)

## --- form matrix ---

str(as.dist(fit))
str(stack(as.dist(fit)))

head(confint(fit))
head(get_fm(fit))
head(get_fm(fit, sort=TRUE, decreasing=TRUE))
head(get_fm(fit, sort=TRUE, decreasing=FALSE))

## --- parametric fit ---
if (FALSE) {
read_pattern(system.file("extdata/example.csv", package="EDMAinR"))
read_pattern(system.file("extdata/example.xlsx", package="EDMAinR"))

p <- matrix(c(
    "a", NA, NA, NA,
    NA, "a", NA, NA,
    NA,  NA, "b", NA,
    NA,  NA, NA, "b"
), 4, 4, byrow=TRUE)
parm <- c(a=0.25, b=0.35)
print_tb(p)
plot_tb(p)
SigmaK <- make_Sigma(parm, p)
dimnames(p) <- dimnames(SigmaK)

M <- matrix(c(-2.5, 7.5, -2.5, -2.5,
    -7.5, 2.5, 2.5, 4.5),
    4, 2)

sim <- edma_simulate_data(n=500, M, SigmaK)

fit <- edma_fit(sim, B=5)
o <- SigmaK_fit(fit, p)
o

cbind(true=parm, est=o$results$par[names(parm)])
SigmaK(o)
SigmaKfull(o)

sensitivity(o)

plot_ord(o)
plot_clust(o)
}

## --- form difference matrix ---

numerator <- edma_fit(x1, B=B)
denominator <- edma_fit(x2, B=B)

fdm0 <- edma_fdm(numerator, denominator, B=B, mix=TRUE)
fdm <- edma_fdm(numerator, denominator, B=B)
fdm2 <- edma_fdm(numerator, denominator, B=B, ref_denom = FALSE)
fdm
fdm2

head(get_fdm(fdm))
head(get_fdm(fdm, sort=TRUE, decreasing=TRUE))
head(get_fdm(fdm, sort=TRUE, decreasing=FALSE))

T_test(fdm)
T_test(fdm2)

head(confint(fdm))

head(infl <- get_influence(fdm))
plot(infl)

plot_ord(fdm)
plot_clust(fdm)
plot_Ttest(fdm)
plot_ci(fdm)
plot_2d(fdm)
if (interactive()) plot_3d(fdm)

## --- growth matrix ---

file_a1 <- system.file("extdata/purplebook/cebusage1.xyz",
    package="EDMAinR")
a1 <- read_xyz(file_a1)
file_a2 <- system.file("extdata/purplebook/cebusage6.xyz",
    package="EDMAinR")
a2 <- read_xyz(file_a2)
a1
a2

fit_a1 <- edma_fit(a1, B=B)
fit_a2 <- edma_fit(a2, B=B)

gm <- edma_gm(fit_a1, fit_a2, B=B)
gm
T_test(gm)
head(confint(gm))
head(get_gm(gm))
head(get_gm(gm, sort=TRUE, decreasing=TRUE))
head(get_gm(gm, sort=TRUE, decreasing=FALSE))

plot_ord(gm)
plot_clust(gm)
plot_Ttest(gm)
plot_ci(gm)
plot_2d(gm)
if (interactive()) plot_3d(gm)

## --- growth difference matrix ---

gdm <- edma_gdm(a1=fit_a1, a2=fit_a2, b1=fit_a1, b2=fit_a2, B=B)
gdm
T_test(gdm)
head(confint(gdm))
head(get_gdm(gdm))
head(get_gdm(gdm, sort=TRUE, decreasing=TRUE))
head(get_gdm(gdm, sort=TRUE, decreasing=FALSE))

plot_ord(gdm)
plot_clust(gdm)
plot_Ttest(gdm)
plot_ci(gdm)
#plot_2d(gdm) # need real data
#if (interactive()) plot_3d(gdm)

