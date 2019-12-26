## stacked growth matrix
## inputs are edma_fit objects
edma_gm <- function (a1, a2, ...) {
    .compare_objects(a1, a2, "a1 vs a2:")
    out <- edma_fdm(numerator=a2, denominator=a1, ...)
    class(out) <- c("edma_gm", class(out))
    out$call <- match.call()
    out
}
print.edma_gm <- function(x, ...) {
    .print_edma_fdm(x, "EDMA growth matrix", ...)
}

## stacked growth difference matrix
## inputs are edma_fit objects
## NOTE: bootstrap is taken to be Ba x Bb + 1
## using all combinations from GM objects
edma_gdm <- function (a1, a2, b1, b2, ...) {
    .compare_objects(a1, a2, "a1 vs a2:")
    .compare_objects(b1, b2, "b1 vs b2:")
    .compare_objects(a1, b1, "a1 vs b1:")
    .compare_objects(a2, b2, "a2 vs b2:")
    # ref is a1 an b1
    gma <- edma_fdm(numerator=a2, denominator=a1, ...)
    gmb <- edma_fdm(numerator=b2, denominator=b1, ...)
    B <- gma$B * gmb$B
    gdm <- gma$dm
    gdm$dist <- gmb$dm$dist / gma$dm$dist
    ii <- expand.grid(a=seq_len(gma$B), b=seq_len(gmb$B))
    b <- cbind(gdm$dist,
            gmb$boot[,ii$b+1,drop=FALSE] / gma$boot[,ii$a+1,drop=FALSE])
    out <- list(
        call=match.call(),
        a1=a1, a2=a2, b1=b1, b2=b2,
        B=B,
        ref_denom=gma$ref_denom,
        dm=gdm,
        boot=b)
    attr(out$boot, "Tval") <- apply(out$boot, 2, max) / apply(out$boot, 2, min)
    class(out) <- c("edma_gdm", "edma_dm", class(out))
    out
}

print.edma_gdm <- function(x, ...) {
    .print_edma_fdm(x, "EDMA growth difference matrix", ...)
}

T_test.edma_gm <- function (object, ...)
    .T_test(object, DNAME="growth matrix", ...)
T_test.edma_gdm <- function (object, ...)
    .T_test(object, DNAME="growth difference matrix", ...)

get_gm <- function (object, ...) UseMethod("get_gm")
get_gm.edma_gm <- get_fdm.edma_fdm

get_gdm <- function (object, ...) UseMethod("get_gdm")

get_gdm.edma_gdm <- get_fdm.edma_fdm

