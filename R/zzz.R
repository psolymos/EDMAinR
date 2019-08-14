.onAttach <- function(libname, pkgname){
    ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields=c("Version", "Date"))
    packageStartupMessage(paste(pkgname, ver[1], "\t", ver[2]))
    options("pboptions" = list(
        type = if (interactive()) "timer" else "none",
        char = "[=-]",
        txt.width = 50,
        gui.width = 300,
        style = 6,
        initial = 0,
        title = "EDMA in progress...",
        label = "",
        nout = 100L,
        min_time = 2))
    invisible(NULL)
}

.onUnload <- function(libpath){
    invisible(NULL)
}
