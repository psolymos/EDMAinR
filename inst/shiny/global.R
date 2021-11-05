library(EDMAinR)
library(nat)
library(shiny)
library(shinydashboard)
library(ggplot2)
library(plotly)
library(reactable)
library(pbapply)
pboptions(type="shiny", title="Calculating")
options(rgl.useNULL = TRUE)
library(rgl)

ver <- read.dcf(
    file = system.file("DESCRIPTION", package = "EDMAinR"),
    fields = "Version")
