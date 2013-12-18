#################################################################################
#################################################################################
# Erstellen des eiwild-Paketes

library(roxygen2)
library(devtools)


load_all("eiwild/")
load_data(pkg="eiwild/")

document("eiwild/")

dev_help("tuneVars")
dev_example("tuneVars")

build("eiwild")
devtools::build("eiwild",binary = TRUE)
