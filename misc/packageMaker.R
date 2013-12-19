#################################################################################
#################################################################################
# Erstellen des eiwild-Paketes

library(roxygen2)
library(devtools)



load_data(pkg="eiwild/")

document("../pkg/eiwild/")

dev_help("tuneVars")
dev_example("tuneVars")

build("../pkg/eiwild/")
build("../pkg/eiwild/", binary=TRUE)


# During development
dev_mode()
install("../pkg/eiwild/")
library(eiwild) # uses the development version



load_all("../pkg/eiwild/")
document("../pkg/eiwild/")

devtools:::check("../pkg/eiwild")
clean_dll("../pkg/eiwild/")
install("../pkg/eiwild")

show_news("../pkg/eiwild/", FALSE)




