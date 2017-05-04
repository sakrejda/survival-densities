.libPaths("R/library")
options(repos=c(CRAN="https://cran.rstudio.com"))
#install.packages("packrat")
install.packages("drat")

args <- commandArgs(trailingOnly=TRUE)
source_dir <- args[1]

pkgs <- readLines(con=file.path(source_dir, "R-libraries.txt"))
pkgs <- pkgs[pkgs != '']
drat::addRepo("sakrejda")
for (pkg in pkgs) {
  if (!require(basename(pkg), character.only=TRUE)) {
    o <- try(install.packages(basename(pkg), type='source'))
    if (class(o) == 'try-error')
      try(install.packages(devtools::install_github(pkg)))
  }
  if (!require(basename(pkg), character.only=TRUE))
    warning(paste0("Package ", pkg, " failed to load.")) 
}

if (args[2] == 'update')
  update.packages(ask=FALSE, checkBuilt=FALSE)
#install.packages('rgdal', type = "source", configure.args=c('--with-proj-include=/usr/local/include','--with-proj-lib=/usr/local/lib'))
#packrat::set_opts(ignored.packages = c("waitup"))
#packrat::init()
#packrat::snapshot()


