.onAttach <- function(libname, pkgname)
  {
  packageStartupMessage("Package: ",paste(pkgname),"  The package implements the approach proposed in Baghaei and Hohensinn (submitted).")
}