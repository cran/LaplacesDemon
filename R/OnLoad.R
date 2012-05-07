.onLoad <- function(libname, pkgname) {
     version <- utils::packageDescription("LaplacesDemon")$Version
     packageStartupMessage("\nLaplacesDemon : Software for Bayesian Inference")
     packageStartupMessage("\nVersion: ",version)
     packageStartupMessage("\n\n``Probability theory is nothing but common sense reduced to")
     packageStartupMessage("calculation'' (Pierre-Simon Laplace, 1814).\n")     
     packageStartupMessage("Laplace's Demon is ready for you.\n")
     #library.dynam("LaplacesDemon", pkgname, libname)
}
