#' @useDynLib LAVA, .registration = TRUE

.onLoad = function(libname, pkgname) {
	packageStartupMessage(paste("Running LAVA version", packageVersion("LAVA")))
}
