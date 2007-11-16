#.onLoad <- function(lib, pkg) {

#}

#.onAttach <- function(lib, pkg) {
#  ## some preprocessing
#  where <- match(paste("package:", pkg, sep=""), search())
#  groupGOTerms(where)
#}

.onUnload <- function( libpath ) {
   library.dynam.unload( "Rtreemix", libpath )
}
