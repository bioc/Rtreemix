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

.onAttach <- function(libname, pkgname) {
    msg <- sprintf(
        "Package '%s' is deprecated and will be removed from Bioconductor
         version %s", pkgname, "3.22")
    .Deprecated(msg=paste(strwrap(msg, exdent=2), collapse="\n"))
}
