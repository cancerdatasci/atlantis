
# $Id$

.onLoad <- function(lib, pkg) {
    GCtorture <<- FALSE
    .Call("party_init", PACKAGE = "atlantisPartyMod")
    return(TRUE)
}
