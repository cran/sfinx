# Loading options
.onLoad <- function(libname, pkgname) {
  op <- options()
  op.sfinx <- list(
    sfinx.name = "Kevin Titeca",
    sfinx.desc.author = '"Kevin Titeca <sfinxinteractomics@gmail.com> [aut, cre]"',
    sfinx.desc.license = "Apache 2.0",
    sfinx.desc.suggests = NULL,
    sfinx.desc = list()
  )
  toset <- !(names(op.sfinx) %in% names(op))
  if(any(toset)) options(op.sfinx[toset])

  invisible()
}
