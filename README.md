Toward a R package ?

Build depends on:
 - devtools
 - roxygen2


```bash
    # Build Install
    R CMD INSTALL --no-multiarch --with-keep.source .

    # Clean Build Install
    R CMD INSTALL --preclean --no-multiarch --with-keep.source .

    # Check package
    Rscript -e 'library(devtools) ; devtools::check(document = FALSE)'

    # Check package and documentation
    Rscript -e 'library(devtools) ; devtools::check(document = TRUE)'

    # Build source
    Rscript -e 'library(devtools) ; devtools::build()'

    # Build binary
    Rscript -e "library(devtools) ; devtools::build(binary = TRUE, args = c('--preclean'))"

    # Build documentation with Roxygen2
    Rscript -e "library(devtools) ; devtools::document()"
```
