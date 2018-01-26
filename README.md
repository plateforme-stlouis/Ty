Toward a R package ?

```bash
    # Build Install
    R CMD INSTALL --no-multiarch --with-keep.source .

    # Clean Build Install
    R CMD INSTALL --preclean --no-multiarch --with-keep.source .

    # Check package
    Rscript -e 'library(devtools) ; devtools::check(document = FALSE)'

    # Build source
    Rscript -e 'library(devtools) ; devtools::build()'

    # Build binary
    Rscript -e "library(devtools) ; devtools::build(binary = TRUE, args = c('--preclean'))"
```
