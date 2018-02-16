Toward a R package ?

Build depends on:
 - devtools
 - roxygen2

```bash
    # Dependencies for building the package
    Rscript --no-init-file -e \
        'install.packages(c("devtools", "roxygens2"), repos="https://cran.univ-paris1.fr")'
```

Package depends on:
 - readr
 - doBy
(should be manage by the Install process from DECRIPTION file)

```bash
    # Dependencies of the packahe
    Rscript --no-init-file -e \
        'install.packages(c("readr", "doBy"), repos="https://cran.univ-paris1.fr")'

```


Some command lines for the build process
```bash
    # Build Install
    R CMD INSTALL --no-multiarch --with-keep.source .

    # Clean Build Install
    R CMD INSTALL --preclean --no-multiarch --with-keep.source .

    # Check package
    Rscript --no-init-file -e 'library(devtools) ; devtools::check(document = FALSE)'

    # Check package and documentation
    Rscript --no-init-file -e 'library(devtools) ; devtools::check(document = TRUE)'

    # Build source
    Rscript --no-init-file -e 'library(devtools) ; devtools::build()'

    # Build binary
    Rscript --no-init-file -e "library(devtools) ; devtools::build(binary = TRUE, args = c('--preclean'))"

    # Build documentation with Roxygen2
    Rscript --no-init-file -e "library(devtools) ; devtools::document()"
```
