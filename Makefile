

doc:
	Rscript --no-init-file -e "library(devtools) ; devtools::document()"

install: doc
	R CMD INSTALL --preclean --no-multiarch --with-keep.source .


pack: doc
	R CMD build .

check: doc pack
	tar -xzvf Ty_*.tar.gz
	R CMD check Ty

quick-install:
	R CMD INSTALL --no-multiarch --with-keep.source .

install-deps:
	Rscript --no-init-file -e \
			'install.packages(c("devtools", "roxygens2"), repos="https://cran.univ-paris1.fr")'



clean:
	rm -f .RData

	rm -f Ty_*.tar.gz
	rm -fr Ty
	rm -fr Ty.Rcheck
