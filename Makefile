PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: clean build

clean: 
	echo "Clean"

build: doc
	echo $(PKGNAME) $(PKGVERS)
	cd ..;\
	R CMD build $(PKGSRC)

doc:
	R --slave -e 'library(roxygen2); roxygenise()'
	R --slave -e 'library(devtools); build_manual()'
	#-git add --all man/*.Rd

check:
	R --slave -e 'library(devtools); check(error_on="error")'
