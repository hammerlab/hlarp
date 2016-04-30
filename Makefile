
PACKAGES=nonstd re re.posix re.glob cmdliner
PACKAGES_INSTALL=$(PACKAGES)

SOURCE_DIRS=/util /unc /stats /cls /rgr /uns

#install uninstall setup

.PHONY: default clean build setup

default: build

# This should be called something else.
setup:
	opam install $(PACKAGES_INSTALL)

build:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/app -I src/lib hlarp_cli.native

clean:
	ocamlbuild -clean

#install:
#	ocamlfind install META \

#uninstall:
#	ocamlfind remove oml
