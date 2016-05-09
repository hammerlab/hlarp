PACKAGES=nonstd re cmdliner
PACKAGES_INSTALL=$(PACKAGES)

.PHONY: default clean build deps install uninstall

default: build

deps:
	opam install $(PACKAGES_INSTALL)

build:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) \
           -I src/app -I src/lib hlarp_cli.native \
           hlarp.cma hlarp.cmxs hlarp.cmxa
	cp hlarp_cli.native hlarp

clean:
	ocamlbuild -clean
	-rm ./hlarp


install:
	ocamlfind install hlarp META\
	    _build/src/lib/hlarp.a\
	    _build/src/lib/hlarp.o\
	    _build/src/lib/hlarp.cma\
	    _build/src/lib/hlarp.cmi\
	    _build/src/lib/hlarp.cmo\
	    _build/src/lib/hlarp.cmx\
	    _build/src/lib/hlarp.cmxa\
            _build/src/lib/hlarp.cmxs
uninstall:
	ocamlfind remove hlarp
