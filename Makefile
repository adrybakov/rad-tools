.DEFAULT_GOAL := help

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
# help:
# 	@$(SPHINXBUILD) -M help "docs/$(SOURCEDIR)" "docs/$(BUILDDIR)" $(SPHINXOPTS) $(O)

# $(O) is meant as a shortcut for $(SPHINXOPTS).
html: 
	@$(SPHINXBUILD) -M html "docs/$(SOURCEDIR)" "docs/$(BUILDDIR)" $(SPHINXOPTS) $(O)

html-examples: 
	@make examples
	@$(SPHINXBUILD) -M html "docs/$(SOURCEDIR)" "docs/$(BUILDDIR)" $(SPHINXOPTS) $(O)

doctest: 
	@$(SPHINXBUILD) -b doctest "docs/$(SOURCEDIR)" "docs/$(BUILDDIR)" $(SPHINXOPTS) $(O)

clean:
	-@rm -r docs/build
	-@rm -r docs/source/api/generated
	-@rm -r docs/source/api/crystal/generated
	-@rm -r docs/source/api/exchange/generated
	-@rm -r docs/source/api/spinham/generated
	-@rm -r docs/source/api/magnons/generated
	-@rm -r docs/source/api/_autosummary
	-@rm -r rad_tools.egg-info
	-@rm -r build
	-@rm -r dist
	-@rm -r .env3.11/lib/python3.11/site-packages/radtools
	-@rm -r .env3.11/lib/python3.11/site-packages/rad_tools*
	-@rm -r .env3.11/bin/rad-*
	-@rm -r .env3.10/lib/python3.10/site-packages/radtools
	-@rm -r .env3.10/lib/python3.10/site-packages/rad_tools*
	-@rm -r .env3.10/bin/rad-*
	-@rm -r .env3.9/lib/python3.9/site-packages/radtools
	-@rm -r .env3.9/lib/python3.9/site-packages/rad_tools*
	-@rm -r .env3.9/bin/rad-*
	-@rm -r .env3.8/lib/python3.8/site-packages/radtools
	-@rm -r .env3.8/lib/python3.8/site-packages/rad_tools*
	-@rm -r .env3.8/bin/rad-*

install:
	@python3 -m pip install .

test: 
	@pytest -s

test-all: clean install test bravais-pictures examples html doctest
	@echo "Done"


.ONESHELL:
pip: prepare-release
	@read -p "Press Enter to publish to PyPi"
	-@rm -r dist
	-@rm -r build
	-@rm -r radtools.egg-info
	@python3 setup.py sdist bdist_wheel 
	@python3 -m twine upload --repository pypi dist/* --verbose
	@git tag -a "$(VERSION)" -m "Version $(VERSION)"


help:
	@echo "\x1b[31m"
	@echo "Please specify what do you want to do!"
	@echo "\x1b[0m"
	@echo "Available options are:\n"
	@echo "    help - show this message"
	@echo "    html - build the html docs"
	@echo "    html-examples - update examples and build html docs"
	@echo "    doctest - run doctests"
	@echo "    test - execute unit tests"
	@echo "    clean - clean all files from docs and pip routines"
	@echo "    test-all - execute testing suite"
	@echo "    pip - publish the package to the PyPi index"
	@echo "    example-plot-dos - update example of plotting DOS"
	@echo "    example-identify-wannier-centres - update example of identifying wannier centres"
	@echo "    example-make-template - update example of making template"
	@echo "    example-plot-tb2j - update example of plotting tb2j"
	@echo "    example-extract-tb2j - update example of extracting tb2j"
	@echo "    example-plot-dos-gallery - update example of plotting DOS gallery"
	@echo "    examples - update all examples"
	@echo "    bravais-pictures - update pictures of bravais lattices"
	@echo "    check-scripts - check consistency of argument names in scripts"
	@echo "    prepare-release - prepare the package for release"
	@echo

example-plot-dos:
	-@rm -r docs/examples/rad-plot-dos/style-examples/*
	@rad-plot-dos.py -ip docs/examples/rad-plot-dos/collinear-spin-polarized/ -ew -7 -2 --custom "Ni (d)" "I (p)" -op docs/examples/rad-plot-dos/style-examples
	@rad-plot-dos.py -ip docs/examples/rad-plot-dos/collinear-spin-polarized/ -ew -7 -2 --custom "Ni (d)" "I (p)" -op docs/examples/rad-plot-dos/style-examples -n
	@rad-plot-dos.py -ip docs/examples/rad-plot-dos/collinear-spin-polarized/ -ew -7 -2 --custom "Ni (d)" "I (p)" -op docs/examples/rad-plot-dos/style-examples -bt
	@rad-plot-dos.py -ip docs/examples/rad-plot-dos/collinear-spin-polarized/ -ew -7 -2 --custom "Ni (d)" "I (p)" -op docs/examples/rad-plot-dos/style-examples -n -bt
	@rad-plot-dos.py -ip docs/examples/rad-plot-dos/collinear-spin-polarized/ -ew -7 -2 --custom "Ni (d)" "I (p)" -op docs/examples/rad-plot-dos/style-examples -r
	@rad-plot-dos.py -ip docs/examples/rad-plot-dos/collinear-spin-polarized/ -ew -7 -2 --custom "Ni (d)" "I (p)" -op docs/examples/rad-plot-dos/style-examples -r -n
	@rad-plot-dos.py -ip docs/examples/rad-plot-dos/collinear-spin-polarized/ -ew -7 -2 --custom "Ni (d)" "I (p)" -op docs/examples/rad-plot-dos/style-examples -r -bt
	@rad-plot-dos.py -ip docs/examples/rad-plot-dos/collinear-spin-polarized/ -ew -7 -2 --custom "Ni (d)" "I (p)" -op docs/examples/rad-plot-dos/style-examples -r -n -bt

example-identify-wannier-centres:
	@rad-identify-wannier-centres.py docs/examples/rad-identify-wannier-centres/example_centres.xyz --span 0.11 --output-name example_centres.xyz_bigger_span

example-make-template:
	@rad-make-template.py -on docs/examples/rad-make-template/template_demo.txt
	@rad-make-template.py -if docs/examples/rad-make-template/exchange.out -on docs/examples/rad-make-template/full_template.txt
	@rad-make-template.py -if docs/examples/rad-make-template/exchange.out -on docs/examples/rad-make-template/filtered_template.txt -maxd 8

example-plot-tb2j:
	@rad-plot-tb2j.py -if docs/examples/rad-plot-tb2j/exchange.out -op docs/examples/rad-plot-tb2j/
	@rad-plot-tb2j.py -if docs/examples/rad-plot-tb2j/exchange.out -op docs/examples/rad-plot-tb2j/  -on exchange_filtered -wtp iso -maxd 5 -dc -sa 1.2 -sd 1.2 -t "First neighbour exchange"
	@rad-plot-tb2j.py -if docs/examples/rad-plot-tb2j/exchange.out -op docs/examples/rad-plot-tb2j/  -on exchange_template -wtp iso -tf docs/examples/rad-plot-tb2j/template.txt -dc -sa 1.2 -sd 1.2 -t "First neighbour exchange"
	@rad-plot-tb2j.py -if docs/examples/rad-plot-tb2j/exchange.out -op docs/examples/rad-plot-tb2j/  -on exchange_forced_symmetry -tf docs/examples/rad-plot-tb2j/template.txt -fs -dc -sa 1.2 -sd 1.2 -t "Forced symmetry exchange"
	@rad-plot-tb2j.py -if docs/examples/rad-plot-tb2j/exchange.out -op docs/examples/rad-plot-tb2j/  -on exchange_R -wtp iso -R 1 0 0 1 1 0 0 1 0 -1 0 0 -1 -1 0 0 -1 0 -dc -sa 1.2 -sd 1.2 -t "First neighbour exchange"
	
example-extract-tb2j:
	@rad-extract-tb2j.py -if docs/examples/rad-extract-tb2j/exchange.out -tf docs/examples/rad-extract-tb2j/template.txt -on docs/examples/rad-extract-tb2j/summary_forced_symmetry.txt -all -fs
	@rad-extract-tb2j.py -if docs/examples/rad-extract-tb2j/exchange.out -tf docs/examples/rad-extract-tb2j/template.txt -on docs/examples/rad-extract-tb2j/summary.txt -all

example-plot-dos-gallery:
	@rad-plot-dos.py -ip docs/examples/rad-plot-dos/collinear/ -ew -6.5 6.5 -ef -1.7806 -op docs/examples/rad-plot-dos/collinear/ 
	@rad-plot-dos.py -ip docs/examples/rad-plot-dos/collinear/ -ew -6.5 6.5 -ef -1.7806 -op docs/examples/rad-plot-dos/collinear/ -r -n
	@rad-plot-dos.py -ip docs/examples/rad-plot-dos/collinear-spin-polarized/ -ew -6.5 6.5 -ef -1.7810 -op docs/examples/rad-plot-dos/collinear-spin-polarized/ 
	@rad-plot-dos.py -ip docs/examples/rad-plot-dos/collinear-spin-polarized/ -ew -6.5 6.5 -ef -1.7810 -op docs/examples/rad-plot-dos/collinear-spin-polarized/ -r -n
	@rad-plot-dos.py -ip docs/examples/rad-plot-dos/noncollinear-nonso/ -ew -6.5 6.5 -ef -1.7810 -op docs/examples/rad-plot-dos/noncollinear-nonso/ 
	@rad-plot-dos.py -ip docs/examples/rad-plot-dos/noncollinear-nonso/ -ew -6.5 6.5 -ef -1.7810 -op docs/examples/rad-plot-dos/noncollinear-nonso/ -r -n
	@rad-plot-dos.py -ip docs/examples/rad-plot-dos/noncollinear-so/ -ew -6.5 6.5 -ef -1.6372 -op docs/examples/rad-plot-dos/noncollinear-so/ 
	@rad-plot-dos.py -ip docs/examples/rad-plot-dos/noncollinear-so/ -ew -6.5 6.5 -ef -1.6372 -op docs/examples/rad-plot-dos/noncollinear-so/ -r -n

examples: example-plot-dos example-identify-wannier-centres example-make-template example-plot-tb2j example-extract-tb2j example-plot-dos-gallery
	@echo "Done"

bravais-pictures:
	@python3 docs/source/user-guide/module/crystal/bravais-lattices/plot_all.py -op docs/source/user-guide/module/crystal/bravais-lattices/

check-scripts:
	@python3 tools/check-scripts.py

VERSION="undefined"
prepare-release:
	@python3 -u tools/prepare-release.py -v $(VERSION)

docs-pictures:
	@python3 tools/plot-data-structure.py