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
	-@rm -r docs/source/api/_autosummary
	-@rm -r rad_tools.egg-info
	-@rm -r build
	-@rm -r dist

test:
	@pip3 install . --upgrade
	@pip3 install pytest
	@pytest -s

VERSION := $(shell grep __version__ rad_tools/__init__.py | tr -d 'a-zA-Z =_":')
.ONESHELL:
pip:
#	@echo "\x1b[33m"
#	@echo "pip is disabled for your own safety"
#	@echo "\x1b[0m"
	@echo "\x1b[33m"
	@echo "Have you done this?"
	@echo "  * Change version in __init__.py?"
	@echo "\x1b[0m"
	@grep "__version__" rad_tools/__init__.py
	@echo "\x1b[33m"
	@echo "  * Release note?"
	@echo "  * Commit all the changes?"
	@echo "\x1b[0m"
	@git status
	@echo "\x1b[33m"
	@echo "  * Merge to the stable?"
	@echo "  * Push to GitHub?"
	@echo "  * Create a tag?"
	@echo "\x1b[0m"
	@git log --oneline --all --graph --decorate -10
	@echo "\x1b[33m"
	@echo "\x1b[31m"
	@echo "\x1b[33m"
	@echo Are all script arguments ok?
	@echo "\x1b[0m"
	@make check-script-names
	@read -p "Press Enter if yes:"
	@echo "\x1b[0m"
	-@rm -r dist
	-@rm -r build
	-@rm -r rad_tools.egg-info
	@python3 setup.py sdist bdist_wheel 
	@python3 -m twine upload --repository pypi dist/* --verbose
	@git tag -a "$(VERSION)" -m "Version $(VERSION)"


help:
	@echo "\x1b[31m"
	@echo "Please specify what do you want to do!"
	@echo "\x1b[32m"
	@echo "Available options are:"
	@echo "    html - build the html docs"
	@echo "    html-examples - update examples and build html docs"
	@echo "    doctest - run doctests"
	@echo "    clean - clean all files from docs and pip routines"
	@echo "    test - execute testing suite"
	@echo "    pip - publish the package to the PyPi index"
	@echo "    examples - update all examples"
	@echo "    push - update examples and git push"
	@echo "    check-script-names - check consistency of argument names in scripts"
	@echo "\x1b[0m"

examples:
	@pip3 install . --upgrade
	@identify-wannier-centres.py docs/examples/identify-wannier-centres/example_centres.xyz -nc > docs/examples/identify-wannier-centres/console_output.txt
	@identify-wannier-centres.py docs/examples/identify-wannier-centres/example_centres.xyz --span 0.11 --output-name example_centres.xyz_bigger_span
	@rad-make-template.py -on docs/examples/rad-make-template/template_demo
	@rad-make-template.py -if docs/examples/rad-make-template/exchange.out -on docs/examples/rad-make-template/full_template
	@rad-make-template.py -if docs/examples/rad-make-template/exchange.out -on docs/examples/rad-make-template/filtered_template -maxd 8
	@tb2j-plotter.py -if docs/examples/tb2j-plotter/exchange.out -op docs/examples/tb2j-plotter/
	@tb2j-plotter.py -if docs/examples/tb2j-plotter/exchange.out -op docs/examples/tb2j-plotter/  -on exchange_filtered -wtp iso -maxd 5 -dc -sa 1.2 -sd 1.2 -t "First neighbour exchange"
	@tb2j-plotter.py -if docs/examples/tb2j-plotter/exchange.out -op docs/examples/tb2j-plotter/  -on exchange_template -wtp iso -tf docs/examples/tb2j-plotter/template.txt -dc -sa 1.2 -sd 1.2 -t "First neighbour exchange"
	@tb2j-plotter.py -if docs/examples/tb2j-plotter/exchange.out -op docs/examples/tb2j-plotter/  -on exchange_forced_symmetry -tf docs/examples/tb2j-plotter/template.txt -fs -dc -sa 1.2 -sd 1.2 -t "Forced symmetry exchange"
	@tb2j-plotter.py -if docs/examples/tb2j-plotter/exchange.out -op docs/examples/tb2j-plotter/  -on exchange_R -wtp iso -R 1 0 0 1 1 0 0 1 0 -1 0 0 -1 -1 0 0 -1 0 -dc -sa 1.2 -sd 1.2 -t "First neighbour exchange"
	@tb2j-extractor.py -if docs/examples/tb2j-extractor/exchange.out -tf docs/examples/tb2j-extractor/template.txt -op docs/examples/tb2j-extractor/ -on summary_forced_symmetry -all -fs
	@tb2j-extractor.py -if docs/examples/tb2j-extractor/exchange.out -tf docs/examples/tb2j-extractor/template.txt -op docs/examples/tb2j-extractor/ -on summary -all

push: examples
	@git push

check-script-names:
	@python3 check-script-names.py
