LATEXMK=latexmk
PYTHON=python3

.PHONY: all verify regenerate clean

all: paper.pdf

paper.pdf: paper.tex references.bib
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error paper.tex

verify:
	$(PYTHON) artifact/verify_final_results.py

regenerate:
	$(PYTHON) artifact/compute_final_results.py

clean:
	$(LATEXMK) -C paper.tex
