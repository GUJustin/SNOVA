LATEXMK=latexmk

.PHONY: all clean

all: paper.pdf

paper.pdf: paper.tex references.bib
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error paper.tex

clean:
	$(LATEXMK) -C paper.tex
