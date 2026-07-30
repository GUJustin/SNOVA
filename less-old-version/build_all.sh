#!/bin/sh
set -eu
latexmk -pdf -interaction=nonstopmode -halt-on-error snova-symmetry-quotient.tex
cp snova-symmetry-quotient.pdf snova-symmetry-quotient-public.pdf
printf '%s\n' '\def\ANONYMOUS{1}' '\input{snova-symmetry-quotient.tex}' > build-anonymous.tex
pdflatex -interaction=nonstopmode -halt-on-error build-anonymous.tex >/dev/null
biber build-anonymous >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error build-anonymous.tex >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error build-anonymous.tex >/dev/null
cp build-anonymous.pdf snova-symmetry-quotient-anonymous.pdf
printf '%s\n' '\def\DISCLOSURE{1}' '\input{snova-symmetry-quotient.tex}' > build-disclosure.tex
pdflatex -interaction=nonstopmode -halt-on-error build-disclosure.tex >/dev/null
biber build-disclosure >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error build-disclosure.tex >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error build-disclosure.tex >/dev/null
cp build-disclosure.pdf snova-symmetry-quotient-disclosure.pdf
