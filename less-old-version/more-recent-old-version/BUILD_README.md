# Build instructions

The paper is self-contained in:

- `snova-paper-final-mq-accessible.tex`
- `references.bib`

Build with a current TeX Live installation containing `biblatex` and Biber:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  snova-paper-final-mq-accessible.tex
```

Optional source branches:

```tex
\def\ANONYMOUS{}
\input{snova-paper-final-mq-accessible.tex}
```

or

```tex
\def\DISCLOSURE{}
\input{snova-paper-final-mq-accessible.tex}
```

The delivered PDF was built with pdfTeX 1.40.26 and Biber 2.20.
