# The Symmetry-Quotient Attack on Odd-Characteristic SNOVA

This directory contains the complete standalone concrete cryptanalysis paper.

Build with:

```sh
make
```

Run the independent deterministic checks with:

```sh
make verify
```

Regenerate the machine-readable numerical ledgers with:

```sh
make regenerate
```

The manuscript has exactly one LaTeX source file, `paper.tex`. The `artifact/`
directory contains the numerical and circuit certificates referenced by the
paper.  The checker uses only Python's standard library and requires no network
access.
