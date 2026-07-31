LATEXMK=latexmk
PYTHON=python3

.PHONY: all verify verify-official-inputs verify-ledgers verify-correspondence verify-reduced-e2e verify-reduced-stress regenerate clean

all: paper.pdf

paper.pdf: paper.tex references.bib
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error paper.tex

verify: verify-official-inputs verify-ledgers verify-correspondence verify-reduced-e2e verify-reduced-stress

verify-official-inputs:
	cd artifact/official && sha256sum -c SHA256SUMS

verify-ledgers:
	$(PYTHON) artifact/verify_artifact.py

verify-correspondence:
	$(PYTHON) artifact/correspondence/official_key_reduction_harness.py --check

verify-reduced-e2e:
	$(PYTHON) artifact/correspondence/reduced_parameter_end_to_end_forgery.py --check

verify-reduced-stress:
	$(PYTHON) artifact/correspondence/verify_reduced_parameter_stress.py

regenerate:
	$(PYTHON) artifact/generate_ledger.py
	$(PYTHON) artifact/field_circuits.py
	$(PYTHON) artifact/research/verify_l2_channel_conditional_atoms.py
	$(PYTHON) artifact/research/jg_expected_ledger.py
	$(PYTHON) artifact/research/jg_capped_tree_ledger.py
	$(PYTHON) artifact/correspondence/official_key_reduction_harness.py --write
	$(PYTHON) artifact/correspondence/reduced_parameter_end_to_end_forgery.py --write

clean:
	$(LATEXMK) -C paper.tex
