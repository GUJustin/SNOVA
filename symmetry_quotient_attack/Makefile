LATEXMK=latexmk
PYTHON=python3

.PHONY: all verify regenerate clean

all: paper.pdf

paper.pdf: paper.tex references.bib
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error paper.tex

verify:
	$(PYTHON) artifact/verify_referee_revision.py
	$(PYTHON) artifact/reduced_attack/verify_transcript.py

regenerate:
	$(PYTHON) artifact/generate_ledger.py
	$(PYTHON) artifact/field_circuits.py
	$(PYTHON) artifact/research/verify_l2_channel_conditional_atoms.py
	$(PYTHON) artifact/research/jg_revised_ledger.py
	$(PYTHON) artifact/research/jg_capped_tree_ledger.py
	$(PYTHON) artifact/reduced_attack/attack.py --out artifact/reduced_attack/transcript.json
	$(PYTHON) artifact/reduced_attack/run_batch.py

clean:
	$(LATEXMK) -C paper.tex
