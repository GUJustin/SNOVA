# Reproducibility artifact for the SNOVA symmetry-quotient attack

This artifact accompanies **Symmetry-Quotient Forgery Attacks on Odd-Characteristic SNOVA**. It is deliberately scoped to the claims in the public paper: the exact symmetry quotient, the Version 2.3 rank and affine-elimination certificates, the Level-I two-column 50-in-52 certificate, and the candidate-relative cost tables. Exploratory negative results and cross-scheme research are excluded.

## Requirements

- Python 3.10 or later
- NumPy 1.24 or later

Install the sole Python dependency with:

```bash
python3 -m pip install -r requirements.txt
```

## Checks that do not require a KAT file

```bash
./run_without_kat.sh
```

This reconstructs the fixed public ABQ constants, checks the symmetry-quotient ranks for all nine Version 2.3 q=19 shapes, evaluates the small l=2 format-rejection tails, and regenerates both the Version 2.3 and conditional Version 2.4-preview cost tables.

Expected headline Version 2.3 values under the stricter `p1=1` convention are:

```text
138.94, 130.30, 138.94,
184.29, 184.29, 184.29,
227.95, 238.77, 227.95 bits.
```

## Official-KAT checks

Generate `PQCsignKAT_SNOVA_28_5_19_4.txt` from the official SNOVA repository at the audited commit:

```bash
git clone https://github.com/PQCLAB-SNOVA/SNOVA.git
cd SNOVA
git checkout 9da14981336ede257c41ef53cc069989051e8181
cd dist
make
cd ref
make kat
```

Copy the resulting KAT response file into `source_snapshots/`. Then run:

```bash
./run_with_kat.sh
```

The first script reconstructs the Level-I public key byte-for-byte, verifies the 80-by-50 quotient rank, performs the 30 affine eliminations, checks exact equivalence with the original 80 verifier coordinates, and verifies the deterministic skew-offset format bypass. The second script verifies the four self-column ranks, six invertible cross-column maps, full 80-by-680 feature-map rank, explicit rank-80 `x0` certificate, 52-dimensional affine kernel, 50-dimensional restricted self-quadratic span, and equality with the direct verifier on random two-column assignments.

## Recorded outputs

Files ending in `_recorded.json` are outputs retained from the source-level audit that produced the manuscript. The KAT-independent checks have also been rerun in this release environment. KAT-dependent scripts are included and pinned to the expected public input, but the public KAT response file itself is not redistributed here.

## File map

- `repro/symmetry_attack_validation.py`: quotient, affine elimination, format bypass, and all-parameter ranks.
- `repro/official_estimator.py`: both semi-regular conventions and all cost tables.
- `repro/cross_column_certificate.py`: official Level-I two-column certificate.
- `results/validation_results_recorded.json`: recorded end-to-end KAT audit.
- `results/cross_column_certificate_recorded.json`: positive certificate data used in Section 6.
- `results/full_feature_subset_ranks_recorded.json`: complete official feature-subset rank table.
- `results/affine_rank_survey_recorded.json`: affine-rank survey.
- `results/snova24_symmetry_profiles_recorded.json`: Version 2.4-preview calculations.
- `source_snapshots/SNOVA_2.3.pdf`: audited specification snapshot.
