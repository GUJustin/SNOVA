# Official input snapshots

`SNOVA_2.3.pdf` is the Version 2.3 draft used by the audit. Its SHA-256 digest is

```
963fa43d7e6ed8ecf1e77556e7cf47f62d4c34177efece2cd9a9517d0fbaf783
```

The Level-I KAT response file is not redistributed in this package. Obtain or generate the public file

```
PQCsignKAT_SNOVA_28_5_19_4.txt
```

from the official SNOVA repository at audited commit

```
9da14981336ede257c41ef53cc069989051e8181
```

Generate the recommended instances and KAT files using the repository's documented commands:

```bash
git clone https://github.com/PQCLAB-SNOVA/SNOVA.git
cd SNOVA
git checkout 9da14981336ede257c41ef53cc069989051e8181
cd dist
make
cd ref
make kat
```

Copy the resulting `PQCsignKAT_SNOVA_28_5_19_4.txt` into this directory before running the KAT-dependent commands in the top-level README. This keeps the artifact small while pinning the exact official input expected by the scripts.
