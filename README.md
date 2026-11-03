# durg_domain_interaction

Exploratory experiments on inferring drug-domain interactions from drug-protein interaction data.
The question: if a drug binds a single-domain protein with domain X, does it bind other proteins that
contain X? The notebook looks for counter-examples ("problem 1") using negative interaction data, and
sketches the harder multi-domain case ("problem 2").

Everything lives in one notebook, `initial_experiments.ipynb`, which walks through:

1. **First attempt (Coelho 2016 dataset).** Reads positive/negative drug-target pairs from
   `DTIPred/` (`drugbank_DTIs_REAL_NEGS.txt`, `yamanishi_DTIs_REAL_NEGS.txt`, `test_data_sc_and_bc.txt`)
   and Pfam annotations exported from UniProt (`DTIPred/uniprotids_annnots.tab`). For each
   single-domain protein P and drug D it bins the other proteins sharing P's domain into pos/neg/unk and
   writes `result_drug_level.csv`. No counter-examples were found, attributed to too few negatives.
2. **Second attempt (BindingDB).** Reads `data/BindingDB_All_2021m0.tsv/BindingDB_All.tsv`, keeps
   single-chain rows with a UniProt ID and PubChem CID, and binarizes Ki / IC50 (nM): below
   `interaction_threshold = 1000` is an interaction, above `noninteraction_threshold = 30000` is a
   non-interaction. Domain layouts are fetched per protein from the Pfam XML API and cached in
   `temp/domain_data/<uniprot>.xml`; a protein counts as single-domain if it has one Pfam match and no
   domain-free stretch longer than `blank_length_threshold = 30` residues. The parsed interaction
   tables are pickled to `temp/pickled_data/data.pickle` so the slow BindingDB pass runs once
   (`reload_dataset = True` forces a rerun).
3. **Conflict removal.** `rc1` drops drug-domain pairs whose negative protein is itself single-domain;
   `rc2` is stricter and drops pairs where a positive and a negative protein have the same
   domain/non-domain segment layout.
4. **Report.** For the surviving pairs it writes `outputs/problem_1_examples_rc2/doc.md` with UniProt
   tables of positive and negative proteins and an SVG domain-layout figure per pair.
5. **Dataset 2.** Extends dataset 1 with all negatives (and positives) that share a drug or protein
   with a dataset-1 positive pair.

## Running

```
pip install -r requirements.txt
jupyter notebook initial_experiments.ipynb
```

Run the cells top to bottom. Note that the two attempts define the same names (`interacts`/`interact`,
`drugsof`, `domainsof`, ...) from different sources; the second attempt supersedes the first.

## Data not included

The repository contains only the notebook. You need to provide:

- `DTIPred/` with the three interaction files above and `uniprotids_annnots.tab` (the notebook writes
  `DTIPred/uniprotids.txt`, the ID list to upload to UniProt's ID mapping to produce that table);
- `data/BindingDB_All_2021m0.tsv/BindingDB_All.tsv`, the full BindingDB TSV export;
- the directories `temp/domain_data/`, `temp/pickled_data/` and `outputs/` must exist before the cells
  that write into them.

The Pfam records are downloaded from `http://pfam.xfam.org/protein/<id>?output=xml`. That legacy Pfam
site has since been retired in favour of InterPro, so the download step will need to be pointed at a
current source if `temp/domain_data/` is not already populated.
