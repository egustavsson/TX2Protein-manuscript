## Protein Group Filtering Scripts

These scripts allow a user to extract relevant entries from a set of MetaMorpheus quantified protein group results files.

### Input Requirements

- `transcrpts.tsv`: A single-column TSV file listing novel transcripts identified prior to running MetaMorpheus.

### Step 1: Filter for Relevant Protein Groups

```bash
python3 filter_proteins.py \
  --id_file ./transcrpts.tsv \
  --output_file output.tsv \
  --tsv_files \
    /screen/MSV000086944-GM12878/Task3SearchTask/AllQuantifiedProteinGroups.tsv \
    /screen/MSV000087506/Task3SearchTask/AllQuantifiedProteinGroups.tsv \
    /screen/MSV000088006/Task3SearchTask/AllQuantifiedProteinGroups.tsv \
    /screen/MSV000091887/Task3SearchTask/AllQuantifiedProteinGroups.tsv \
    /screen/MSV000086944-20151201/Task3SearchTask/AllQuantifiedProteinGroups.tsv
```

### Step 2 (Optional): Filter Non-Unique Peptides

MetaMorpheus may mark a peptide as unique even when its amino acid sequence appears in multiple proteins. This behavior is technically valid but complex. These scripts simplify this by filtering only for strictly unique sequences.

- Input: Reference proteome FASTA file.
- Effect: Removes peptides that match exactly with any known protein sequence.

```bash
python3 filter_known_peptides.py \
  ../lr_pipeline/uniprotkb_proteome_UP000005640_2023_09_12.fasta \
  ./output.tsv \
  ./filtered_known_peptides.tsv
```

### Step 3: Summarize Novel Transcripts

- Input: GTF file containing transcript definitions.

```bash
python3 summarize_transcripts.py \
  ./novel.gtf \
  ./filtered_known_peptides.tsv \
  ./summary.tsv
```
