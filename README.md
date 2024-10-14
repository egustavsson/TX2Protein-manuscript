
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TX2Protein - Transcript to protein

## Background

This repository contains the code used for the analysis in the manuscript *Protein isoform detection from long-read RNA-seq data*. This manuscript describes our software [`TX2P`](https://github.com/MurphyDavid/TX2P.git), a simplified method for intergating protein coding potential of full-length transcript structures in mass spectrometry data. `TX2P` only requires a `GTF/GFF` file and the raw proteomic dataset, typically from The ProteomeXchange Consortium, and uses [`MetaMorpheus`](https://github.com/smith-chem-wisc/MetaMorpheus) for protein identifications.

## Data

### Gene set

We used the the Early onset or syndromic epilepsy (Version 4.77) from Genomics England PanelApp. We only included Green genes, these are the genes with the highest vidence. This panel is used as a virtual panel to analyse genome or exome data in the NHS Genomic Medicine Service; the panel will routinely be applied for clinical indication 'R59 Early onset or syndromic epilepsy' but can also be used as part of the analysis for a broader clinical presentation, where relevant.

This panel is available through <https://panelapp.genomicsengland.co.uk/panels/402/>

### Long-read RNA-seq

This consists of 9 publicly available frontal cortex samples from [ENCODE](https://www.encodeproject.org/rna-seq/long-read-rna-seq/). All samples were sequenced on the PacBio Sequel II platform and processed with the ENCODE DCC deployment of the [TALON pipeline](https://github.com/ENCODE-DCC/long-read-rna-pipeline). We downloaded the TSV and GTF files from the following samples:

<table style="width:100%;">
<colgroup>
<col width="20%" />
<col width="23%" />
<col width="26%" />
<col width="9%" />
<col width="21%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Sample</th>
<th align="left">Disease status</th>
<th align="left">Tissue</th>
<th align="left">Sex</th>
<th align="left">Age</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><a href="https://www.encodeproject.org/experiments/ENCSR462COR/">ENCSR462COR</a></td>
<td align="left">Alzheimer's disease</td>
<td align="left">middle frontal area 46</td>
<td align="left">female</td>
<td align="left">90 or above years</td>
</tr>
<tr class="even">
<td align="left"><a href="https://www.encodeproject.org/experiments/ENCSR169YNI/">ENCSR169YNI</a></td>
<td align="left">Control</td>
<td align="left">middle frontal area 46</td>
<td align="left">female</td>
<td align="left">90 or above years</td>
</tr>
<tr class="odd">
<td align="left"><a href="https://www.encodeproject.org/experiments/ENCSR257YUB/">ENCSR257YUB</a></td>
<td align="left">Alzheimer's disease</td>
<td align="left">middle frontal area 46</td>
<td align="left">male</td>
<td align="left">90 or above years</td>
</tr>
<tr class="even">
<td align="left"><a href="https://www.encodeproject.org/experiments/ENCSR690QHM/">ENCSR690QHM</a></td>
<td align="left">Alzheimer's disease</td>
<td align="left">middle frontal area 46</td>
<td align="left">female</td>
<td align="left">90 or above years</td>
</tr>
<tr class="odd">
<td align="left"><a href="https://www.encodeproject.org/experiments/ENCSR316ZTD/">ENCSR316ZTD</a></td>
<td align="left">Alzheimer's disease</td>
<td align="left">middle frontal area 46</td>
<td align="left">female</td>
<td align="left">90 or above years</td>
</tr>
<tr class="even">
<td align="left"><a href="https://www.encodeproject.org/experiments/ENCSR697ASE/">ENCSR697ASE</a></td>
<td align="left">Alzheimer's disease</td>
<td align="left">middle frontal area 46</td>
<td align="left">female</td>
<td align="left">90 or above years</td>
</tr>
<tr class="odd">
<td align="left"><a href="https://www.encodeproject.org/experiments/ENCSR094NFM/">ENCSR094NFM</a></td>
<td align="left">Control</td>
<td align="left">middle frontal area 46</td>
<td align="left">female</td>
<td align="left">90 or above years</td>
</tr>
<tr class="even">
<td align="left"><a href="https://www.encodeproject.org/experiments/ENCSR463IDK/">ENCSR463IDK</a></td>
<td align="left">Control</td>
<td align="left">middle frontal area 46</td>
<td align="left">female</td>
<td align="left">79</td>
</tr>
<tr class="odd">
<td align="left"><a href="https://www.encodeproject.org/experiments/ENCSR205QMF/">ENCSR205QMF</a></td>
<td align="left">Control</td>
<td align="left">middle frontal area 46</td>
<td align="left">female</td>
<td align="left">88</td>
</tr>
</tbody>
</table>

## Analysis steps

### 1. Get gene set

Downloading the desired panel from PanelApp, and generation of the `genes_set.tsv` used as input later, can be done by running [`00a_get_gene_set.sh`](./scripts/00a_get_gene_set.sh). Usage:

    bash 00a_get_gene_set.sh /path/to/output/directory "https://panelapp.genomicsengland.co.uk/panels/402/download/34/"

This will store the gene panel and a `genes_set.tsv` file which looks like this:

``` bash
ENSG00000078369
ENSG00000187730
ENSG00000198793
```

**IMPORTANT!** This script filters on gene names and assumes this information is contained in column 3.

### 2. Get and process long-read RNA-seq data

This data can be downloaded using this [`00b_download_ENCODE_files.sh`](./scripts/00b_download_ENCODE_files.sh). It will calculate the total size of the download and allow a Y/N choice prior to downloading. Usage:

    # first give the execute permission to the script:
    chmod +x 00b_download_ENCODE_files.sh

    # then run it:
    00b_download_ENCODE_files.sh /path/to/directory/

To filter the data to only include the genes of interest (output from `genes_set.tsv` generated by `00a_get_gene_set.sh`) run [`00c_filter_ENCODE_long_read.sh`](./scripts/00c_filter_ENCODE_long_read.sh). Usage:

    bash 00c_filter_ENCODE_long_read.sh --directory /path/to/input/directory/ --output /path/to/output/directory/ --gene-file genes_set.tsv

### 3. Get unique transcripts across samples

When retrieving novel transcripts through long-read RNA-seq, their IDs are generated in the order they are called. Consequently, IDs for the same transcript can vary among different samples. To ensure a consistent unique ID for all transcripts across samples (in `GTF` files), we execute run [`01a_process_ENCODE_transcript_structures.R`](./scripts/01a_process_ENCODE_transcript_structures.R) to get a set of unique transcripts.

## Contents

Within this repository you will find:

<table>
<colgroup>
<col width="11%" />
<col width="88%" />
</colgroup>
<thead>
<tr class="header">
<th>Directory</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><a href="raw_data" class="uri">raw_data</a></td>
<td>Data used for the analysis. Most will not be available due to size.</td>
</tr>
<tr class="even">
<td><a href="results" class="uri">results</a></td>
<td>Results from all analyses.</td>
</tr>
<tr class="odd">
<td><a href="scripts" class="uri">scripts</a></td>
<td>Contains analysis scripts. Each script contains a one-line description and is also referenced in its corresponding <code>.Rmd</code>.</td>
</tr>
</tbody>
</table>
