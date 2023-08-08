
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TX2Protein - Transcript to protein

## Background

This repository contains the code used for the analysis in the manuscript *Protein isoform detection from long-read RNA-seq data*. This manuscript describes our software \[NAME OF SOFTWARE\], a simplified method for intergating protein coding potential of full-length transcript structures in mass spectrometry data. \[NAME OF SOFTWARE\] only requires a `GTF/GFF` file and the raw proteomic dataset, typically from The ProteomeXchange Consortium, and uses [`MetaMorpheus`](https://github.com/smith-chem-wisc/MetaMorpheus) for protein identifications.

## Data

The data used in the

### 

These analyses were performed using publicly available long-read RNA seq data, consisting of PacBio IsoSeq from 9 frontal cortex samples.

### ENCODE long-read data

This consists of publicly available frontal cortex samples from [ENCODE](https://www.encodeproject.org/rna-seq/long-read-rna-seq/). All samples were sequenced on the PacBio Sequel II platform and processed with the ENCODE DCC deployment of the [TALON pipeline](https://github.com/ENCODE-DCC/long-read-rna-pipeline). We downloaded the TSV and GTF files from the following samples:

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

This data can be downloaded using this [script](./scripts/00a_download_ENCODE_files.sh). It will calculate the total size of the download and allow a Y/N choice prior to downloading. Usage:

    # first give the execute permission to the script:
    chmod +x 00a_download_ENCODE_files.sh

    # then run it:
    ./00a_download_ENCODE_files.sh /path/to/directory/

This data can be filtered using this [script](./scripts/00b_filter_ENCODE_long_read.sh). Usage:

    $ bash 00b_filter_ENCODE_long_read.sh --directory /path/to/directory/ --output /path/to/output/ --gene ENSG00000136156.13

## Code contents

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
<td><a href="docs" class="uri">docs</a></td>
<td>Contains all <code>.Rmd</code>s and their corresponding <code>.html</code>s describing analyses performed for this project.</td>
</tr>
<tr class="even">
<td><a href="logs" class="uri">logs</a></td>
<td>For any scripts that were run outside of an <code>.Rmd</code> (e.g. scripts from the <a href="scripts" class="uri">scripts</a> directory), a log file was recorded and can be accessed here.</td>
</tr>
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
