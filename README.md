# Nanopore Variant Analysis Pipeline

This pipeline is designed for analyzing amplicon-based nanopore sequencing for CRISPR mediated mutations. It supports read alignment (using minimap2), variant calling (via bcftools), and basic variant summarization, providing insights into reference matches, insertions, deletions, and substitutions. Note, for the summary statistics output file, the reference percentages assume a single variant per read in the target site, and thus if multiple variants are present in a single read that impact the target site, the percent of reference bases will be underestimated.

## Installation Guide

### System Requirements

This pipeline has been tested on Ubuntu Linux with the following external tool dependencies:

- Minimap2 (2.26-r1175 tested)
- BCFtools (1.10.2 tested)
- Samtools (1.10 tested)
- Biopython
- pysam

Ensure Minimap2, BCFtools, and Samtools are installed and accessible in your system's PATH. This pipeline has been tested on Ubuntu Linux.

### Setting Up Python Environment

Create a virtual environment and install Python dependencies:

```sh
python3 -m venv nanopore_env
source nanopore_env/bin/activate  # On Windows use `nanopore_env\Scripts\activate`
pip install biopython pysam
```
### Installing Command-line Tools

- **Minimap2**, **BCFtools**, and **Samtools**: These can be installed via package managers on most Linux distributions. For Ubuntu, use:

```bash
sudo apt-get update
sudo apt-get install minimap2 bcftools samtools
```
## Usage

To run the pipeline, use the following command:

```bash
python Nanopore_Variant_analysis.py --fastq_path <path_to_fastq> --ref_fasta <path_to_ref_fasta> --output_dir <output_directory> --targets <path_to_targets> --min_read_num <min_read> --min_read_percentage <min_percentage> --min_read_quality <min_quality> --min_alignment_quality <min_align_quality> --max_depth <max_depth>
```
### Parameters

- `--fastq_path <path>`: Specify the path to the FASTQ file or directory containing FASTQ files.
- `--ref_fasta <path>`: Path to the reference genome FASTA file.
- `--output_dir <path>`: Directory where output files will be saved.
- `--targets <path>`: Path to a TSV file containing target sites for variant analysis.
- `--min_read_num <number>`: Minimum number of reads supporting a variant (default: 3).
- `--min_read_percentage <percentage>`: Minimum percentage of reads that support a variant, expressed as a percentage (e.g., for 10%, use 0.1 as the value; default: 0.1).
- `--min_read_quality <quality score>`: Minimum quality score for reads to be considered in the analysis (default: 20).
- `--min_alignment_quality <quality score>`: Minimum alignment quality score for alignments to be considered (default: 10).
- `--max_depth <depth>`: Maximum read depth to consider in variant calling (default: 10000).

### Example

```sh

python Nanopore_Variant_analysis.py --fastq_path barcode01 --ref_fasta SRSF2_offtargets.fasta --output_dir barcode01_out --targets Targets.txt --min_read_num 3 --min_read_percentage 0.1 --min_read_quality 16 --min_alignment_quality 16 --max_depth 10000
```

## Contributing

Contributions to improve the pipeline or extend its functionality are welcome. Please submit pull requests or report issues via GitHub.
