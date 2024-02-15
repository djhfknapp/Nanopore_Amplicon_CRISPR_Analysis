# export PATH="$PATH:/mnt/c/minimap2-2.26_x64-linux/" ###Export the path to your minimap install if it isn't already in path
import argparse
import glob
import gzip
import os
import pysam
from Bio import SeqIO
import subprocess
from numpy import median
def pool_fastq_files(fastq_files, output_file):
    """ Pool all sequences from fastq.gz files into a single fastq file """
    with gzip.open(output_file, 'wt') as pooled_file:
        for file in fastq_files:
            with gzip.open(file, 'rt') as f:
                for record in SeqIO.parse(f, "fastq"):
                    SeqIO.write(record, pooled_file, "fastq")

def align_reads(pooled_fastq, fasta_file, sorted_bam_output):
    """ Align reads to reference sequences using Minimap2 and pipe to Samtools for sorting and indexing """
    # Command to run Minimap2 and pipe to Samtools for conversion to BAM, sorting, and indexing
    minimap2_cmd = f"minimap2 -ax map-ont {fasta_file} {pooled_fastq}"
    samtools_view_cmd = "samtools view -b"
    samtools_sort_cmd = f"samtools sort -o {sorted_bam_output}"
    samtools_index_cmd = f"samtools index {sorted_bam_output}"

    # Chain the commands with a pipe
    full_cmd = f"{minimap2_cmd} | {samtools_view_cmd} | {samtools_sort_cmd}; {samtools_index_cmd}"

    # Execute the command
    subprocess.run(full_cmd, shell=True, check=True)

def call_variants(bam_file, ref_fasta, output_vcf, min_base_quality, min_map_quality, max_depth):
    """
    Calls variants using bcftools.
    :param bam_file: Path to the sorted and indexed BAM file.
    :param ref_fasta: Path to the reference fasta file.
    :param output_vcf: Path to the output VCF file.
    """
    cmd = f"bcftools mpileup -Ov -f {ref_fasta} -Q {min_base_quality} -a FORMAT/DP,FORMAT/AD -d {max_depth} --min-MQ {min_map_quality} {bam_file} > {output_vcf} "
    subprocess.run(cmd, shell=True)

def load_targeted_sites(target_sites_file):
    """
    Loads targeted sites from a TSV file.

    :param target_sites_file: Path to the TSV file containing target sites.
    :return: Dictionary with reference name as key and target sequence as value.
    """
    target_sites = {}
    with open(target_sites_file, 'r') as file:
        for line in file:
            ref_name, target_seq = line.strip().split('\t')
            target_sites[ref_name] = target_seq
    return target_sites

def find_target_site_positions(ref_fasta, target_sites):
    """
    Finds the positions of target sites in the reference fasta and adjusts them to 1-based indexing.
    :param ref_fasta: Path to the reference fasta file.
    :param target_sites: Dictionary of target sequences keyed by reference name.
    :return: Dictionary with reference name as key and target positions (start, end) as value.
    """
    fasta = pysam.FastaFile(ref_fasta)
    target_positions = {}
    for ref_name, target_seq in target_sites.items():
        sequence = fasta.fetch(ref_name)
        start_position = sequence.find(target_seq)
        if start_position != -1:
            # Adjusting from 0-based to 1-based indexing
            start_position += 1
            end_position = start_position + len(target_seq) - 1
            target_positions[ref_name] = (start_position, end_position)
    return target_positions

def analyze_variants(vcf_file, ref_fasta, target_positions, min_read_num, min_read_percentage):
    """
    Analyzes variants in a VCF file, focusing on target sites and calculating specific statistics.

    :param vcf_file: Path to the VCF file.
    :param ref_fasta: Path to the reference fasta file.
    :param target_sites_file: Path to the TSV file containing target sites.
    :param min_read_num: Minimum number of reads supporting a variant.
    :param min_read_percentage: Minimum percentage of total reads for a variant.
    :return: Summary statistics and path to the filtered VCF file.
    """

    # Load the VCF file
    vcf = pysam.VariantFile(vcf_file)

    # Prepare to write filtered variants to a new VCF file
    tsv_file = vcf_file.replace('.vcf', '_filtered.tsv')

    # Initialize summary statistics
    summary_stats = {
        ref: {
            'on_target_ref': [], 'on_target_ins': [], 'on_target_del': [], 'on_target_sub': [],
            'off_target_ref': [], 'off_target_ins': [], 'off_target_del': [], 'off_target_sub': [],
            'total_depth_on_target': [], 'total_depth_off_target': []
        } for ref in target_positions
    }
    with open(tsv_file, 'w') as tsv:
        tsv.write("ReferenceGene\tSite\tRefBase\tAltAlleles\tOnTarget\tTotalReadDepth\tAltAlleleDepths\n")
        for record in vcf:
            ref_name = record.chrom
            if ref_name not in target_positions:
                continue

            target_start, target_end = target_positions[ref_name]
            is_in_target = target_start <= record.pos <= target_end
            is_deletion_overlap_target = 'DEL' in record.alts and record.pos < target_start and (record.pos + len(record.ref) - 1) >= target_start

            # Calculate total read depth and allele depths
            total_depth = record.samples[0]['DP']  # Sample-specific DP
            allele_depths = record.samples[0]['AD']
            ref_depth = allele_depths[0]

            if is_in_target or is_deletion_overlap_target:
                summary_stats[ref_name]['total_depth_on_target'].append(total_depth)
            else:
                summary_stats[ref_name]['total_depth_off_target'].append(total_depth)

            # Check if variant passes thresholds and update statistics
            for alt_allele, alt_depth in zip(record.alts, allele_depths[1:]):
                if alt_depth >= min_read_num and (alt_depth / total_depth) * 100 >= min_read_percentage:
                    percentage = (alt_depth / total_depth) * 100
                    tsv.write(f"{ref_name}\t{record.pos}\t{record.ref}\t{alt_allele}\t"
                            f"{'YES' if is_in_target else 'NO'}\t{total_depth}\t{alt_depth}\n")
                    if is_in_target or is_deletion_overlap_target:
                        if len(alt_allele) > len(record.ref):  # Insertion
                            summary_stats[ref_name]['on_target_ins'].append(percentage)
                        elif len(alt_allele) < len(record.ref):  # Deletion
                            summary_stats[ref_name]['on_target_del'].append(percentage)
                        elif alt_allele != '<*>':  # Substitution
                            summary_stats[ref_name]['on_target_sub'].append(percentage)
                    # else:
                    #     if len(alt_allele) > len(record.ref):  # Insertion
                    #         summary_stats[ref_name]['off_target_ins'].append(percentage)
                    #     elif len(alt_allele) < len(record.ref):  # Deletion
                    #         summary_stats[ref_name]['off_target_del'].append(percentage)
                    #     elif alt_allele != '<*>':  # Substitution
                    #         summary_stats[ref_name]['off_target_sub'].append(percentage)
    
    # Adjust summary statistics to percentages
    for stats in summary_stats.values():
        stats['on_target_ref'] = 100 - (sum(stats['on_target_ins']) + sum(stats['on_target_del']) + sum(stats['on_target_sub']))
        # stats['off_target_ref'] = 100 - (sum(stats['off_target_ins']) + sum(stats['off_target_del']) + sum(stats['off_target_sub']))

    # Write summary statistics to TSV file
    summary_tsv_file = vcf_file.replace('.vcf', '_summary.tsv')
    with open(summary_tsv_file, 'w') as tsv:
        tsv.write("Reference\tMedianDepthOnTarget\tMedianDepthOffTarget\tOnTargetRef%\tOnTargetIns%\tOnTargetDel%\tOnTargetSub%\t"
                  "OffTargetRef%\tOffTargetIns%\tOffTargetDel%\tOffTargetSub%\n")
        for ref, stats in summary_stats.items():
            tsv.write(f"{ref}\t{median(stats['total_depth_on_target'])}\t{median(stats['total_depth_off_target'])}\t{stats['on_target_ref']}\t{sum(stats['on_target_ins'])}\t"
                      f"{sum(stats['on_target_del'])}\t{sum(stats['on_target_sub'])}\n")
    
    return summary_stats

def variant_summary(vcf_file):
    """
    Provide a per-site summary of a VCF file.

    :param vcf_file: Path to the VCF file.
    :return: Number of high quality reads supporting reference, insertion, deletion or substitution from the VCF file.
    """

    # Load the VCF file
    vcf = pysam.VariantFile(vcf_file)

    # Prepare to write filtered variants to a new VCF file
    tsv_file = vcf_file.replace('.vcf', '_PerSiteSummary.tsv')

    with open(tsv_file, 'w') as tsv:
        tsv.write("ReferenceGene\tSite\tRefBase\tTotalHighQualityReads\tReferenceReads\tInsertions\tDeletions\tSubstitutions\n")
        for record in vcf:
            ref_name = record.chrom
            # Calculate total read depth and allele depths
            total_depth = record.samples[0]['DP']  # Sample-specific DP
            allele_depths = record.samples[0]['AD']
            ref_depth = allele_depths[0]
            insertion=0
            deletion=0
            substitution=0
            for alt_allele, alt_depth in zip(record.alts, allele_depths[1:]):
                if len(alt_allele) > len(record.ref):  # Insertion
                    insertion+=alt_depth
                elif len(alt_allele) < len(record.ref):  # Deletion
                    deletion+=alt_depth
                elif alt_allele != '<*>':  # Substitution
                    substitution+=alt_depth
            tsv.write(f"{record.chrom}\t{record.pos}\t{record.ref}\t{total_depth}\t{ref_depth}\t{insertion}\t{deletion}\t{substitution}\n")

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Oxford Nanopore Sequencing Data Analysis Pipeline')
    parser.add_argument('--fastq_path', type=str, help='Path to the FASTQ file or directory containing FASTQ files')
    parser.add_argument('--ref_fasta', type=str, help='Path to the reference fasta file')
    parser.add_argument('--output_dir', type=str, help='Directory to store output files')
    parser.add_argument('--targets', type=str, help='Path to tsv containing the target sites for each gene in the reference')
    parser.add_argument('--min_read_num', type=int, default=3, help='Minimum read number to support a variant')
    parser.add_argument('--min_read_percentage', type=float, default=0.1, help='Minimum percent of reads that support a given variant, 0-100 scale')
    parser.add_argument('--min_read_quality', type=float, default=18.0, help='Minimum read quality score')
    parser.add_argument('--min_alignment_quality', type=float, default=10.0, help='Minimum alignment quality score')
    parser.add_argument('--max_depth', type=int, default=10000, help='Maximum read depths used in bcftools')
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_arguments()
    if not os.path.exists(args.fastq_path):
        raise FileNotFoundError(f"FASTQ path {args.fastq_path} does not exist.")
    if args.targets and not os.path.exists(args.targets):
        raise FileNotFoundError(f"Targets file {args.targets} does not exist.")
    if args.ref_fasta and not os.path.exists(args.ref_fasta):
        raise FileNotFoundError(f"Reference FASTA file {args.ref_fasta} does not exist.")
    try:
        # Create output directory if it doesn't exist
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)
        # Check if the fastq_path is a directory or a file
        if os.path.isdir(args.fastq_path):
            fastq_files = [os.path.join(args.fastq_path, f) for f in os.listdir(args.fastq_path) if f.endswith('.fastq.gz')]
        else:
            fastq_files = [args.fastq_path]
        print("Pooling FASTQ files...")
        pooled_fastq = os.path.join(args.output_dir, os.path.basename(args.fastq_path) + '_pooled.fastq.gz')
        pool_fastq_files(fastq_files, pooled_fastq)
        
        print("Aligning reads...")
        # Run alignment for each FASTQ file and process for variants
        output_bam_path = os.path.join(args.output_dir, os.path.basename(args.fastq_path) + '.bam')
        align_reads(pooled_fastq, args.ref_fasta, output_bam_path)
        print("Calling variants...")
        target_positions=find_target_site_positions(args.ref_fasta, load_targeted_sites(args.targets))
        vcf_file=os.path.join(args.output_dir, os.path.basename(args.fastq_path) + '.vcf')
        call_variants(output_bam_path, args.ref_fasta, vcf_file, args.min_read_quality, args.min_alignment_quality, args.max_depth)
        analyze_variants(vcf_file,args.ref_fasta,target_positions,args.min_read_num, args.min_read_percentage)
        variant_summary(vcf_file)
        print("Analysis complete. Results saved in:", args.output_dir)
    except Exception as e:
        print(f"Error encountered: {e}")

if __name__ == "__main__":
    main()
