# ChIP-seq Analysis Pipeline Tutorial

## Overview

This tutorial provides step-by-step instructions for running the complete ChIP-seq analysis pipeline using Nextflow. The pipeline processes single-end ChIP-seq data from raw reads through quality control, alignment, peak calling, and annotation.

## Pipeline Overview

The pipeline performs the following steps:
1. Quality control with FastQC
2. Adapter trimming with Trimmomatic
3. Read alignment with Bowtie2
4. BAM processing (sorting, indexing, statistics)
5. BigWig coverage track generation
6. Sample correlation analysis
7. Peak calling with HOMER
8. Reproducible peak identification
9. Blacklist filtering
10. Peak annotation
11. Gene signal visualization across gene bodies
12. Motif enrichment analysis

## Prerequisites

### Software Requirements

The pipeline uses Singularity containers, so you need:
- Nextflow (version 20.0 or later)
- Singularity (for container execution)
- Access to compute cluster (optional but recommended)

### Input Data

The pipeline expects:
- Sample sheet CSV file with columns: `name`, `path`
- Reference genome FASTA file
- GTF annotation file
- Adapter sequences file
- ENCODE blacklist BED file
- UCSC genes BED file (TSS/TTS positions)

### Analysis Environment

For results visualization and interpretation:
- Conda environment (specified in `environment.yml`)
- Jupyter notebook with Python 3.10, pandas, matplotlib, seaborn
- Interactive analysis notebook (`analysis.ipynb`)

## Quick Start

### 1. Basic Test Run (Stub Mode)

Test the workflow logic without running actual computations:

```bash
nextflow run main.nf -stub-run -profile local
```

This will:
- Validate all process definitions
- Test channel logic
- Create empty stub output files
- Complete in seconds

### 2. Development Run (Subsampled Data)

Run with subsampled data for faster testing:

```bash
nextflow run main.nf -profile singularity,local
```

This uses `params.subsampled_samplesheet` by default, which points to smaller FASTQ files.

### 3. Production Run (Full Data)

For the complete analysis with full datasets:

```bash
nextflow run main.nf -profile singularity,cluster
```

## Configuration

### Sample Sheet Format

The pipeline expects a CSV file with the following format:

```csv
name,path
INPUT_rep1,/path/to/INPUT_rep1.fastq.gz
INPUT_rep2,/path/to/INPUT_rep2.fastq.gz
IP_rep1,/path/to/IP_rep1.fastq.gz
IP_rep2,/path/to/IP_rep2.fastq.gz
```

**Important naming conventions:**
- Samples must start with either `INPUT` or `IP`
- Must include replicate identifier: `rep1`, `rep2`, etc.
- Format: `{TYPE}_rep{N}` or `{TYPE}_rep{N}_subset`

### Switching Between Datasets

To switch from subsampled to full data, edit `main.nf` line 25:

```groovy
# Change from:
Channel.fromPath(params.subsampled_samplesheet)

# To:
Channel.fromPath(params.samplesheet)
```

### Key Parameters

Edit `nextflow.config` to customize:

```groovy
params {
    // Reference files
    genome = "/path/to/genome.fa"
    gtf = "/path/to/annotation.gtf"
    adapter_fa = "/path/to/adapters.fa"
    blacklist = "/path/to/blacklist.bed"
    
    // Sample sheets
    samplesheet = "$projectDir/full_samplesheet.csv"
    subsampled_samplesheet = "$projectDir/subsampled_samplesheet.csv"
    
    // Output directory
    outdir = "$projectDir/results/"
}
```

## Execution Profiles

### Local Profile

For running on a single machine:

```bash
nextflow run main.nf -profile singularity,local
```

- Uses local executor
- Single CPU per process
- Suitable for testing

### Cluster Profile

For running on SGE cluster:

```bash
nextflow run main.nf -profile singularity,cluster
```

- Uses SGE executor
- Allocates resources based on process labels:
  - `process_single`: 1 CPU
  - `process_low`: 4 CPUs
  - `process_medium`: 8 CPUs
  - `process_high`: 16 CPUs

## Pipeline Workflow Details

### Stage 1: Quality Control & Preprocessing

```
Raw FASTQ → FastQC (QC) → Trimmomatic (trim adapters) → Trimmed FASTQ
```

**Outputs:**
- `results/fastqc/` - FastQC HTML and ZIP reports
- `results/trimmomatic/` - Trimmed FASTQ files and logs

### Stage 2: Alignment

```
Trimmed FASTQ → Bowtie2 (align) → SAM → Samtools (convert) → BAM
                    ↑
            Bowtie2 Index (built once)
```

**Outputs:**
- `refs/bowtie2_index/` - Genome index files
- `results/bowtie2_align/` - BAM files and alignment logs

### Stage 3: BAM Processing

```
BAM → Samtools Sort → Sorted BAM → Samtools Index → Indexed BAM
                           ↓
                    Samtools Flagstat (statistics)
```

**Outputs:**
- `results/samtools_sort/` - Sorted BAM files
- `results/samtools_idx/` - BAM index files (.bai)
- `results/samtools_flagstat/` - Alignment statistics

### Stage 4: Coverage Tracks & Correlation

```
Indexed BAM → BamCoverage → BigWig files → MultiBigwigSummary → Matrix
                                                                    ↓
                                                            PlotCorrelation → Heatmap
```

**Outputs:**
- `results/bamcoverage/` - BigWig coverage tracks (.bw)
- `results/multibwsummary/` - Correlation matrix (NPZ format)
- `results/plotcorrelation/` - Correlation heatmap (PNG) and matrix (TSV)

### Stage 5: Peak Calling

```
Indexed BAM → makeTagDirectory → Tag Directory
                                       ↓
                     IP + INPUT paired by replicate
                                       ↓
                              findPeaks (HOMER)
                                       ↓
                              Peak TXT files → pos2bed → BED files
```

**Outputs:**
- `results/homer_tagdir/` - HOMER tag directories
- `results/homer_findpeaks/` - Peak files in HOMER format (.txt)
- `results/homer_pos2bed/` - Peak files in BED format (.bed)

### Stage 6: Reproducible Peaks & Filtering

```
rep1_peaks.bed ──┐
                 ├─→ Bedtools Intersect → Reproducible peaks
rep2_peaks.bed ──┘                              ↓
                                    Bedtools Remove (blacklist)
                                                ↓
                                    Filtered reproducible peaks
```

**Outputs:**
- `results/bedtools_intersect/` - Reproducible peaks (repr_peaks.bed)
- `results/bedtools_remove/` - Filtered peaks (repr_peaks_filtered.bed)

### Stage 7: Peak Annotation

```
Filtered peaks + Genome + GTF → annotatePeaks → Annotated peaks
```

**Outputs:**
- `results/homer_annotate/` - Annotated peaks (annotated_peaks.txt)

### Stage 8: Gene Signal Visualization

```
BigWig (IP only) + Genes BED → computeMatrix → Matrix → plotProfile → Signal plot
```

**Outputs:**
- `results/computematrix/` - Read coverage matrices for IP samples
- `results/plotprofile/` - Gene signal profile plots (PNG)

**Note:** Only IP samples are processed for gene signal visualization, as INPUT samples represent background noise.

### Stage 9: Motif Enrichment Analysis

```
Filtered peaks + Genome → findMotifsGenome → Motif results
```

**Outputs:**
- `results/homer_motifs/motifs/` - HOMER motif analysis results
  - `homerResults.html` - Main results page
  - `knownResults.html` - Known motif matches
  - Motif logos and statistics

### Stage 10: QC Aggregation

```
FastQC reports ──┐
Trimmomatic logs ├─→ MultiQC → Comprehensive QC report
Bowtie2 logs    ─┤
Flagstat files ──┘
```

**Outputs:**
- `results/multiqc/` - MultiQC HTML report and data directory

## Monitoring Pipeline Execution

### View Progress

Nextflow provides real-time progress monitoring:

```bash
nextflow run main.nf -profile singularity,cluster
```

Output shows:
```
executor >  sge (12)
[ab/cd1234] process > FASTQC (1)           [100%] 4 of 4 ✔
[ef/gh5678] process > TRIM (2)             [100%] 4 of 4 ✔
[ij/kl9012] process > BOWTIE2_BUILD        [100%] 1 of 1 ✔
...
```

### Resume Failed Runs

If the pipeline fails, resume from the last successful step:

```bash
nextflow run main.nf -profile singularity,cluster -resume
```

Nextflow uses caching to skip completed processes.

## Output Structure

```
project-3-Luckyfruit88/
├── main.nf                      # Main workflow file
├── nextflow.config              # Configuration
├── modules/                     # Process definitions
├── results/                     # Pipeline outputs
│   ├── fastqc/
│   ├── trimmomatic/
│   ├── bowtie2_align/
│   ├── samtools_sort/
│   ├── samtools_idx/
│   ├── samtools_flagstat/
│   ├── bamcoverage/
│   ├── multibwsummary/
│   ├── plotcorrelation/
│   ├── homer_tagdir/
│   ├── homer_findpeaks/
│   ├── homer_pos2bed/
│   ├── bedtools_intersect/
│   ├── bedtools_remove/
│   ├── homer_annotate/
│   └── multiqc/
├── work/                        # Temporary working directory
└── .nextflow/                   # Nextflow metadata
```

## Troubleshooting

### Common Issues

**1. Container pull errors**
```bash
# Pre-pull containers
singularity pull docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
```

**2. Memory errors**
- Increase memory allocation in `nextflow.config`
- Use subsampled data for testing

**3. File path errors**
- Use absolute paths in configuration
- Verify file permissions

**4. Sample pairing issues**
- Check sample naming conventions
- Ensure `IP_rep1` pairs with `INPUT_rep1`

### Debugging

**Check process logs:**
```bash
# Find failed process directory
ls -lt work/
cd work/ab/cd1234...

# View log files
cat .command.log
cat .command.err
```

**Validate sample sheet:**
```bash
# Test CSV parsing
nextflow run main.nf -stub-run -profile local
```

## Advanced Usage

### Custom Resource Allocation

Edit process labels in `nextflow.config`:

```groovy
process {
    withLabel: process_high {
        cpus = 32
        memory = '64 GB'
    }
}
```

### Custom Correlation Method

Change correlation method in `modules/deeptools_plotcorrelation/main.nf`:

```groovy
plotCorrelation -in ${matrix} \
    --corMethod pearson \  // Change from spearman
    --whatToPlot heatmap \
    --plotFile correlation_plot.png
```

### Custom Peak Calling Parameters

Modify HOMER findPeaks in `modules/homer_findpeaks/main.nf`:

```groovy
findPeaks ${ip_tagdir} \
    -style factor \
    -i ${input_tagdir} \
    -F 3 \              // Fold enrichment threshold
    -P 0.001 \          // P-value threshold
    -o ${rep}_peaks.txt
```

## Performance Tips

1. **Use cluster profile** for production runs
2. **Enable resume** with `-resume` flag
3. **Monitor work directory** size and clean periodically
4. **Use subsampled data** during development
5. **Adjust process labels** based on data size

## Expected Runtime

### Subsampled Data (local profile)
- Total: ~30-60 minutes
- Alignment: ~5-10 minutes per sample
- Peak calling: ~2-5 minutes per replicate

### Full Data (cluster profile)
- Total: ~4-8 hours
- Alignment: ~30-60 minutes per sample
- Peak calling: ~15-30 minutes per replicate

## Results Interpretation

### Key Outputs to Review

1. **MultiQC Report** (`results/multiqc/multiqc_report.html`)
   - Overall data quality
   - Alignment rates
   - Sample statistics

2. **Correlation Heatmap** (`results/plotcorrelation/correlation_plot.png`)
   - Sample similarity
   - Expected: IP samples cluster together, INPUT samples cluster together

3. **Annotated Peaks** (`results/homer_annotate/annotated_peaks.txt`)
   - Genomic location of peaks
   - Nearest genes
   - Peak characteristics

4. **Final Peak Set** (`results/bedtools_remove/repr_peaks_filtered.bed`)
   - Reproducible peaks between replicates
   - Blacklist-filtered
   - Ready for downstream analysis

5. **Gene Signal Profiles** (`results/plotprofile/IP_rep*_signal_coverage.png`)
   - Signal distribution across gene bodies
   - Binding patterns at TSS, gene body, and TTS
   - 2kb flanking regions shown

6. **Motif Enrichment** (`results/homer_motifs/motifs/homerResults.html`)
   - Enriched sequence motifs
   - Known motif matches
   - Potential co-factors

### Using the Analysis Notebook

The pipeline includes a comprehensive Jupyter notebook for results visualization:

```bash
# Set up conda environment
conda env create -f environment.yml
conda activate chipseq-analysis

# Launch Jupyter notebook
jupyter notebook analysis.ipynb
```

**Notebook sections:**
- Quality control summaries with MultiQC integration
- Correlation analysis interpretation (Spearman justification)
- Peak calling statistics and visualizations
- Reproducible peak strategy explanation
- Peak annotation analysis
- Gene signal profile displays
- Motif enrichment results
- Summary and conclusions

## Citation

If you use this pipeline, please cite the following tools:

- Nextflow: Di Tommaso et al., Nature Biotechnology 2017
- FastQC: Andrews S., Babraham Bioinformatics 2010
- Trimmomatic: Bolger et al., Bioinformatics 2014
- Bowtie2: Langmead and Salzberg, Nature Methods 2012
- SAMtools: Li et al., Bioinformatics 2009
- deepTools: Ramírez et al., Nucleic Acids Research 2016
- HOMER: Heinz et al., Molecular Cell 2010
- BEDTools: Quinlan and Hall, Bioinformatics 2010

## Support

For issues with the pipeline:
1. Check the troubleshooting section above
2. Review process logs in `work/` directory
3. Consult individual tool documentation
4. Verify configuration parameters

## License

This pipeline is provided for educational purposes as part of BF528 coursework.
