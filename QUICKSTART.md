# Quick Start Guide

## Running the Pipeline

### Test Mode (Fast)
```bash
nextflow run main.nf -stub-run -profile local
```
Validates workflow logic without running actual analysis.

### Development Mode (Subsampled Data)
```bash
nextflow run main.nf -profile singularity,local
```
Runs with small subsampled files for quick testing.

### Production Mode (Full Data)
```bash
nextflow run main.nf -profile singularity,cluster
```
Full analysis on complete dataset using cluster resources.

## Common Commands

### Resume After Failure
```bash
nextflow run main.nf -profile singularity,cluster -resume
```

### View Pipeline DAG
```bash
nextflow run main.nf -with-dag flowchart.png
```

### Generate Report
```bash
nextflow run main.nf -profile singularity,cluster -with-report report.html
```

## Switching to Full Data

Edit `main.nf` line 25:
```groovy
Channel.fromPath(params.samplesheet)  // instead of params.subsampled_samplesheet
```

## Key Output Files

| File | Location | Description |
|------|----------|-------------|
| MultiQC Report | `results/multiqc/multiqc_report.html` | Comprehensive QC summary |
| Correlation Plot | `results/plotcorrelation/correlation_plot.png` | Sample similarity heatmap |
| Final Peaks | `results/bedtools_remove/repr_peaks_filtered.bed` | Reproducible, filtered peaks |
| Annotated Peaks | `results/homer_annotate/annotated_peaks.txt` | Peaks with gene annotations |
| Gene Signal Plots | `results/plotprofile/IP_rep*_signal_coverage.png` | Signal across gene bodies (IP samples) |
| Motif Results | `results/homer_motifs/motifs/homerResults.html` | Enriched motifs in peaks |
| BigWig Files | `results/bamcoverage/*.bw` | Coverage tracks for visualization |
| Analysis Notebook | `analysis.ipynb` | Jupyter notebook with all figures and discussion |

## Troubleshooting

### Check process logs
```bash
# Go to work directory
cd work/
# Find latest directory
ls -lt | head
# View logs
cd <hash>/<hash>/
cat .command.log
cat .command.err
```

### Verify sample sheet
```bash
head subsampled_samplesheet.csv
head full_samplesheet.csv
```

### Clean work directory
```bash
rm -rf work/
rm -rf .nextflow/
```

## Sample Naming Requirements

Samples must follow this pattern:
- `INPUT_rep1` - Control replicate 1
- `INPUT_rep2` - Control replicate 2  
- `IP_rep1` - IP sample replicate 1
- `IP_rep2` - IP sample replicate 2

Subsampled files may have `_subset` suffix:
- `INPUT_rep1_subset`
- `IP_rep1_subset`

## Pipeline Stages

1. **QC & Preprocessing** - FastQC, Trimmomatic
2. **Alignment** - Bowtie2
3. **BAM Processing** - Sort, Index, Flagstat
4. **Coverage** - bamCoverage, Correlation
5. **Peak Calling** - HOMER (makeTagDirectory, findPeaks)
6. **Peak Processing** - pos2bed, Intersect, Filter
7. **Annotation** - annotatePeaks
8. **Gene Signal** - computeMatrix, plotProfile (IP only)
9. **Motif Analysis** - findMotifsGenome
10. **QC Report** - MultiQC

## Expected Runtime

| Dataset | Profile | Time |
|---------|---------|------|
| Stub run | local | < 1 min |
| Subsampled | local | 30-60 min |
| Full | cluster | 4-8 hours |

## Analysis Notebook

### Set Up Environment
```bash
conda env create -f environment.yml
conda activate chipseq-analysis
jupyter notebook analysis.ipynb
```

### Notebook Features
- QC summaries and correlation interpretation
- Peak statistics and reproducibility justification
- Gene signal profile visualization
- Motif enrichment results
- Complete analysis workflow with figures

## Getting Help

1. Read `TUTORIAL.md` for detailed instructions
2. Check `TUTORIAL_CN.md` for Chinese instructions
3. Review `analysis.ipynb` for results interpretation
4. Review process logs in `work/` directory
5. Verify configuration in `nextflow.config`
