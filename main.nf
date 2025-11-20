// Include your modules here
include { FASTQC } from './modules/fastqc/main.nf'
include { TRIM } from './modules/trimmomatic/main.nf'
include { BOWTIE2_BUILD } from './modules/bowtie2_build/main.nf'
include { BOWTIE2_ALIGN } from './modules/bowtie2_align/main.nf'
include { SAMTOOLS_SORT } from './modules/samtools_sort/main.nf'
include { SAMTOOLS_IDX } from './modules/samtools_idx/main.nf'
include { SAMTOOLS_FLAGSTAT } from './modules/samtools_flagstat/main.nf'
include { MULTIQC } from './modules/multiqc/main.nf'
include { BAMCOVERAGE } from './modules/deeptools_bamcoverage/main.nf'
include { MULTIBWSUMMARY } from './modules/deeptools_multibwsummary/main.nf'
include { PLOTCORRELATION } from './modules/deeptools_plotcorrelation/main.nf'
include { TAGDIR } from './modules/homer_maketagdir/main.nf'
include { FINDPEAKS } from './modules/homer_findpeaks/main.nf'
include { POS2BED } from './modules/homer_pos2bed/main.nf'
include { BEDTOOLS_INTERSECT } from './modules/bedtools_intersect/main.nf'
include { BEDTOOLS_REMOVE } from './modules/bedtools_remove/main.nf'
include { ANNOTATE } from './modules/homer_annotatepeaks/main.nf'
include { COMPUTEMATRIX } from './modules/deeptools_computematrix/main.nf'
include { PLOTPROFILE } from './modules/deeptools_plotprofile/main.nf'
include { FIND_MOTIFS_GENOME } from './modules/homer_findmotifsgenome/main.nf'

workflow {

    
    //Here we construct the initial channels we need
    
    Channel.fromPath(params.subsampled_samplesheet)
    | splitCsv( header: true )
    | map{ row -> tuple(row.name, file(row.path)) }
    | set { read_ch }

    // Quality control on raw reads
    FASTQC(read_ch)
    
    // Trim adapters and low quality bases
    TRIM(read_ch)
    
    // Build bowtie2 index for the reference genome
    BOWTIE2_BUILD(file(params.genome))
    
    // Align trimmed reads to the reference genome
    BOWTIE2_ALIGN(TRIM.out.trimmed_reads, BOWTIE2_BUILD.out.index)
    
    // Sort BAM files
    SAMTOOLS_SORT(BOWTIE2_ALIGN.out.bam)
    
    // Index sorted BAM files
    SAMTOOLS_IDX(SAMTOOLS_SORT.out.sorted_bam)
    
    // Calculate alignment statistics
    SAMTOOLS_FLAGSTAT(SAMTOOLS_SORT.out.sorted_bam)
    
    // Generate BigWig coverage files
    BAMCOVERAGE(SAMTOOLS_IDX.out.indexed_bam)
    
    // Collect all QC outputs for MultiQC
    FASTQC.out.zip
        .map { sample_id, zip -> zip }
        .mix(
            TRIM.out.log.map { sample_id, log -> log },
            BOWTIE2_ALIGN.out.log.map { sample_id, log -> log },
            SAMTOOLS_FLAGSTAT.out.flagstat.map { sample_id, flagstat -> flagstat }
        )
        .collect()
        .set { multiqc_ch }
    
    // Run MultiQC
    MULTIQC(multiqc_ch)
    
    // Collect bigWig files and compute correlation matrix
    BAMCOVERAGE.out.bigwig
        .map { sample_id, bw -> bw }
        .collect()
        .set { bw_ch }
    
    MULTIBWSUMMARY(bw_ch)
    PLOTCORRELATION(MULTIBWSUMMARY.out.matrix)
    
    // Create HOMER tag directories from indexed BAM files
    TAGDIR(SAMTOOLS_IDX.out.indexed_bam)
    
    // Pair IP and INPUT samples for peak calling
    // Extract replicate number and sample type from sample_id
    TAGDIR.out.tagdir
        .map { sample_id, tagdir -> 
            def parts = sample_id.tokenize('_')
            def type = parts[0]  // INPUT or IP
            def rep = parts[1]   // rep1 or rep2 (may have _subset after)
            // Extract just the replicate identifier (e.g., rep1, rep2)
            tuple(type, rep, sample_id, tagdir)
        }
        .branch {
            ip: it[0] == 'IP'
            input: it[0] == 'INPUT'
        }
        .set { tagdir_branched }
    
    // Combine IP and INPUT by replicate
    tagdir_branched.ip
        .map { type, rep, sample_id, tagdir -> tuple(rep, tagdir) }
        .set { ip_tagdir }
    
    tagdir_branched.input
        .map { type, rep, sample_id, tagdir -> tuple(rep, tagdir) }
        .set { input_tagdir }
    
    ip_tagdir
        .join(input_tagdir)
        .set { paired_tagdir }
    
    // Call peaks using HOMER
    FINDPEAKS(paired_tagdir)
    
    // Convert HOMER peaks to BED format
    POS2BED(FINDPEAKS.out.peaks)
    
    // Find reproducible peaks between replicates
    POS2BED.out.bed
        .map { rep, bed -> bed }
        .collect()
        .set { bed_files }
    
    BEDTOOLS_INTERSECT(bed_files)
    
    // Remove blacklisted regions
    BEDTOOLS_REMOVE(BEDTOOLS_INTERSECT.out.peaks, file(params.blacklist))
    
    // Annotate filtered peaks
    ANNOTATE(BEDTOOLS_REMOVE.out.filtered_peaks, file(params.genome), file(params.gtf))
    
    // Filter IP samples only for gene signal visualization
    BAMCOVERAGE.out.bigwig
        .filter { sample_id, bw -> sample_id.startsWith('IP') }
        .set { ip_bigwig }
    
    // Compute matrix for gene signal across IP samples
    COMPUTEMATRIX(ip_bigwig, file(params.ucsc_genes))
    
    // Plot profile for gene signal visualization
    PLOTPROFILE(COMPUTEMATRIX.out.matrix)
    
    // Find enriched motifs in reproducible filtered peaks
    FIND_MOTIFS_GENOME(BEDTOOLS_REMOVE.out.filtered_peaks, file(params.genome))


}