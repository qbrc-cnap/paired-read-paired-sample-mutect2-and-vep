import "bwa_align.wdl" as bwa_align
import "gatk_tools.wdl" as gatk_tools
import "vep.wdl" as vep

workflow SingleSampleHaplotypecallerWorkflow {
    # This workflow operates on a single "sample" and
    # assumes that all the sequences reads are contained in
    # 2 FASTQ files, paired-end sequencing.
    #
    # Output are 
    
    File r1_fastq
    File r2_fastq
    Boolean use_dedup
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_bwt
    File ref_sa
    File ref_amb
    File ref_ann
    File ref_pac
    File dbsnp
    File dbsnp_index
    File known_indels
    File known_indels_index

    File vep_cache_tar
    String vep_species

    Array[String] contigs

    # Extract the samplename from the fastq filename
    String sample_name = basename(r1_fastq, "_R1.fastq.gz")

    # Perform the alignment, which output a sorted BAM and BAM index
    call bwa_align.perform_align as alignment {
        input:
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq,
            sample_name = sample_name,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            ref_bwt = ref_bwt,
            ref_amb = ref_amb,
            ref_ann = ref_ann,
            ref_pac = ref_pac,
            ref_sa = ref_sa
    }

    # Perform alignment QC
    call gatk_tools.run_alignment_metrics as aln_metrics {
        input:
            input_bam = alignment.sorted_bam,
            input_bam_index = alignment.sorted_bam_index,
            sample_name = sample_name,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict
    }

    # Identify points of recalibration of the BAM
    call gatk_tools.base_recalibrator as base_recal {
        input:
            input_bam = alignment.sorted_bam,
            input_bam_index = alignment.sorted_bam_index,
            sample_name = sample_name,
            dbsnp = dbsnp,
            dbsnp_index = dbsnp_index,
            known_indels = known_indels,
            known_indels_index = known_indels_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict
    }

    # Applys the recalibration of reads to the BAM
    call gatk_tools.apply_recalibration as apply_recal {
        input:
            input_bam = alignment.sorted_bam,
            input_bam_index = alignment.sorted_bam_index,
            recalibration_report = base_recal.recalibration_report,
            sample_name = sample_name,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict
    }

    # Deduplicate the BAM file
    call gatk_tools.deduplicate_bam as dedup_bam {
        input:
            input_bam = apply_recal.recalibrated_bam,
            input_bam_index = apply_recal.recalibrated_bam_index,
            sample_name = sample_name
    }

    # Scatters over the contig intervals
    scatter (scatter_interval in contigs) {
        # Identifies variants with GATK's HaplotypeCaller
        call gatk_tools.haplotypecaller as haplotypecaller {
            input:
                input_bam = apply_recal.recalibrated_bam,
                input_bam_index = apply_recal.recalibrated_bam_index,
                input_dedup_bam = dedup_bam.sorted_bam,
                input_dedup_bam_index = dedup_bam.sorted_bam_index,
                use_dedup = use_dedup,
                interval = scatter_interval,
                sample_name = sample_name,
                ref_dict = ref_dict,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
        }
    }
    
    # Merges the scattered VCFs together
    call gatk_tools.merge_vcf as merge_vcf {
        input:
            input_vcfs = haplotypecaller.output_vcf,
            sample_name = sample_name,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index
    }

    # Annotated the VCF with VEP
    call vep.vep_annotate as vep_annotate {
        input:
            input_vcf = merge_vcf.output_vcf,
            sample_name = sample_name,
            species = vep_species,
            vep_cache_tar = vep_cache_tar
    }

    output {
        File vcf = merge_vcf.output_vcf
        File annotated_vcf = vep_annotate.output_vcf
        File annotated_vcf_stats = vep_annotate.vcf_stats
        File deduplication_metrics = dedup_bam.deduplication_metrics
        File alignment_metrics = aln_metrics.alignment_metrics
    }
}