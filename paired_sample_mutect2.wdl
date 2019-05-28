import "bwa_align.wdl" as bwa_align
import "gatk_tools.wdl" as gatk_tools
import "vep.wdl" as vep

workflow PairedSampleMutect2Workflow {
    # This workflow operates on a single "sample" and
    # assumes that all the sequences reads are contained in
    # 2 FASTQ files, paired-end sequencing.
    #
    # Output are 
    
    File tumor_r1_fastq
    File tumor_r2_fastq
    File normal_r1_fastq
    File normal_r2_fastq
    Boolean use_dedup
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_bwt
    File ref_sa
    File ref_amb
    File ref_ann
    File ref_pac
    File ref_exon_intervals
    File dbsnp
    File dbsnp_index
    File known_indels
    File known_indels_index
    File gnomad
    File gnomad_index

    File vep_cache_tar
    String vep_species

    Array[String] contigs

    # Extract the samplename from the fastq filename
    String tumor_sample_name = basename(tumor_r1_fastq, "_R1.fastq.gz")
    String normal_sample_name = basename(normal_r1_fastq, "_R1.fastq.gz")

    # Perform the alignment, which output a sorted BAM and BAM index
    call bwa_align.perform_align as tumor_alignment {
        input:
            r1_fastq = tumor_r1_fastq,
            r2_fastq = tumor_r2_fastq,
            sample_name = tumor_sample_name,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            ref_bwt = ref_bwt,
            ref_amb = ref_amb,
            ref_ann = ref_ann,
            ref_pac = ref_pac,
            ref_sa = ref_sa
    }

    call bwa_align.perform_align as normal_alignment {
        input:
            r1_fastq = normal_r1_fastq,
            r2_fastq = normal_r2_fastq,
            sample_name = normal_sample_name,
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
    call gatk_tools.run_alignment_metrics as tumor_aln_metrics {
        input:
            input_bam = tumor_alignment.sorted_bam,
            input_bam_index = tumor_alignment.sorted_bam_index,
            sample_name = tumor_sample_name,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict
    }

    call gatk_tools.run_alignment_metrics as normal_aln_metrics {
        input:
            input_bam = normal_alignment.sorted_bam,
            input_bam_index = normal_alignment.sorted_bam_index,
            sample_name = normal_sample_name,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict
    }

    # Identify points of recalibration of the BAM
    call gatk_tools.base_recalibrator as tumor_base_recal {
        input:
            input_bam = tumor_alignment.sorted_bam,
            input_bam_index = tumor_alignment.sorted_bam_index,
            sample_name = tumor_sample_name,
            dbsnp = dbsnp,
            dbsnp_index = dbsnp_index,
            known_indels = known_indels,
            known_indels_index = known_indels_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict
    }

    call gatk_tools.base_recalibrator as normal_base_recal {
        input:
            input_bam = normal_alignment.sorted_bam,
            input_bam_index = normal_alignment.sorted_bam_index,
            sample_name = normal_sample_name,
            dbsnp = dbsnp,
            dbsnp_index = dbsnp_index,
            known_indels = known_indels,
            known_indels_index = known_indels_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict
    }

    # Applys the recalibration of reads to the BAM
    call gatk_tools.apply_recalibration as tumor_apply_recal {
        input:
            input_bam = tumor_alignment.sorted_bam,
            input_bam_index = tumor_alignment.sorted_bam_index,
            recalibration_report = tumor_base_recal.recalibration_report,
            sample_name = tumor_sample_name,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict
    }

    call gatk_tools.apply_recalibration as normal_apply_recal {
        input:
            input_bam = normal_alignment.sorted_bam,
            input_bam_index = normal_alignment.sorted_bam_index,
            recalibration_report = normal_base_recal.recalibration_report,
            sample_name = normal_sample_name,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict
    }

    # Deduplicate the BAM file
    call gatk_tools.deduplicate_bam as tumor_deduplicate_bam {
        input:
            input_bam = tumor_apply_recal.recalibrated_bam,
            input_bam_index = tumor_apply_recal.recalibrated_bam_index,
            sample_name = tumor_sample_name
    }

    call gatk_tools.deduplicate_bam as normal_deduplicate_bam {
        input:
            input_bam = normal_apply_recal.recalibrated_bam,
            input_bam_index = normal_apply_recal.recalibrated_bam_index,
            sample_name = normal_sample_name
    }

    call gatk_tools.conpair_pileup as conpair_tumor_pileup {
        input:
            input_bam = tumor_apply_recal.recalibrated_bam,
            input_bam_index = tumor_apply_recal.recalibrated_bam_index,
            input_dedup_bam = tumor_deduplicate_bam.sorted_bam,
            input_dedup_bam_index = tumor_deduplicate_bam.sorted_bam_index,
            use_dedup = use_dedup,
            sample_name = tumor_sample_name,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index
    }

    call gatk_tools.conpair_pileup as conpair_normal_pileup {
        input:
            input_bam = normal_apply_recal.recalibrated_bam,
            input_bam_index = normal_apply_recal.recalibrated_bam_index,
            input_dedup_bam = normal_deduplicate_bam.sorted_bam,
            input_dedup_bam_index = normal_deduplicate_bam.sorted_bam_index,
            use_dedup = use_dedup,
            sample_name = normal_sample_name,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index
    }

    call gatk_tools.conpair_concordance as conpair_concordance {
        input:
            tumor_pileup = conpair_tumor_pileup.output_pileup,
            normal_pileup = conpair_normal_pileup.output_pileup,
            sample_name = tumor_sample_name
    }

    call gatk_tools.conpair_contamination as conpair_contamination {
        input:
            tumor_pileup = conpair_tumor_pileup.output_pileup,
            normal_pileup = conpair_normal_pileup.output_pileup,
            sample_name = tumor_sample_name
    }

    call gatk_tools.coverage_metrics as tumor_coverage {
        input:
            input_bam = tumor_apply_recal.recalibrated_bam,
            input_bam_index = tumor_apply_recal.recalibrated_bam_index,
            input_dedup_bam = tumor_deduplicate_bam.sorted_bam,
            input_dedup_bam_index = tumor_deduplicate_bam.sorted_bam_index,
            use_dedup = use_dedup,
            sample_name = tumor_sample_name,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_exon_intervals = ref_exon_intervals
    }

    call gatk_tools.coverage_metrics as normal_coverage {
        input:
            input_bam = normal_apply_recal.recalibrated_bam,
            input_bam_index = normal_apply_recal.recalibrated_bam_index,
            input_dedup_bam = normal_deduplicate_bam.sorted_bam,
            input_dedup_bam_index = normal_deduplicate_bam.sorted_bam_index,
            use_dedup = use_dedup,
            sample_name = normal_sample_name,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_exon_intervals = ref_exon_intervals
    }

    # Scatters over the contig intervals
    scatter (scatter_interval in contigs) {
        # Identifies variants with GATK's HaplotypeCaller
        call gatk_tools.mutect as mutect {
            input:
                tumor_bam = tumor_apply_recal.recalibrated_bam,
                tumor_bam_index = tumor_apply_recal.recalibrated_bam_index,
                tumor_dedup_bam = tumor_deduplicate_bam.sorted_bam,
                tumor_dedup_bam_index = tumor_deduplicate_bam.sorted_bam_index,
                normal_bam = normal_apply_recal.recalibrated_bam,
                normal_bam_index = normal_apply_recal.recalibrated_bam_index,
                normal_dedup_bam = normal_deduplicate_bam.sorted_bam,
                normal_dedup_bam_index = normal_deduplicate_bam.sorted_bam_index,
                use_dedup = use_dedup,
                interval = scatter_interval,
                tumor_sample_name = tumor_sample_name,
                normal_sample_name = normal_sample_name,
                ref_dict = ref_dict,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                gnomad = gnomad,
                gnomad_index = gnomad_index
        }
    }
    
    # Merges the scattered VCFs together
    call gatk_tools.merge_vcf as merge_vcf {
        input:
            input_vcfs = mutect.output_vcf,
            sample_name = tumor_sample_name,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index
    }

    # Annotated the VCF with VEP
    call vep.vep_annotate as vep_annotate {
        input:
            input_vcf = merge_vcf.output_vcf,
            sample_name = tumor_sample_name,
            species = vep_species,
            vep_cache_tar = vep_cache_tar,
            is_vcf_empty = merge_vcf.is_vcf_empty
    }

    output {
        File vcf = merge_vcf.output_vcf
        File annotated_vcf = vep_annotate.output_vcf
        File annotated_vcf_stats = vep_annotate.vcf_stats
        File tumor_deduplication_metrics = tumor_deduplicate_bam.deduplication_metrics
        File normal_deduplication_metrics = normal_deduplicate_bam.deduplication_metrics
        File tumor_alignment_metrics = tumor_aln_metrics.alignment_metrics
        File normal_alignment_metrics = normal_aln_metrics.alignment_metrics
        File concordance_metrics = conpair_concordance.concordance_metrics
        File contamination_metrics = conpair_contamination.contamination_metrics
        File tumor_coverage_metrics = tumor_coverage.coverage_metrics
        File normal_coverage_metrics = normal_coverage.coverage_metrics
    }
}