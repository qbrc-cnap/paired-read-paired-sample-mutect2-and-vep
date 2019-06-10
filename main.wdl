import "paired_sample_mutect2.wdl" as paired_sample_mutect2
import "fastqc.wdl" as fastqc
import "report.wdl" as reporting

workflow PairedMatchedMutect2AndVepWorkflow {
    # This workflow is a 'super' workflow that parallelizes
    # HaplotypeCaller and VEP analysis over multiple samples.

    # Input files
    Array[File] r1_files
    Array[File] r2_files
    File match_annotations
    Array[Array[File]] matched_samples = read_tsv(match_samples_to_tsv.output_annotations)
    Boolean use_dedup

    #Array[Pair[File, File]] fastq_pairs = zip(r1_files, r2_files)
    
    # Reference files
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_bwt
    File ref_sa
    File ref_amb
    File ref_ann
    File ref_pac
    File ref_exon_intervals

    # Inputs for GATK BQSR
    File dbsnp
    File dbsnp_index
    File known_indels
    File known_indels_index

    # Inputs for HaplotypeCaller and Mutect2 scatter
    File contig_list
    Array[String] contigs = read_lines(contig_list)

    # Inputs for GATK Mutect2
    File gnomad
    File gnomad_index

    # Inputs for VEP
    String vep_species
    File vep_cache_tar
    
    # Other
    String output_zip_name
    String genome
    String git_repo_url
    String git_commit_hash

    call match_samples as match_samples_to_tsv {
        input:
            r1_fastq = r1_files,
            r2_fastq = r2_files,
            match_annotations = match_annotations
    }

    scatter(item in matched_samples){
        call fastqc.run_fastqc as fastqc_for_tumor_read1 {
            input:
                fastq = item[0]
        }
        
        call fastqc.run_fastqc as fastqc_for_tumor_read2 {
            input:
                fastq = item[1]
        }

        call fastqc.run_fastqc as fastqc_for_normal_read1 {
            input:
                fastq = item[2]
        }
        
        call fastqc.run_fastqc as fastqc_for_normal_read2 {
            input:
                fastq = item[3]
        }

        call paired_sample_mutect2.PairedSampleMutect2Workflow as paired_sample_process {
            input:
                tumor_r1_fastq = item[0],
                tumor_r2_fastq = item[1],
                normal_r1_fastq = item[2],
                normal_r2_fastq = item[3],
                use_dedup = use_dedup,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                ref_pac = ref_pac,
                ref_bwt = ref_bwt,
                ref_sa = ref_sa,
                ref_amb = ref_amb,
                ref_ann = ref_ann,
                ref_exon_intervals = ref_exon_intervals,
                dbsnp = dbsnp,
                dbsnp_index = dbsnp_index,
                known_indels = known_indels,
                known_indels_index = known_indels_index,
                gnomad = gnomad,
                gnomad_index = gnomad_index,
                vep_cache_tar = vep_cache_tar,
                vep_species = vep_species,
                contigs = contigs
        }
    }

    call reporting.create_multi_qc as multiqc {
        input:
            tumor_alignment_metrics = paired_sample_process.tumor_alignment_metrics,
            tumor_dedup_metrics = paired_sample_process.tumor_deduplication_metrics,
            normal_alignment_metrics = paired_sample_process.normal_alignment_metrics,
            normal_dedup_metrics = paired_sample_process.normal_deduplication_metrics,
            tumor_r1_fastqc_zips = fastqc_for_tumor_read1.fastqc_zip,
            tumor_r2_fastqc_zips = fastqc_for_tumor_read2.fastqc_zip,
            normal_r1_fastqc_zips = fastqc_for_normal_read1.fastqc_zip,
            normal_r2_fastqc_zips = fastqc_for_normal_read2.fastqc_zip,
            concordance_metrics = paired_sample_process.concordance_metrics,
            contamination_metrics = paired_sample_process.contamination_metrics
    }

    call reporting.generate_report as generate_report {
        input:
            r1_files = r1_files,
            r2_files = r2_files,
            genome = genome,
            dedup_bool = use_dedup,
            git_commit_hash = git_commit_hash,
            git_repo_url = git_repo_url,
    }

    call zip_results {
        input:
            zip_name = output_zip_name,
            vcf_files = paired_sample_process.vcf,
            tsv_files = paired_sample_process.annotated_vcf,
            vep_stats_files = paired_sample_process.annotated_vcf_stats,
            multiqc_results = multiqc.report,
            analysis_report = generate_report.report
    }

    output {
        File zip_out = zip_results.zip_out
    }

    meta {
        workflow_title: "Somatic exome variant calling"
        workflow_short_description: "A pipline for variant calling and variant annotation from tumor / normal matched exome DNASeq libraries"
        workflow_long_description: "Use this workflow for aligning matched tumor-normal paired-end Illumina NGS libraries with BWA, optional deduplication, base quality score recalibration with GATK, variant calling with GATK Mutect2, and variant annotation with Ensembl's VEP."
    }
}

task match_samples {
    Array[String] r1_fastq
    Array[String] r2_fastq
    File match_annotations

    # runtime parameters
    Int disk_size = 20

    command {
        python3 \
            /usr/local/bin/match_annotations.py \
            -r1 ${sep=" " r1_fastq} \
            -r2 ${sep=" " r2_fastq} \
            -annot ${match_annotations} \
        > matched_annotations.tsv;
    }

    output {
        File output_annotations = "matched_annotations.tsv"
    }

    runtime {
        docker: "docker.io/hsphqbrc/gatk-mutect2-workflow-tools:1.0"
        cpu: 1
        memory: "2 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task zip_results {
    String zip_name
    Array[File] vcf_files
    Array[File] tsv_files
    Array[File] vep_stats_files
    File multiqc_results
    File analysis_report

    # runtime parameters
    Int disk_size = 250

    command {
        mkdir output
        mkdir output/VCFs
        mkdir output/TSVs
        mkdir output/variant_stats
        mv ${multiqc_results} output
        mv ${analysis_report} output
        mv -t output/VCFs ${sep=" " vcf_files}
        mv -t output/TSVs ${sep=" " tsv_files}
        mv -t output/variant_stats ${sep=" " vep_stats_files}
        zip -r "${zip_name}.report_and_output.zip" output
    }

    output {
        File zip_out = "${zip_name}.report_and_output.zip"
    }

    runtime {
        docker: "docker.io/hsphqbrc/gatk-mutect2-workflow-tools:1.0"
        cpu: 2
        memory: "6 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}