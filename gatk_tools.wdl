workflow gatk_test {
    File input_bam
    File input_bam_index
    String sample_name
    Boolean use_dedup
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File dbsnp
    File dbsnp_index
    File known_indels
    File known_indels_index
    File gnomad
    File gnomad_index
    Array[String] contigs

    call base_recalibrator {
        input:
            input_bam = input_bam,
            input_bam_index = input_bam_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            dbsnp = dbsnp,
            dbsnp_index = dbsnp_index,
            known_indels = known_indels,
            known_indels_index = known_indels_index
    }

    call run_alignment_metrics {
        input:
            input_bam = input_bam,
            input_bam_index = input_bam_index,
            sample_name = sample_name,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict
    }

    call apply_recalibration {
        input:
            input_bam = input_bam,
            input_bam_index = input_bam_index,
            recalibration_report = base_recalibrator.recalibration_report,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict
    }

    call deduplicate_bam {
        input:
            input_bam = input_bam,
            input_bam_index = input_bam_index,
            sample_name = sample_name
    }

    scatter (scatter_interval in contigs) {
        # Identifies variants with GATK's HaplotypeCaller
        call haplotypecaller {
            input:
                input_bam = apply_recalibration.recalibrated_bam,
                input_bam_index = apply_recalibration.recalibrated_bam_index,
                input_dedup_bam = deduplicate_bam.sorted_bam,
                input_dedup_bam_index = deduplicate_bam.sorted_bam_index,
                use_dedup = use_dedup,
                interval = scatter_interval,
                sample_name = sample_name,
                ref_dict = ref_dict,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
        }
    }
    
    # Merges the scattered VCFs together
    call merge_vcf {
        input:
            input_vcfs = haplotypecaller.output_vcf,
            sample_name = sample_name,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index
    }

    output {
        File vcf = merge_vcf.output_vcf
    }
}

task run_alignment_metrics {
    File input_bam
    File input_bam_index
    String sample_name

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    # Runtime parameters
    Int disk_size = 100

    command {
        mkdir tmp
        java -Djava.io.tmpdir=./tmp -jar $PICARD_JAR \
            CollectAlignmentSummaryMetrics \
            R=${ref_fasta} \
            I=${input_bam} \
            O=${sample_name}.alignment_metrics.txt;
    }

    output {
        File alignment_metrics = "${sample_name}.alignment_metrics.txt"
    }

    runtime {
        docker: "docker.io/hsphqbrc/gatk-mutect2-workflow-tools:1.0"
        cpu: 2
        memory: "4 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task deduplicate_bam {
    File input_bam
    File input_bam_index
    String sample_name

    # runtime commands
    Int disk_size = 150

    command {
        mkdir tmp
        java -Djava.io.tmpdir=./tmp -Xmx3000m -jar $PICARD_JAR \
            MarkDuplicates \
            INPUT=${input_bam} \
            OUTPUT=${sample_name}.bam \
            ASSUME_SORTED=TRUE \
            TMP_DIR=/tmp \
            REMOVE_DUPLICATES=TRUE \
            METRICS_FILE=${sample_name}.dedup_metrics.txt \
            VALIDATION_STRINGENCY=LENIENT;
        samtools index ${sample_name}.bam;
    }

    output {
        File sorted_bam = "${sample_name}.bam"
        File sorted_bam_index = "${sample_name}.bam.bai"
        File deduplication_metrics = "${sample_name}.dedup_metrics.txt"
    }

    runtime {
        docker: "docker.io/hsphqbrc/gatk-mutect2-workflow-tools:1.0"
        cpu: 2
        memory: "4 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task base_recalibrator {
    File input_bam
    File input_bam_index
    String sample_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File dbsnp
    File dbsnp_index
    File known_indels
    File known_indels_index

    # runtime commands
    Int disk_size = 150

    command {
        mkdir tmp
        java -Djava.io.tmpdir=./tmp -Xmx4000m -jar $GATK_JAR \
            BaseRecalibrator \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -known-sites ${dbsnp} \
            -known-sites ${known_indels} \
            -O recal_data.table;
    }
    
    output {
        File recalibration_report = "recal_data.table"
    }

    runtime {
        docker: "docker.io/hsphqbrc/gatk-mutect2-workflow-tools:1.0"
        cpu: 4
        memory: "6 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task apply_recalibration {
    File input_bam
    File input_bam_index
    File recalibration_report
    String sample_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    # runtime commands
    Int disk_size = 150

    command {
        mkdir tmp
        java -Djava.io.tmpdir=./tmp -Xmx4000m -jar $GATK_JAR \
            ApplyBQSR \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${sample_name}.bam \
            -bqsr-recal-file ${recalibration_report};
        samtools index ${sample_name}.bam;
    }

    output {
        File recalibrated_bam = "${sample_name}.bam"
        File recalibrated_bam_index = "${sample_name}.bam.bai"
    }

    runtime {
        docker: "docker.io/hsphqbrc/gatk-mutect2-workflow-tools:1.0"
        cpu: 4
        memory: "6 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task conpair_pileup {
    File input_bam
    File input_bam_index
    File input_dedup_bam
    File input_dedup_bam_index
    String use_dedup
    String sample_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    Int disk_size = 250

    command {
        mkdir tmp
        TMP_DIR=`pwd`/tmp
        if [ "${use_dedup}" = "true" ]
        then
            /opt/software/Conpair-0.2/scripts/run_gatk_pileup_for_sample.py \
                -B ${input_dedup_bam} \
                -O ${sample_name}.pileup \
                --reference ${ref_fasta} \
                --gatk /opt/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
        else
            /opt/software/Conpair-0.2/scripts/run_gatk_pileup_for_sample.py \
                -B ${input_bam} \
                -O ${sample_name}.pileup \
                --reference ${ref_fasta} \
                --gatk /opt/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
        fi
    }

    runtime {
        docker: "docker.io/hsphqbrc/gatk-mutect2-workflow-tools:1.0"
        cpu: 4
        memory: "16 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }

    output {
        File output_pileup = "${sample_name}.pileup"
    }
}

task conpair_concordance {
    File tumor_pileup
    File normal_pileup
    String sample_name
    
    # runtime commands
    Int disk_size = 250

    command {
        mkdir tmp
        TMP_DIR=`pwd`/tmp
        /opt/software/Conpair-0.2/scripts/verify_concordance.py \
            -T ${tumor_pileup} \
            -N ${normal_pileup} \
            --outfile ${sample_name}.concordance_metrics.txt > stdout_log.txt;
        if [ ! -f ${sample_name}.concordance_metrics ]; then
            mv stdout_log.txt ${sample_name}.concordance_metrics.txt;
        fi
    }

    runtime {
        docker: "docker.io/hsphqbrc/gatk-mutect2-workflow-tools:1.0"
        cpu: 2
        memory: "8 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }

    output {
        File concordance_metrics = "${sample_name}.concordance_metrics.txt"
    }
}

task conpair_contamination {
    File tumor_pileup
    File normal_pileup
    String sample_name
    
    # runtime commands
    Int disk_size = 250

    command {
        mkdir tmp
        TMP_DIR=`pwd`/tmp
        /opt/software/Conpair-0.2/scripts/estimate_tumor_normal_contamination.py \
            -T ${tumor_pileup} \
            -N ${normal_pileup} \
            --outfile ${sample_name}.contamination_metrics.txt 2> catch_stderr.txt;
        if grep -n "ValueError" catch_stderr.txt; then
            echo "Error caught in conpair contamination. Likely cause is lack of shared depth between normal and tumor samples." > ${sample_name}.contamination_metrics.txt
        fi
    }

    runtime {
        docker: "docker.io/hsphqbrc/gatk-mutect2-workflow-tools:1.0"
        cpu: 2
        memory: "8 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }

    output {
        File contamination_metrics = "${sample_name}.contamination_metrics.txt"
    }
}

task coverage_metrics {
    File input_bam
    File input_bam_index
    File input_dedup_bam
    File input_dedup_bam_index
    String use_dedup
    String sample_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_exon_intervals

    Int disk_size = 250

    command {
        mkdir tmp
        if [ "${use_dedup}" = "true" ]
        then
            java -Djava.io.tmpdir=./tmp -jar -Xmx6000m -jar $PICARD_JAR \
                CollectWgsMetrics \
                R=${ref_fasta} \
                I=${input_dedup_bam} \
                INTERVALS=${ref_exon_intervals} \
                O=${sample_name}.coverage_metrics.txt \
                COVERAGE_CAP=500;
        else
            java -Djava.io.tmpdir=./tmp -jar -Xmx6000m -jar $PICARD_JAR \
                CollectWgsMetrics \
                R=${ref_fasta} \
                I=${input_bam} \
                INTERVALS=${ref_exon_intervals} \
                O=${sample_name}.coverage_metrics.txt \
                COVERAGE_CAP=500;
        fi
    }

    runtime {
        docker: "docker.io/hsphqbrc/gatk-mutect2-workflow-tools:1.0"
        cpu: 2
        memory: "8 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }

    output {
        File coverage_metrics = "${sample_name}.coverage_metrics.txt"
    }
}

task haplotypecaller {
    File input_bam
    File input_bam_index
    File input_dedup_bam
    File input_dedup_bam_index
    String use_dedup
    String sample_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    String interval

    # runtime commands
    Int disk_size = 250

    command {
        mkdir tmp
        if [ "${use_dedup}" = "true" ]
        then
            java -Djava.io.tmpdir=./tmp -Xmx8000m -jar $GATK_JAR \
                HaplotypeCaller \
                -R ${ref_fasta} \
                -I ${input_dedup_bam} \
                -O ${sample_name}.vcf \
                -L ${interval};
        else
            java -Djava.io.tmpdir=./tmp -Xmx8000m -jar $GATK_JAR \
                HaplotypeCaller \
                -R ${ref_fasta} \
                -I ${input_bam} \
                -O ${sample_name}.vcf \
                -L ${interval};
        fi
    }

    output {
        File output_vcf = "${sample_name}.vcf"
    }

    runtime {
        docker: "docker.io/hsphqbrc/gatk-mutect2-workflow-tools:1.0"
        cpu: 8
        memory: "12 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task mutect {
    File normal_bam
    File normal_bam_index
    File tumor_bam
    File tumor_bam_index
    File normal_dedup_bam
    File normal_dedup_bam_index
    File tumor_dedup_bam
    File tumor_dedup_bam_index
    String use_dedup
    String normal_sample_name
    String tumor_sample_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File gnomad
    File gnomad_index
    String interval

    # runtime commands
    Int disk_size = 500

    command {
        mkdir tmp
        if [ "${use_dedup}" = "true" ]
        then
            java -Djava.io.tmpdir=./tmp -Xmx8000m -jar $GATK_JAR \
                Mutect2 \
                -R ${ref_fasta} \
                -I ${tumor_dedup_bam} \
                -tumor ${tumor_sample_name} \
                -I ${normal_dedup_bam} \
                -normal ${normal_sample_name} \
                --germline-resource ${gnomad} \
                -L ${interval} \
                -O ${tumor_sample_name}.vcf;
        else
            java -Djava.io.tmpdir=./tmp -Xmx8000m -jar $GATK_JAR \
                Mutect2 \
                -R ${ref_fasta} \
                -I ${tumor_bam} \
                -tumor ${tumor_sample_name} \
                -I ${normal_bam} \
                -normal ${normal_sample_name} \
                --germline-resource ${gnomad} \
                -L ${interval} \
                -O ${tumor_sample_name}.vcf;
        fi
    }

    output {
        File output_vcf = "${tumor_sample_name}.vcf"
    }

    runtime {
        docker: "docker.io/hsphqbrc/gatk-mutect2-workflow-tools:1.0"
        cpu: 4
        memory: "12 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task merge_vcf {
    Array[File] input_vcfs
    String sample_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    # runtime commands
    Int disk_size = 500

    command {
        mkdir tmp
        java -Djava.io.tmpdir=./tmp -Xmx3000m -jar $PICARD_JAR \
            MergeVcfs \
            INPUT=${sep=' INPUT=' input_vcfs} \
            OUTPUT=${sample_name}.vcf;
        grep -v "^#" ${sample_name}.vcf > is_vcf_empty.txt
    }

    output {
        File output_vcf = "${sample_name}.vcf"
        File is_vcf_empty = "is_vcf_empty.txt"
    }

    runtime {
        docker: "docker.io/hsphqbrc/gatk-mutect2-workflow-tools:1.0"
        cpu: 2
        memory: "4 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
