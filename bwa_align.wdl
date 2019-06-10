workflow align_test {
    File r1_fastq
    File? r2_fastq
    String sample_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_bwt
    File ref_sa
    File ref_amb
    File ref_ann
    File ref_pac

    call perform_align {
        input:
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq,
            sample_name = sample_name,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            ref_bwt = ref_bwt,
            ref_sa = ref_sa,
            ref_amb = ref_amb,
            ref_ann = ref_ann,
            ref_pac = ref_pac
    }

    output {
        File bam = perform_align.sorted_bam
        File bam_index = perform_align.sorted_bam_index
    }
}

task perform_align {
    # align to the reference genome using the bwa aligner
    # The bwa alignment process produces a position-sorted BAM file
    # and index file

    # Input params passed by a parent Workflow:
    # We require that there is a single R1 fastq and possibly a single
    # R2 fastq for paired sequencing.
    # Genome is a string that matches one of the keys in the Map below
    # sample_name helps with naming files for eventual concatenation.

    File r1_fastq
    File? r2_fastq
    String sample_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_bwt
    File ref_sa
    File ref_amb
    File ref_ann
    File ref_pac

    # runtime commands
    Int disk_size = 500

    command {
        RGID="N000.1"
        RGPU="X0000XXXX000000.1.NNNN"
        RGPL="illumina"
        RGLB="XXX"
        RGSM="${sample_name}"
        RGCN="unknown"
        bwa mem -t 8 ${ref_fasta} ${r1_fastq} ${r2_fastq} \
        | samtools view -bht ${ref_fasta} - \
        | samtools sort -o ${sample_name}.preid.bam -;
        samtools index ${sample_name}.preid.bam;
        java -Xmx2500m -jar $PICARD_JAR \
            AddOrReplaceReadGroups \
            I=${sample_name}.preid.bam \
            O=${sample_name}.bam \
            RGID=$RGID \
            RGPU=$RGPU \
            RGPL=$RGPL \
            RGLB=$RGLB \
            RGSM=$RGSM \
            RGCN=$RGCN;
        samtools index ${sample_name}.bam;
    }

    output {
        File sorted_bam = "${sample_name}.bam"
        File sorted_bam_index = "${sample_name}.bam.bai"
    }

    runtime {
        docker: "docker.io/hsphqbrc/gatk-mutect2-workflow-tools:1.0"
        cpu: 8
        memory: "12 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}