workflow PairedMatchedMutect2AndVepWorkflow {
    # This workflow is a 'super' workflow that parallelizes
    # HaplotypeCaller and VEP analysis over multiple samples.

    # Input files
    Array[File] r1_files
    Array[File] r2_files
    File match_annotations
    Boolean use_dedup
    
    # Reference files
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_bwt
    File ref_sa
    File ref_amb
    File ref_ann
    File ref_pac

    # Inputs for GATK BQSR
    File dbsnp
    File dbsnp_index
    File known_indels
    File known_indels_index
    File gnomad
    File gnomad_index

    # Inputs for HaplotypeCaller scatter
    File contig_list
    Array[String] contigs = read_lines(contig_list)

    # Inputs for VEP
    String vep_species
    File vep_cache_tar
    
    # Other
    String output_zip_name
    String genome
    String git_repo_url
    String git_commit_hash


    Array[Pair[File, File]] fastq_pairs = zip(r1_files, r2_files)
    scatter(item in fastq_pairs){

        call assert_valid_fastq {
            input:
                r1_file = item.left,
                r2_file = item.right
        }
    }

    call assert_valid_annotation {
        input:
            r1_files = r1_files,
            r2_files = r2_files,
            match_annotations = match_annotations
    }
}

task assert_valid_fastq {

    File r1_file
    File r2_file
    Int disk_size = 100

    command <<<
        python3 /opt/software/precheck/check_fastq.py -r1 ${r1_file} -r2 ${r2_file}
    >>>

    runtime {
        docker: "docker.io/hsphqbrc/gatk-mutect2-workflow-tools:1.0"
        cpu: 2
        memory: "6 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task assert_valid_annotation {
    Array[String] r1_files
    Array[String] r2_files
    File match_annotations
    
    Int disk_size = 20

    command <<<
        python3 /opt/software/precheck/check_annotations.py \
            -r1 ${sep=" " r1_files} \
            -r2 ${sep=" " r2_files} \
            -annot ${match_annotations}
    >>>

    runtime {
        docker: "docker.io/hsphqbrc/gatk-mutect2-workflow-tools:1.0"
        cpu: 1
        memory: "2 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}