workflow test_qc {
    Array[File] tumor_alignment_metrics
    Array[File] tumor_dedup_metrics
    Array[File] tumor_r1_fastqc_zips
    Array[File] tumor_r2_fastqc_zips
    Array[File] normal_alignment_metrics
    Array[File] normal_dedup_metrics
    Array[File] normal_r1_fastqc_zips
    Array[File] normal_r2_fastqc_zips

    call create_multi_qc{
        input:
            tumor_alignment_metrics = tumor_alignment_metrics,
            tumor_dedup_metrics = tumor_dedup_metrics,
            normal_alignment_metrics = normal_alignment_metrics,
            normal_dedup_metrics = normal_dedup_metrics,
            tumor_r1_fastqc_zips = tumor_r1_fastqc_zip,
            tumor_r2_fastqc_zips = tumor_r2_fastqc_zip,
            normal_r1_fastqc_zips = normal_r1_fastqc_zip,
            normal_r2_fastqc_zips = normal_r2_fastqc_zip
    }

    output {
        File multiqc_zip = create_multi_qc.report
    }
}

task create_multi_qc {
    Array[File] tumor_alignment_metrics
    Array[File] tumor_dedup_metrics
    Array[File] tumor_r1_fastqc_zips
    Array[File] tumor_r2_fastqc_zips
    Array[File] normal_alignment_metrics
    Array[File] normal_dedup_metrics
    Array[File] normal_r1_fastqc_zips
    Array[File] normal_r2_fastqc_zips
    Array[File] concordance_metrics
    Array[File] contamination_metrics

    # Runtime parameters
    Int disk_size = 100

    command {
        multiqc .
    }

    output {
        File report = "multiqc_report.html"
    }

    runtime {
        zones: "us-east4-c"
        docker: "docker.io/hsphqbrc/gatk-mutect2-workflow-tools:1.0"
        cpu: 2
        memory: "4 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task generate_report {
    Array[String] r1_files
    Array[String] r2_files
    String genome
    String dedup_bool
    String git_repo_url
    String git_commit_hash

    # Runtime parameters
    Int disk_size = 50

    command <<<
        # make a json file with various parameters:
        echo "{" >> config.json
        echo '"genome": "${genome}",' >>config.json
        echo '"dedup_bool": "${dedup_bool}",' >> config.json
        echo '"git_repo": "${git_repo_url}",' >>config.json
        echo '"git_commit": "${git_commit_hash}"' >>config.json
        echo "}" >> config.json

        generate_report.py \
          -r1 ${sep=" " r1_files} \
          -r2 ${sep=" " r2_files} \
          -j config.json \
          -t /opt/report/report.md \
          -o completed_report.md

        pandoc \
            -H /opt/report/report.css \
            -s completed_report.md \
            -o analysis_report.html
    >>>

    output {
        File report = "analysis_report.html"
    }

    runtime {
        zones: "us-east4-c"
        docker: "docker.io/hsphqbrc/gatk-mutect2-workflow-tools:1.0"
        cpu: 2
        memory: "6 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}