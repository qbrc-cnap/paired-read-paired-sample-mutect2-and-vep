workflow vep_test {
    File input_vcf
    String sample_name
    String species
    File vep_cache_tar

    call vep_annotate {
        input:
            input_vcf = input_vcf,
            sample_name = sample_name,
            species = species,
            vep_cache_tar = vep_cache_tar
    }

    output {
        File output_vcf = vep_annotate.output_vcf
        File vcf_stats = vep_annotate.vcf_stats
    }
}

task vep_annotate {
    File input_vcf
    String sample_name
    String species
    File vep_cache_tar
    File is_vcf_empty

    # runtime commands
    Int disk_size = 200

    command {
        is_vcf_empty=$(cat ${is_vcf_empty})
        if [ ${is_vcf_empty} -eq 0 ]
        then
            echo "The Mutect2 VCF output was empty, so VEP would fail if applied to it. No annotated VCF can be provided then." > ${sample_name}.vep.tsv;
            echo ""The Mutect2 VCF output was empty, so VEP would fail if applied to it. No annotated VCF can be provided then." > ${sample_name}.vep_stats.html;
        else
            tar xf ${vep_cache_tar};
            #working_dir=$(pwd)
            /opt/vep/src/ensembl-vep/vep \
                -i ${input_vcf} \
                -o ${sample_name}.vep.tsv \
                --stats_file ${sample_name}.vep_stats.html \
                --species ${species} \
                --dir ./vep_data \
                --tab \
                --cache \
                --offline \
                --everything;
        fi
    }

    output {
        File output_vcf = "${sample_name}.vep.tsv"
        File vcf_stats = "${sample_name}.vep_stats.html"
    }

    runtime {
        zones: "us-east4-c"
        docker: "docker.io/ensemblorg/ensembl-vep:release_95.0"
        cpu: 4
        memory: "12 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
