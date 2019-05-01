## GRCh38.95
* primary assembly
    * ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
* FASTA index
    * samtools: Version: 1.9 (using htslib 1.9)
    ```bash
    samtools index Homo_sapiens.GRCh38.dna.primary_assembly.fa
    ```
* FASTA sequence dictionary
    * java version: java version "1.8.0_192"
        * Java(TM) SE Runtime Environment (build 1.8.0_192-b12)
        * Java HotSpot(TM) 64-Bit Server VM (build 25.192-b12, mixed mode)
    * picard version: 2.18.17
    ```bash
    java -jar picard.jar \
        CreateSequenceDictionary \
        R=Homo_sapiens.GRCh38.dna.primary_assembly.fa
    ```
* bwa index
    * bwa: Version: 0.7.17-r1188
    ```bash
    bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa
    ```
* broad resources
    * ftp://ftp.broadinstitute.org/bundle/hg38/
        * Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
        * Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
        * dbsnp_146.hg38.vcf.gz
        * dbsnp_146.hg38.vcf.gz.tbi
    * ftp://ftp.broadinstitute.org/bundle/Mutect2/
        * af-only-gnomad.hg38.vcf
## GRCh37.95
* primary assembly
    * ftp://ftp.ensembl.org/pub/grch37/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
* broad resources
    * ftp://ftp.broadinstitute.org/bundle/hg19/
        * Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
        * Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx
        * dbsnp_138.hg19.vcf
        * dbsnp_138.hg19.vcf.idx
    * http://bioinfo5pilm46.mit.edu/software/GATK/resources/
        * af-only-gnomad.hg19.vcf