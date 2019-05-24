Report for alignment and differential expression analysis
---

This document discusses the steps that were performed in the analysis pipeline.  It also describes the format of the output files and some brief interpretation.  For more detailed questions about interpretation of results, consult the documentation of the various tools.


## Outputs:

This section describes the contents of the delivered results.

#### Variant calls

Individual variant files (in VCF format, ending with ".vcf") are available for download, but are provided separately. You can find an breakdown of the VCF format at <https://samtools.github.io/hts-specs/VCFv4.2.pdf>.

#### Main results

The main results are contained in a zip-archive and should be downloaded an "unzipped" on your local computer.  It contains several sub-directories which contain files produced in each step of the pipeline.

- **QC**
    - This directory contains an interactive HTML-based QC report which summarizes read quality, alignment quality, and other metrics.  It was produced by MultiQC, and information can be found at <https://multiqc.info/>.
        - QC of the FASTQ themselves with FASTQC: <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>
        - QC of the read alignment with Picard: <https://broadinstitute.github.io/picard/>
        - QC of the read duplication content of the libraries with Picard: <https://broadinstitute.github.io/picard/>
        - QC of the per-base depth coverage of the exons
            - exons are defined Ensembl's gff3 exon regions
        - QC of the base quality score recalibration as determined by GATK: <https://software.broadinstitute.org/gatk/>
- **Annotated TSV**
    - This directory contains all the individual annotated tab delimited variant files produced by the annotation of the VCF with Ensembl's VEP.
- **Logs**
    - This contains logs and summaries produced by the various tools.  These can be used for troubleshooting, if necessary.


## Methods:

Input fastq-format files are aligned to the {{genome}} reference genome using the bwa aligner ({{bwa_version}}) [1].  BAM-format alignment files were filtered to retain only the primary-aligned reads using samtools ({{samtools_version}}) [2].  Additionally, "de-duplicated" versions of the primary-filtered BAM files were created using PicardTools' MarkDuplicates software ({{picard_mark_duplicates_version}}) [3].

Quality-control software included FastQC ({{fastqc_version}}), Picard ({{picard_alignment_metrics_version}}), and MultiQC ({{multiqc_version}}).  Please see the respective references for interpretation of output information and figures.

Depending on the library type used, it makes sense to use the "unfiltered" BAM files. **For low diversity libraries (such as amplicons) or enzymatically produced library fragmentation (such as Agilent's Haloplex), it is REQUIRED to use the "unfiltered" BAM for variant calling.** The more typical high diversity libraries such as Agilent's SureSelect are simply recommended to use the "deduplicated" libraries. **The chosen deduplication option will be applied to all samples.** By default, we only perform variant calling on the library type picked in the pre-analysis screen, so it is important to have known prior to analysis which method to use. However, it is helpful to see the the amount of duplication present in high diversity libraries to ascertain library quality, so we provide the deduplication QC plots regardless of the choice picked in the pre-analysis screen.

Further manipulation of the read data is the base quality score recalibration of the libraries. The read quality scores are subject to systemic biases is the quality score in the data. Base quality score recalibration as applied by GATK [4-7] uses machine learning to model these errors empirically and adjust the quality scores accordingly. This results in improved variant calling. Variant calling is done with GATK's Mutect2. You can find further information on how Mutect2 identified variants at: <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.4/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php>. Additionally, we use the gnomad allele frequency resources as additionaly data for Mutect2 to imporove variant calling accuracy. (Please note that we do not use a panel of normals as an additional resource for variant calling.) **IMPORTANT: we do not perform any filtering on the variant calls in the VCF.** The variant caller by default is permissive to maximize sensitivity. VCF filtering is a case by case matter, and **we recommend looking over the resultant variants with a critical eye and applying at least a qualitative approach to ascertaining the validity of the variant call.**

Variant annotation is done with Ensembl's VEP software package [8]. Among the annotation included are:
- SIFT amino acid substitution predictions: <https://sift.bii.a-star.edu.sg/> [9]
- PolyPhen prediction of functional effects of amino acid substitution: <http://genetics.bwh.harvard.edu/pph2/> [10]
- gnomAD population frequencies <https://gnomad.broadinstitute.org/> [11]

## Inputs:

Samples and sequencing fastq-format files:

{% for obj in file_display %}
  - {{obj.sample_name}}
    - R1 fastq: {{obj.r1}}
    - R2 fastq: {{obj.r2}}
{% endfor %}

Deduplication used for input into variant caller: {{dedup_bool}}

## Version control:
To facilitate reproducible analyses, the analysis pipeline used to process the data is kept under git-based version control.  The repository for this workflow is at 

<{{git_repo}}>

and the commit version was {{git_commit}}.

This allows us to run the *exact* same pipeline at any later time, discarding any updates or changes in the process that may have been added. 


#### References:

[1] Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168] 

[2] Li H. and Handsaker B. and Wysoker A. and Fennell T. and Ruan J. and Homer N. and Marth G. and Abecasis G. and Durbin R. and 1000 Genome Project Data Processing Subgroup.  The Sequence alignment/map (SAM) format and SAMtools.  Bioinformatics. 2009.

[3] <http://broadinstitute.github.io/picard/>

[4] <https://software.broadinstitute.org/gatk/documentation/article?id=11081>

[5] The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA, 2010 GENOME RESEARCH 20:1297-303 

[6] A framework for variation discovery and genotyping using next-generation DNA sequencing data DePristo M, Banks E, Poplin R, Garimella K, Maguire J, Hartl C, Philippakis A, del Angel G, Rivas MA, Hanna M, McKenna A, Fennell T, Kernytsky A, Sivachenko A, Cibulskis K, Gabriel S, Altshuler D, Daly M, 2011 NATURE GENETICS 43:491-498 

[7] From FastQ Data to High-Confidence Variant Calls: The Genome Analysis Toolkit Best Practices Pipeline Van der Auwera GA, Carneiro M, Hartl C, Poplin R, del Angel G, Levy-Moonshine A, Jordan T, Shakir K, Roazen D, Thibault J, Banks E, Garimella K, Altshuler D, Gabriel S, DePristo M, 2013 CURRENT PROTOCOLS IN BIOINFORMATICS 43:11.10.1-11.10.33 

[8] McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F. The Ensembl Variant Effect Predictor. Genome Biology Jun 6;17(1):122. (2016) doi:10.1186/s13059-016-0974-4 

[9] SIFT missense predictions for genomes. Nat Protocols 2016; 11:1-9

[10] Adzhubei IA, Schmidt S, Peshkin L, Ramensky VE, Gerasimova A, Bork P, Kondrashov AS, Sunyaev SR. Nat Methods 7(4):248-249 (2010).

[11] Variation across 141,456 human exomes and genomes reveals the spectrum of loss-of-function intolerance across human protein-coding genes.  View ORCID ProfileKonrad J Karczewski, Laurent C Francioli, Grace Tiao, Beryl B Cummings, Jessica Alföldi, Qingbo Wang, Ryan L Collins, Kristen M Laricchia, Andrea Ganna, Daniel P Birnbaum, Laura D Gauthier, Harrison Brand, Matthew Solomonson, Nicholas A Watts, Daniel Rhodes, Moriel Singer-Berk, Eleanor G Seaby, Jack A Kosmicki, Raymond K Walters, Katherine Tashman, Yossi Farjoun, Eric Banks, Timothy Poterba, Arcturus Wang, Cotton Seed, Nicola Whiffin, Jessica X Chong, Kaitlin E Samocha, Emma Pierce-Hoffman, Zachary Zappala, Anne H O'Donnell-Luria, Eric Vallabh Minikel, Ben Weisburd, Monkol Lek, James S Ware, Christopher Vittal, Irina M Armean, Louis Bergelson, Kristian Cibulskis, Kristen M Connolly, Miguel Covarrubias, Stacey Donnelly, Steven Ferriera, Stacey Gabriel, Jeff Gentry, Namrata Gupta, Thibault Jeandet, Diane Kaplan, Christopher Llanwarne, Ruchi Munshi, Sam Novod, Nikelle Petrillo, David Roazen, Valentin Ruano-Rubio, Andrea Saltzman, Molly Schleicher, Jose Soto, Kathleen Tibbetts, Charlotte Tolonen, Gordon Wade, Michael E Talkowski, The Genome Aggregation Database Consortium, Benjamin M Neale, Mark J Daly, Daniel G MacArthur. 2019. doi: https://doi.org/10.1101/531210 