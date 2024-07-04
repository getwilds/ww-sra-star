version 1.0
# Pulls down paired fastq's for SRA ID's provided 
# and aligns them using the STAR two-pass methodology

#### STRUCT DEFINITIONS

struct RefGenome {
    File star_tar
    String ref_name
}

#### WORKFLOW DEFINITION

workflow SRA_STAR2Pass {
  input { 
    Array[String] sra_id_list
    RefGenome ref_genome
  }

  scatter ( id in sra_id_list ){
    call fastqdump {
      input:
        sra_id = id,
        ncpu = 12
    }

    call STARalignTwoPass {
      input:
        ref_genome = ref_genome,
        r1fastq = fastqdump.r1_end,
        r2fastq = fastqdump.r2_end,
        base_file_name = id,
        cpu = 12
    }
  } # End scatter 

  # Outputs that will be retained when execution is complete
  output {
    Array[File] output_bam = STARalignTwoPass.bam
    Array[File] output_bai = STARalignTwoPass.bai
    Array[File] output_gene_counts = STARalignTwoPass.gene_counts
    Array[File] output_log_final = STARalignTwoPass.log_final
    Array[File] output_log_progress = STARalignTwoPass.log_progress
    Array[File] output_log = STARalignTwoPass.log
    Array[File] output_SJ = STARalignTwoPass.SJout
  }

  parameter_meta {
    sra_id_list: "list of SRA sample ID's to be pulled down and aligned"
    ref_genome: "struct describing the reference genome to be used during alignment"

    output_bam: "array of aligned bam files for each sample"
    output_bai: "array of corresponding index files for each aligned bam file"
    output_gene_counts: "array of text files containing the number of reads in each gene for each sample"
    output_log_final: "array of text files containing an overarching summary of the analysis performed for each sample"
    output_log_progress: "array of text files containing a detailed progress report for each sample"
    output_log: "array of text files containing STAR's raw command line output for each sample"
    output_SJ: "array of text files containing splice junction details for each sample being analyzed"
  }
} # End Workflow

#### TASK DEFINITIONS

task fastqdump {
  input {
    String sra_id
    Int ncpu = 12
  }

  command <<<
    set -eo pipefail
    # check if paired ended
    numLines=$(fastq-dump -X 1 -Z --split-spot "~{sra_id}" | wc -l)
    paired_end="false"
    if [ $numLines -eq 8 ]; then
      paired_end="true"
    fi
    # perform fastqdump
    if [ $paired_end == 'true' ]; then
      echo true > paired_file
      parallel-fastq-dump \
        --sra-id ~{sra_id} \
        --threads ~{ncpu} \
        --outdir ./ \
        --split-files \
        --gzip
    else
      touch paired_file
      parallel-fastq-dump \
        --sra-id ~{sra_id} \
        --threads ~{ncpu} \
        --outdir ./ \
        --gzip
    fi
  >>>

  output {
    File r1_end = "~{sra_id}_1.fastq.gz"
    File r2_end = "~{sra_id}_2.fastq.gz"
    String paired_end = read_string('paired_file')
  }

  runtime {
    memory: 2 * ncpu + " GB"
    docker: "getwilds/pfastqdump:0.6.7"
    cpu: ncpu
  }

  parameter_meta {
    sra_id: "SRA ID of the sample to be downloaded via parallel-fastq-dump"
    ncpu: "number of cpus to use during download"

    r1_end: "R1 fastq file downloaded for the sample in question"
    r2_end: "R2 fastq file downloaded for the sample in question"
    paired_end: "string indicating whether the sample in question used paired-end read sequencing"
  }
}

task STARalignTwoPass {
  input {
    RefGenome ref_genome
    File r1fastq
    File r2fastq
    String base_file_name
    Int cpu
  }

  String star_db_dir = basename(ref_genome.star_tar, ".tar.gz")

  command <<<
    set -eo pipefail
    tar -xzf ~{ref_genome.star_tar}
    STAR \
      --genomeDir ~{star_db_dir} \
      --readFilesIn ~{r1fastq} ~{r2fastq} \
      --runThreadN ~{cpu} \
      --readFilesCommand zcat \
      --sjdbOverhang 100 \
      --outSAMtype BAM SortedByCoordinate \
      --twopassMode Basic \
      --quantMode GeneCounts \
      --quantTranscriptomeBAMcompression 5 
    mv Aligned.sortedByCoord.out.bam "~{base_file_name}.~{ref_genome.ref_name}.Aligned.sortedByCoord.out.bam"
    mv ReadsPerGene.out.tab "~{base_file_name}.~{ref_genome.ref_name}.ReadsPerGene.out.tab"
    mv Log.final.out "~{base_file_name}.~{ref_genome.ref_name}.Log.final.out"
    samtools index "~{base_file_name}.~{ref_genome.ref_name}.Aligned.sortedByCoord.out.bam"
  >>>

  output {
    File bam = "~{base_file_name}.~{ref_genome.ref_name}.Aligned.sortedByCoord.out.bam"
    File bai = "~{base_file_name}.~{ref_genome.ref_name}.Aligned.sortedByCoord.out.bam.bai"
    File gene_counts = "~{base_file_name}.~{ref_genome.ref_name}.ReadsPerGene.out.tab"
    File log_final = "~{base_file_name}.~{ref_genome.ref_name}.Log.final.out"
    File log_progress = "Log.progress.out"
    File log = "Log.out"
    File SJout = "SJ.out.tab"
  }

  runtime {
    docker: "getwilds/star:2.7.6a"
    memory: 2 * cpu + "GB"
    cpu: cpu
  }

  parameter_meta {
    ref_genome: ""
    r1fastq: ""
    r2fastq: ""
    base_file_name: ""
    cpu: ""

    bam: "aligned bam file for the sample in question"
    bai: "corresponding index file for the aligned bam"
    gene_counts: "text file containing the number of reads in each gene"
    log_final: "text file containing an overarching summary of the analysis performed"
    log_progress: "text file containing a detailed progress report for the analysis"
    log: "text file containing STAR's raw command line output"
    SJout: "text file containing splice junction details"
  }
}
