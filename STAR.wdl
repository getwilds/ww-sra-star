version 1.0
## Note:  to keep this workflow simple, we have designed this as written to perform STAR analysis
## for hg38 datasets using a specific module avilable at Fred Hutch.

#### WORKFLOW DEFINITION

workflow SRA_STAR2Pass {
  input { 
    Array[String] sra_id_list = ["SRR10724344"]
    ## This sample takes about 11 hours to complete depending on the download speed from SRA
  }

  String referenceGenome = "hg38"

  scatter ( id in sra_id_list ){
    call fastqdump {
      input:
        sra_id = id,
        ncpu = 12
    }

    call STARalignTwoPass {
      input:
        r1fastq = fastqdump.R1end,
        r2fastq = fastqdump.R2end,
        base_file_name = id,
        referenceGenome = referenceGenome,
        cpu = 12
    }
  } # End scatter 

  # Outputs that will be retained when execution is complete
  output {
    Array[File] output_bam = STARalignTwoPass.bam
    Array[File] output_bai = STARalignTwoPass.bai
    Array[File] output_geneCounts = STARalignTwoPass.geneCounts
    Array[File] output_log_final = STARalignTwoPass.log_final
    Array[File] output_log_progress = STARalignTwoPass.log_progress
    Array[File] output_log = STARalignTwoPass.log
    Array[File] output_SJ = STARalignTwoPass.SJout
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
      # pfastq-dump \
      #   -t ~{ncpu} \
      #   --gzip \
      #   --split-files \
      #   -s "~{sra_id}" -O ./
    else
      touch paired_file
      parallel-fastq-dump \
        --sra-id ~{sra_id} \
        --threads ~{ncpu} \
        --outdir ./ \
        --gzip
      # pfastq-dump \
      #   -t ~{ncpu} \
      #   --gzip \
      #   -s "~{sra_id}" -O ./
    fi
  >>>

  output {
    File R1end = "~{sra_id}_1.fastq.gz"
    File R2end = "~{sra_id}_2.fastq.gz"
    String paired_end = read_string('paired_file')
  }

  runtime {
    memory: 2 * ncpu + " GB"
    # docker: "ncbi/sra-tools:3.0.0"
    # docker: 'ghcr.io/stjude/abralab/sratoolkit:v3.0.0'
    modules: "parallel-fastq-dump/0.6.7-GCCcore-11.2.0"
    cpu: ncpu
  }
}

task STARalignTwoPass {
  input {
    File r1fastq
    File r2fastq
    String base_file_name
    String referenceGenome
    Int cpu
  }

  command <<<
    set -eo pipefail
    STAR \
      --genomeDir /shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/STAR2Index \
      --readFilesIn ~{r1fastq} ~{r2fastq} \
      --runThreadN ~{cpu} \
      --readFilesCommand zcat \
      --sjdbOverhang 100 \
      --outSAMtype BAM SortedByCoordinate \
      --twopassMode Basic \
      --quantMode GeneCounts \
      --quantTranscriptomeBAMcompression 5 
    mv Aligned.sortedByCoord.out.bam "~{base_file_name}.~{referenceGenome}.Aligned.sortedByCoord.out.bam"
    mv ReadsPerGene.out.tab "~{base_file_name}.~{referenceGenome}.ReadsPerGene.out.tab"
    mv Log.final.out "~{base_file_name}.~{referenceGenome}.Log.final.out"
    samtools index "~{base_file_name}.~{referenceGenome}.Aligned.sortedByCoord.out.bam"
  >>>

  output {
    File bam = "~{base_file_name}.~{referenceGenome}.Aligned.sortedByCoord.out.bam"
    File bai = "~{base_file_name}.~{referenceGenome}.Aligned.sortedByCoord.out.bam.bai"
    File geneCounts = "~{base_file_name}.~{referenceGenome}.ReadsPerGene.out.tab"
    File log_final = "~{base_file_name}.~{referenceGenome}.Log.final.out"
    File log_progress = "Log.progress.out"
    File log = "Log.out"
    File SJout = "SJ.out.tab"
  }

  runtime {
    modules: "STAR/2.7.6a-GCC-10.2.0 SAMtools/1.11-GCC-10.2.0"
    memory: 2 * cpu + "GB"
    cpu: cpu
  }
}
