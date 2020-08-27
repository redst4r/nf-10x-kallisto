//
//  Turning a bam file from a 10x sequencing dataset back into
//  the fastq records!
//  Usage:
//    nextflow run mkfastq.nf
//      --bamfile <.bam>
//      --outdir <directory>
//      --publish_mode <copy|symlink>   <- copy creates copies of the files in outdir, symlink just links from the tmp files (saving diskspace).
//  If using --publish_mode symlink,  make sure to copy the files outside of nextflow before doing nextflow cleanup

BAM2FASTQ_BINARY='/home/mstrasse/bamtofastq-1.2.0'
params.cpus = 4
if (!params.outdir){
  exit 1, "--outdir not set!"
}
if (!params.bamfile){
  exit 1, "--bamfile not set!"
}

if (!params.publish_mode){
  params.publish_mode = 'copy'   // or 'symlink'
}

if (!(params.publish_mode == 'copy' || params.publish_mode == 'symlinnk')){
  exit 1, "--publish_mode must be either copy or symlink"
}



process mkfastq {
  publishDir "${params.outdir}", mode: "${params.publish_mode}"

  input:
  file 'input.bam' from file(params.bamfile)

  output:
  // // this will emit a single item (a list) on r1files and a single list on r2files
  // file 'out/gemgroup001/bamtofastq*_R1_*.fastq.gz' into r1files
  // file 'out/gemgroup001/bamtofastq*_R2_*.fastq.gz' into r2files

  // this will emit multiple items!
  file 'out/*/bamtofastq*_R1_*.fastq.gz' into r1files mode flatten
  file 'out/*/bamtofastq*_R2_*.fastq.gz' into r2files mode flatten

  script:
  """
  ${BAM2FASTQ_BINARY} input.bam . --reads-per-fastq 50000000 --nthreads ${params.cpus}
  """
}
