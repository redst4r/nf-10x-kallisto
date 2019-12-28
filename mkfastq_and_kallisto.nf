#!/usr/bin/env nextflow

// Some required setup:
// conda install -c bioconda bustools kallisto
// wget http://cf.10xgenomics.com/misc/bamtofastq-1.2.0 && chmod u+x bamtofastq-1.2.0


params.kallisto_index = '/home/mstrasse/resources/Homo_sapiens.GRCh38.cdna.all.idx'
params.kallisto_gene_map = '/home/mstrasse/resources/transcripts_to_genes.txt'

// V2 chemistr
// params.chemistry = '10xv2'
// params.barcode_whitelist = '/home/mstrasse/resources/10xv2_whitelist.txt'

// v3 chem
params.chemistry = '10xv3'
params.barcode_whitelist = '/home/mstrasse/resources/3M-february-2018.txt'

params.bustools_correct = true
params.cpus = 8
params.mem = '30G'


Channel
    .fromPath(params.kallisto_index)
    .ifEmpty { exit 1, "kallisto index not found: ${params.kallisto_index}" }
    .set {kallisto_index}

Channel
        .fromPath(params.kallisto_gene_map)
        .ifEmpty { exit 1, "kallisto gene map not found: ${params.kallisto_gene_map}" }
        .set {kallisto_gene_map}

Channel.fromPath(params.barcode_whitelist)
       .ifEmpty{ exit 1, "Cannot find ${params.type} barcode whitelist: $barcode_filename" }
       .set{barcode_whitelist_kallisto }



process mkfastq {
  // publishDir "./pipe_out", mode: 'copy'

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
  /home/mstrasse/bamtofastq-1.2.0 input.bam out --reads-per-fastq 50000000 --nthreads ${params.cpus}
  """
}

/*
combined emits pairs of (R1,R2), (S1,S2), ...
we need to put that into a single item (R1,R2,S1,S2)
*/
combined = r1files.merge(r2files)
combined_flat = combined.flatten().collect()
// combined_flat.subscribe onNext: { println it.name }, onComplete: { println 'Done' }
// r1files.subscribe onNext: { println it.name }, onComplete: { println 'Done' }


 process kallisto {
     // label 'low_memory'
     // publishDir "${params.outdir}/kallisto/raw_bus", mode: 'copy'

     input:
     // set val(name), file(reads) from read_files_kallisto
     file(reads) from combined_flat
     file index from kallisto_index.collect()

     output:
     file "bus_output" into kallisto_bus_to_sort
     file "kallisto.log" into kallisto_log_for_multiqc

     script:
     """
     echo $index
     kallisto bus \\
         -i $index \\
         -o bus_output/ \\
         -x ${params.chemistry} \\
         -t ${params.cpus} \\
         $reads | tee kallisto.log
     """
 }
 
 /*
 * Run BUSTools Correct / Sort on Kallisto Output
 */
 process bustools_correct_sort{
     // tag "$bus"
     // publishDir "${params.outdir}/kallisto/sort_bus", mode: 'copy'

     input:
     file bus from kallisto_bus_to_sort
     file whitelist from barcode_whitelist_kallisto

     output:
     file bus into (kallisto_corrected_sort_to_count, kallisto_corrected_sort_to_metrics)


     script:
     if(params.bustools_correct) {
       correct = "bustools correct -w $whitelist -o ${bus}/output.corrected.bus ${bus}/output.bus"
       sort_file = "${bus}/output.corrected.bus"
     } else {
       correct = ""
       sort_file = "${bus}/output.bus"
     }
     """
     $correct
     mkdir -p tmp
     bustools sort -T tmp/ -t ${params.cpus} -m ${params.mem} -o ${bus}/output.corrected.sort.bus $sort_file
     """
 }

 /*
 * Run BUSTools count on sorted/corrected output from Kallisto|Bus
 */
 process bustools_count{
     // tag "$bus"
     publishDir "${params.outdir}/kallisto/bustools_counts", mode: "copy"

     input:
     file bus from kallisto_corrected_sort_to_count
     file t2g from kallisto_gene_map.collect()

     output:
     file "${bus}_eqcount"
     file "${bus}_genecount"

     script:
     """
     mkdir -p ${bus}_eqcount
     mkdir -p ${bus}_genecount
     bustools count -o ${bus}_eqcount/tcc -g $t2g -e ${bus}/matrix.ec -t ${bus}/transcripts.txt ${bus}/output.corrected.sort.bus
     bustools count -o ${bus}_genecount/gene -g $t2g -e ${bus}/matrix.ec -t ${bus}/transcripts.txt --genecounts ${bus}/output.corrected.sort.bus
     """
 }

 /*
 * Run Bustools inspect
 */

 process bustools_inspect{
     // tag "$bus"
     publishDir "${params.outdir}/kallisto/bustools_metrics", mode: "copy"

     input:
     file bus from kallisto_corrected_sort_to_metrics

     output:
     file "${bus}.json"

     script:
     """
     bustools inspect -o ${bus}.json ${bus}/output.corrected.sort.bus
     """
 }
