#!/usr/bin/env nextflow

// Some required setup:
// conda install -c bioconda bustools kallisto

// Example call
//
// nextflow run kallisto_from_fastq.nf --outdir /tmp/day2 --readsglob 'day_2_S1_L00*_R{1,2}_001.fastq.gz'

// checking input args
if (!params.chemistry){
  exit 1, "--chemistry not set!: 10xv3 or 10xv2"
}

if (!(params.chemistry == '10xv2' || params.chemistry == '10xv3')){
  exit 1, "--chemistrymust be either 10xv3 or 10xv2"
}

if (!params.outdir){
  exit 1, "--outdir not set!"
}

RESOURCE_DIR = '/home/michi/mounts/TB4drive/kallisto_resources'
// RESOURCE_DIR = '/run/media/michi/42506642-b470-4238-be14-bb0c303b3682/kallisto_resources'

//params.kallisto_index = RESOURCE_DIR + '/Homo_sapiens.GRCh38.cdna.all.idx'
//params.kallisto_gene_map = RESOURCE_DIR + '/transcripts_to_genes.txt'

//params.kallisto_index = '/home/michi/COVID/Homo_sapiens_COVID.idx'
//params.kallisto_gene_map = '/home/michi/COVID/transcripts_to_genes.txt'

if ((params.chemistry == '10xv2') && (!params.barcode_whitelist)){
  println "warning, using default v2 barcode list"
  params.barcode_whitelist = RESOURCE_DIR + '/10xv2_whitelist.txt'
}

if ((params.chemistry == '10xv3') && (!params.barcode_whitelist)){
  println "warning, using default v3 barcode list"
  params.barcode_whitelist = RESOURCE_DIR + '/3M-february-2018.txt'
}

params.bustools_correct = true
params.cpus = 4
params.mem = '10G'


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


 /*
 combined emits pairs of (R1,R2), (S1,S2), ...
 we need to put that into a single item (R1,R2,S1,S2)
 */
 Channel
    .fromFilePairs( params.readsglob )
    // fromFile pair emits a (name, [f1, f2]), but closures dont do the unpacking, so just take the second element which is [f1 ,f2]
    .flatMap {bla -> bla[1]}
    .collect()
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.readsglob}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\n" }
    .into { combined_flat }

// .subscribe onNext: { println it }, onComplete: { println 'Done' }


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
     publishDir "${params.outdir}/kallisto/sort_bus", mode: 'copy'

     input:
     file bus from kallisto_bus_to_sort
     file whitelist from barcode_whitelist_kallisto

     output:
     file bus into (kallisto_corrected_sort_to_count, kallisto_corrected_sort_to_metrics)


     script:
     if(params.bustools_correct) {
       correct = "bustools correct -w $whitelist -o ${bus}/output.corrected.bus ${bus}/output.bus"
       sort_file = "${bus}/output.corrected.bus"
       // output.corrected.sort.bus is the only important output, so
       // cleanup the intermediates
       cleanup = "rm ${bus}/output.corrected.bus && rm ${bus}/output.bus"
     } else {
       correct = ""
       sort_file = "${bus}/output.bus"
       cleanup = "rm ${bus}/output.bus"

     }
     """
     $correct
     mkdir -p tmp
     bustools sort -T tmp/ -t ${params.cpus} -m ${params.mem} -o ${bus}/output.corrected.sort.bus $sort_file
     $cleanup
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
