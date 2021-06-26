#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ===================================
    nf-sarseq2
    ===================================

    Example Commandline:

    /groups/bioinfo/shared/collaborations/ages.grp/nf/1.0/run.sh  \
                --inputdir /path/fastqdir \
                --bc2sample metadata.xlsx \
                -resume

    ---------------------------------- Required arguments ----------------------------------------

    --inputdir  Path to raw data folder. It should include:
        *R1* file forward reads (typically named Undetermined..._R1_001.fastq.gz)
        *R2* file reverse reads (typically named Undetermined..._R2_001.fastq.gz)
        *I1* i5 indices as separate fastq.gz file (typically named Undetermined..._I1_001.fastq.gz)
        *I2* i7 indices as separate fastq.gz file (typically named Undetermined..._I2_001.fastq.gz)
       or:
    --redo Path to demuxed dir from previous run


    ---------------------------------- Optional arguments ----------------------------------------

    --outdir        basename of the output directory where the results will be saved [default: $params.outdir]

    Options used for Step1- demultiplexing of fastq to samples in SARSeq:

    --W           [file]    text file with well barcodes [default: $params.W]
    --P           [file]    text file with plate barcodes [default: $params.P]
    --platemin    [int]     min plate index to report [default: $params.platemin]
    --platemax    [int]     max plate index to report [default: $params.platemax]
    --MM          [0/1]     mismatches in plate and well indices (0 or 1) [default: $params.MM]

    Options used for Step2- assigning reads of a sample to amplicons (read2tile):
    --primerk    how many residues of primer to use for assigning reads to tile [default: $params.primerk]
    --A          tsv file containing amplicons [default: $params.A]

    Options used for SARSeq2 analysis:
    --bc2sample   [file]    bc2sample file with format EB\tPrimerBC\tSampleID\trest
    --bqfilter    [int]     mpileup bqfilter only applied to per tile analysis (not for ivar- see manual) [default: $params.bqfilter]


    """.stripIndent()
}



min_aln_reads = params.sample_minreads


////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

params.help = false
if (params.help) exit 0, helpMessage()


////////////////////////////////////////////////////
/* --          GET/VALIDATE INPUTS             -- */
////////////////////////////////////////////////////


ch_multiqc_config = file(params.multiqc_config, checkIfExists: true)

skipped_poor_alignment = []
def check_log(logs) {
    def prop_paired = 0;
    logname = logs.getBaseName() - 'flagstat'
    logs.eachLine { line ->
    log.info "$logname"
    log.info line
        if ((matcher = line =~ /^([\d\.]+) .* properly paired/)) {
            prop_paired = matcher[0][1]
        }
    }
    if(prop_paired.toInteger() < min_aln_reads.toInteger() ){
        log.info "####################  IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($logname)    >> ${prop_paired} prop paired <<"
        skipped_poor_alignment << logname
        return false
    } else {
        return true
    }
}


////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary = [:]
if(params.inputdir){
    summary['inputdir'] = params.inputdir
}
if(params.redo){
    summary['demuxdir'] = params.redo
}
summary['column/well BC file'] = params.W
summary['row/plate BC file'] =  params.P
summary['Min Plate Index'] = params.platemin
summary['Max Plate Index'] = params.platemax
summary['Max Mismatches'] = params.MM

summary['amplicon file'] = params.A
summary['primerk for tile matching'] = params.primerk

summary['Aln/Fq Filter (minreads;min prop paired)'] = params.sample_minreads
summary['min positional RC for cons'] = params.min_poscoverage
summary['bqfilter'] = params.bqfilter

summary['IVAR.Min Read Depth']        = params.ivar_min_coverage
summary['IVAR.Max Allele Freq']       = params.ivar_max_allele_freq

summary['Command line']  = workflow.commandLine
summary['Container']      = workflow.container

full_initial_summary = summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
full_initial_summary = full_initial_summary + "\n"
log.info full_initial_summary

full_initial_summary_file = file(params.outdirExt + "/" + "run_parameters.txt" )
full_initial_summary_file.text=full_initial_summary

if(params.bc2sample){
    ch_bc2sample = file(params.bc2sample, checkIfExists: true)
}else{
   ch_bc2sample = Channel.from(file ("empty"))
}
ch_spikefa = file(params.spikefa, checkIfExists: true)
ch_genomefa = file(params.genomefa, checkIfExists: true)
ch_amplicons = file(params.A, checkIfExists: true)
ch_wellbc = file(params.W, checkIfExists: true)

parameter_summary =  Channel
    .fromPath(full_initial_summary_file)
parameter_summary.into{parameter_summary1;parameter_summary2}




////////////////////////////////////////////////////
/* --            GENERATE INDICES              -- */
////////////////////////////////////////////////////



process index_tiles {

    input:
     file amplicons from ch_amplicons

    output:
     file("indices") into ch_tile_indices

    script:
    """

    ${baseDir}/bin/tsv2index.sh ${amplicons}

    """
}

process index_genome {

    input:
    file("indices/sequence.fasta") from ch_genomefa

    output:
    file("indices") into ch_genome_indices

    script:
    """
    samtools faidx indices/sequence.fasta
    minimap2 -d indices/sequence.mmi indices/sequence.fasta
    """
}


process make_primerbed {

    input:
    file genomeseq from ch_genomefa
    file amplicons from ch_amplicons

    output:
    file("primer.bed") into ch_primerbed

    script:
    """
    #!/usr/bin/env Rscript
    library(Biostrings)
    library(dplyr)


    genomeF <- "${genomeseq}"
    primerF <- "${amplicons}"
    genome <- readDNAStringSet(genomeF)
    primer <- readr::read_tsv(primerF,col_names=F) %>% setNames(c("tilename","primFWD","primREV","tileseq"))
    primerBed <- NULL
    for(i in 1:nrow(primer)){
       print(i)
       tilename <- primer[i,]\$tilename
       primerFWD <- DNAString(primer[i,]\$primFWD)
       primerREV <- reverseComplement(DNAString(primer[i,]\$primREV))

       primerFWDmatch <- data.frame(unlist(vmatchPattern(primerFWD, genome))) %>%
                         mutate(strand="+",tile=tilename,score=60) %>%
                         dplyr::select(names,start,end,tile,score,strand)
       primerREVmatch <- data.frame(unlist(vmatchPattern(primerREV, genome)))%>%
                         mutate(strand="-",tile=tilename,score=60) %>%
                         dplyr::select(names,start,end,tile,score,strand)
       primermatch <- rbind(primerFWDmatch,primerREVmatch)

       primerBed <- rbind(primerBed,primermatch)

     }
     primerBed <- primerBed %>% distinct() %>% mutate(start=start-1)

     readr::write_tsv(primerBed,"primer.bed",col_names=F)
    """
}





////////////////////////////////////////////////////
/* --    SPLIT FASTQ AS in SARSeq              -- */
////////////////////////////////////////////////////




if(params.inputdir){

    Channel
    .from(
    file("${params.inputdir}/*R1*", checkIfExists: true, maxDepth: 1),
    file("${params.inputdir}/*R2*", checkIfExists: true, maxDepth: 1),
    file("${params.inputdir}/*I1*", checkIfExists: true, maxDepth: 1),
    file("${params.inputdir}/*I2*", checkIfExists: true, maxDepth: 1)
    )
    .collect()
    .ifEmpty { error "Cannot find any files matching R1 R2 I1 I2 in inputdir ${params.inputdir}/" }
    .set{ch_input}

    process demux {
      publishDir path: "${params.outdirExt}/QC", mode: 'copy',
      saveAs: {filename->
        if (filename.lastIndexOf(".fastq") > 0) null
        else filename
    }

    input:
    set file(r1),file(r2),file(i1),file(i2) from ch_input
    file platebc from platebc_file
    file wellbc from ch_wellbc


    output:
    file ("*fastq") into ch_dmx
    file ("*tsv") into ch_dmx_rmd

    script:
    """
    cp ${baseDir}/bin/SARSeq_fixedOFF_demux.sh .
    ./SARSeq_fixedOFF_demux.sh -1 ${r1} -2 ${r2} -3 ${i1} -4 ${i2} -M ${params.MM}  -t sarseq1  -i ${params.platemin} -a ${params.platemax} -I ${platebc} -W ${wellbc}
    rm not_demuxed.R*
    """
    }

}else if (params.redo){
    ch_dmx = Channel.fromPath("${params.redo}/*fastq")
    ch_dmx_rmd = Channel.fromPath("${params.redo}/QC__demux.tsv")
    file("${params.redo}/QC__demux.tsv").copyTo(params.outdirExt + "/QC/QC__demux.tsv")

}else{
    exit 1, helpMessage() + "Provide either --inputdir with the raw fastq- not demuxed; or --redo with demuxed data"

}




// process all samples with the same second barcode together in the following
ch_dmx.collect().flatten().map{ l ->
    [ l.getName().split("-")[1] - ".R12.fastq" ,l]
}.groupTuple().into{ch_aln2spike_set;ch_read2tile}






////////////////////////////////////////////////////
/* --           TILE BASED ANALYSIS            -- */
////////////////////////////////////////////////////





process read2tile {

    tag "$bcset"

    input:
    set val(bcset), file(reads) from ch_read2tile
    file amplicons from ch_amplicons


    output:
    file ("*tsv") optional true into ch_r2t_report
    file ("*fastq") optional true into ch_r2t_fastq

    script:
    """
    mkdir tmp

    for infile in $reads; do
       name=\$(basename \$infile | sed 's/.R12.fastq//')

       # only work of files with reads
       if [ \$(wc -l \$infile | cut -f 1 -d " ") -ge "${params.sample_minreads_fastq_paired}" ]; then

          ${baseDir}/bin/deinterleave_fastq.sh < \$infile tmp/\${name}.R1.fastq tmp/\${name}.R2.fastq
          awk 'NR%4==1' tmp/\${name}.R1.fastq  | cut -f 1 -d " " > names1
          awk 'NR%4==1' tmp/\${name}.R2.fastq  | cut -f 1 -d " " > names2
          diff --brief names1 names2 >/dev/null
          comp_value=\$?

          if [ \$comp_value -eq 1 ]
          then
              exit "names of deinterleaved files could be different. is something going wrong?"
          fi

          cp ${baseDir}/bin/SARSeq_read2tiles.sh .

         ./SARSeq_read2tiles.sh -1 tmp/\${name}.R1.fastq -2 tmp/\${name}.R2.fastq  -t \${name}  -A ${amplicons} -P ${params.primerk}
      fi

    done

    """



}


ch_r2t_fastq
    .collect()
    .flatten()
    .map{ l ->
    [  l.getName().split("-")[1] - ".R12.fastq" ,l]
}.groupTuple().set{ch_r2t_fastq_set}


process aln2tile {

    tag "${setname}"

    input:
    set val(setname), file(reads) from ch_r2t_fastq_set
    file index from ch_tile_indices

    output:
    set file("*bam"), file("*flagstat") optional true into ch_aln2tile_bam_log
    file("*flagstat") optional true into ch_aln2tile_multiqc


    script:
    """
    for infile in $reads; do

      if [ \$(wc -l \$infile | cut -f 1 -d " ") -ge "${params.sample_minreads_fastq_paired}" ]; then
        name=\$(basename \$infile | sed 's/.R12.fastq//')

        tile=\$(echo \$infile | sed 's/__.*\$//')

        minimap2 -ax sr -t $task.cpus indices/\$tile/sequence.mmi \$infile | \\
          samtools view -@ $task.cpus -b -h -F 0x0100 - | samtools sort -o \$name.bam -

        samtools index \$name.bam

        samtools flagstat -@ $task.cpus \$name.bam > \$name.flagstat
      fi
    done
    """

}



process aln2tile_multiqc {

    publishDir path: "${params.outdirExt}/QC/multiqc_aln2tile", mode: 'copy'

    input:
    file("*") from ch_aln2tile_multiqc.collect()
    file multiqc_config from ch_multiqc_config

    output:
    file "*multiqc*" into ch_aln2tile_multiqc_report

    script:
    """
    multiqc . -m samtools -n multiqc_aln2tile --config $multiqc_config
    """

}




///////////////// Filter low mapping files per sample


ch_aln2tile_bam_log.collect().flatten().map{ l ->
    [ l.getName() - ".bam" - ".flagstat"  ,l]
}.groupTuple().map{ l ->
    [ l[1][0],l[1][1] ]
    }.set{ch_aln2tile_bam_log_per_sample}

// alignment filter not used currently
ch_aln2tile_bam_log_per_sample
//    .filter { bams,logs -> check_log(logs) }
    .collect()
    .flatten()
    .filter( ~/.*bam/ )
    .map{ l ->
    [ l.getName().split("-")[1] - ".bam" ,l]
}.groupTuple().set{ch_aln2tile_bam_set}



process aln2tile_mpile {

    tag "$setname"

    publishDir path: "${params.outdirExt}/aln2tile", mode: 'copy',
    saveAs: {filename->
      if (filename.lastIndexOf(".mpileup") > 0) "mpileup/$filename"
      else if (filename.lastIndexOf(".csv") > 0) "mpileup2csv/$filename"
      else null
    }

    input:
    set val(setname), file (bams) from ch_aln2tile_bam_set
    file index from ch_tile_indices

    output:
    file("*fractions.csv") optional true into ch_aln2tile_fractions_consensus

    script:
    bqfilter = params.bqfilter > 0 ? "--min-BQ ${params.bqfilter}" : '--no-BAQ'

    """
    for infile in $bams; do
      name=\$(basename \$infile | sed 's/.bam//')
      tile=\$(echo \$infile | sed 's/__.*\$//')

      # filter bam
      samtools view -@ $task.cpus -b -f 2 -F4 -F 0x0100  \$infile | samtools sort -o \$name.filt.tmp.bam -
      samtools mpileup -aa --max-depth 0 $bqfilter --min-MQ 1  -f indices/\$tile/sequence.fasta   \$name.filt.tmp.bam > \$name.mpileup.txt
      readstomper.pl \$name.mpileup.txt > \$name.mpileup.fractions.csv

    done

    """
}




////////////////////////////////////////////////////
/* --    OPTIONAL IVAR BASED ANALYSIS     -- */
////////////////////////////////////////////////////




if(params.do_ivar){

  process aln2spike {

      tag "$bcset"

      input:
      set val(bcset), file(reads) from ch_aln2spike_set
      file index from ch_genome_indices

      output:
      set file("*bam"), file("*flagstat") optional true into ch_aln2spike_bam_log
      file("*flagstat") optional true into ch_aln2spike_multiqc

      script:
      """
      for infile in $reads; do
        if [ \$(wc -l \$infile | cut -f 1 -d " ") -ge "${params.sample_minreads_fastq_paired}" ]; then
          name=\$(basename \$infile | sed 's/.R12.fastq//')
          minimap2 -ax sr -t $task.cpus indices/sequence.mmi \$infile | \\
            samtools view -@ $task.cpus -b -h -F 0x0100 - | samtools sort -o \$name.bam -
          samtools index \$name.bam
          samtools flagstat -@ $task.cpus \$name.bam > \$name.flagstat
        fi
      done
      """
  }


  process aln2spike_multiqc {

      publishDir path: "${params.outdirExt}/QC/multiqc_aln2spike", mode: 'copy'

      input:
      file "*" from ch_aln2spike_multiqc.collect()
      file multiqc_config from ch_multiqc_config

      output:
      file "*multiqc*" into ch_aln2spike_multiqc_report

      script:
      """
      multiqc . -m samtools -n multiqc_aln2spike --config $multiqc_config
      """
  }


   ///////////////// Filter low mapping files per sample

  ch_aln2spike_bam_log.collect().flatten().map{ l ->
      [ l.getName() - ".bam" - ".flagstat"  ,l]
  }.groupTuple().map{ l ->
      [ l[1][0],l[1][1] ]
  }.set{ch_aln2spike_bam_log_per_sample}

  // turned off the alignment filter- does not pay off
  ch_aln2spike_bam_log_per_sample
  //    .filter { bams,logs -> check_log(logs) }
    .collect()
    .flatten()
    .filter( ~/.*bam/ )
    .map{ l ->
    [ l.getName().split("-")[1] - ".bam" ,l]
  }.groupTuple().set{ch_aln2spike_bam_set}



///////////////// IVAR procedure

  process aln2spike_ivar_mpile {

      tag "$bcset"

      publishDir path: "${params.outdirExt}/aln2spike", mode: 'copy',
      saveAs: {filename->
        if (filename.lastIndexOf(".mpileup") > 0) "mpileup/$filename"
        else if (filename.lastIndexOf(".csv") > 0) "mpileup2csv/$filename"
        else null
      }

      input:
      set val(bcset), file (bams) from ch_aln2spike_bam_set
      file index  from ch_genome_indices
      file primerbed from ch_primerbed

      output:
      file("*fa") optional true into ch_ivar_consensus
      file("*fractions.csv") optional true into ch_aln2spike_fractions_consensus

      script:
      """
      for infile in $bams; do
        name=\$(basename \$infile | sed 's/.bam//')

        # filter bam
        samtools view -@ $task.cpus -b -f 2 -F4 -F 0x0100  \$infile | samtools sort -o \$name.filt.tmp.bam -

        #  primer trim; excluding reads with no primers
        ivar trim -b ${primerbed} -p \$name.trim.tmp.bam -i \$name.filt.tmp.bam
        samtools sort -@ $task.cpus -o \$name.trim.sort.tmp.bam -T \$name  \$name.trim.tmp.bam
        samtools index \$name.trim.sort.tmp.bam

        # mpileup
        samtools mpileup -aa --max-depth 0 --no-BAQ --min-MQ 1 -f indices/sequence.fasta \$name.trim.sort.tmp.bam > \$name.mpileup

        # mpileup to table
        readstomper.pl \$name.mpileup > \$name.mpileup.fractions.csv

        # mpileup to consensus
        cat \$name.mpileup | ivar consensus  -t ${params.ivar_max_allele_freq} -m ${params.ivar_min_coverage}  -p \$name.ivar.consensus.fa

        header=\$(head -n1 \$name.ivar.consensus.fa | sed 's/>//g')
        sed -i "s/\${header}/\$name/g" \$name.ivar.consensus.fa

        rm *tmp.bam
      done
      """
  }


  process consensus_pangolin {
      publishDir path: "${params.outdirExt}/ivar.consensus.pangolin", mode: 'copy'

      input:
      file "*" from ch_ivar_consensus.collect()

      output:
      file "lineage_report.csv" into ch_lineage
      file "ivar.consensus.fasta"

      script:
      """
      # remove end Ns- for only N sequences will leave a single N
      cat *fa | perl -pe 's/^N+/N/;s/N+\$/N/' > ivar.consensus.fasta
      pangolin -t $task.cpus --min-length 1000 ivar.consensus.fasta
      """
  }
}else{
    ch_lineage=Channel.from(file("empty"))
    ch_aln2spike_multiqc_report=Channel.from(false)
}


process summarize {


    publishDir path: "${params.outdirExt}/", mode: 'copy'


    input:
    file "QC/multiqc_aln2spike/*" from ch_aln2spike_multiqc_report
    file "QC/multiqc_aln2tile/*" from ch_aln2tile_multiqc_report
    file "ReadInTile_mpileup/*" from ch_aln2tile_fractions_consensus.collect()
    file "indices" from ch_tile_indices
    file spikeFa from ch_spikefa
    file lineageF from ch_lineage
    file sampleinfoF from ch_bc2sample


    output:
     file "*xlsx"
     file "*tsv" optional true
     file "*fa" optional true
     file "*rda" optional true


    script:
    def outdir="results"
    """
    #!/usr/bin/env Rscript
    spikeFa="${spikeFa}"
    lineageF="${lineageF}"
    sampleinfoF="${sampleinfoF}"
    minposcoverage=as.numeric(${params.min_poscoverage})



    source("${baseDir}/bin/functions.R")
    source("${baseDir}/bin/summarize_aln2tile.R")



    """
}





process checkContainer {

    publishDir path: "${params.outdirExt}/sw", mode: 'copy'
    cache false

    output:
    file "*txt"

    script:
    """
    set +u
    if [ ! -z \$SINGULARITY_CONTAINER ]; then
      conda list --export > containerPackagelist.txt
    fi
    """
}





workflow.onComplete {
    def subject = "[sarseq2] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[sarseq2] Failed: $workflow.runName"
    }
    log.info subject
    log.info "[sarseq2] Pipeline Complete"
}
