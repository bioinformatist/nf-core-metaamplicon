#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/metaamplicon
========================================================================================
 nf-core/metaamplicon Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/metaamplicon
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info"""
    =======================================================
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__         }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'

     nf-core/metaamplicon v${workflow.manifest.version}
    =======================================================

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/metaamplicon --reads 'data/*.fastq.gz' -profile standard,docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: standard, conda, docker, singularity, test

    Options:
      --singleEnd                   Specifies that the input is single end reads
      --reference                   Specifies the reference alignments file of mothur
	  --taxonomy					Specifies the taxonomy alignments file of mothur

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

if( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
  // Check workDir/outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

// Stage config files
output_docs = file("$baseDir/docs/output.md")


if(params.singleEnd){
    Channel
        .fromPath(params.reads)
        .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
        .set { read_files }
} else {
    Channel
        .fromFilePairs(params.reads, size: -1)
        .ifEmpty { exit 1, "Cannot find any reads in directory: ${params.reads}\nNB: Path needs to be enclosed in quotes!" }
        .set { read_files }
}

Channel
	.fromPath(params.reference)
	.ifEmpty { exit 1, "No reference alignments file found!" }
	.set{ reference }

Channel
	.fromPath(params.taxonomy)
	.ifEmpty { exit 1, "No taxonomy alignments file found!" }
	.set{ taxonomy }


// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/metaamplicon v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/metaamplicon'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Reads']        = params.reads
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


/*
 * Parse software version numbers
 */
// process get_software_versions {

//     output:
//     file 'software_versions_mqc.yaml' into software_versions_yaml

//     script:
//     // TODO nf-core: Get all tools to print their version number here
//     """
//     echo $workflow.manifest.version > v_pipeline.txt
//     echo $workflow.nextflow.version > v_nextflow.txt
//     # fastqc --version > v_fastqc.txt
//     # multiqc --version > v_multiqc.txt
//     scrape_software_versions.py > software_versions_mqc.yaml
//     """
// }
   
// Use a temp dir for mothur's intermediate results
// TODO: add a parameter for choosing not keep this dir when finish
mothur_temp = file("${params.outdir}/mothur_temp")
mothur_temp.mkdirs()


// TODO: add a parameter for skipping this QC step
/*
 * STEP 1 - fastp
 */
if (params.singleEnd) {
    process fastp_single {
        publishDir "${params.outdir}/fastp", mode: 'copy'

        input:
        file reads from read_files

        output:
        file "*.{fastq,html,json}"
        // TODO: simplify below using Groovy structures
        file "*.fasta" into for_mothur0
        file "*.qual" into for_mothur1
        file "*.fasta" into for_group
        
        // Example: reads: DRR067825.fastq.gz; reads.baseName: DRR067825.fastq; reads.simpleName: DRR067825
        // See https://github.com/nextflow-io/nextflow/issues/278
        """
        fastp -i ${reads} -o ${reads.baseName} -h ${reads.simpleName}.html -j ${reads.simpleName}.json       
        mothur "#fastq.info(fastq=${reads.baseName})"
        """
    }
    // We don't know when the intermediate files will be published so publishDir should not be used here
    for_mothur0.subscribe { it.copyTo(mothur_temp) }
    for_mothur1.subscribe { it.copyTo(mothur_temp) }
} else {
    process fastp_paired {
        publishDir "${params.outdir}/fastp", mode: 'copy'

        input:
        set name, file(reads) from read_files

        output:
        file "*.{fastq,html,json}"
        // TODO: simplify below using Groovy structures
        file "*.fastq" into for_mothur0
        // Channel only for labelling this process is all finished
        file '*.html' into done  
        
        """
        fastp -i ${reads[0]} -I ${reads[1]} -o ${reads[0].baseName} -O ${reads[1].baseName} -h ${name}.html -j ${name}.json
        # vsearch of current version (called by mothur's chimera.vsearch() function) cut too long seq names,
        # which will cause "XXX is not in your count table. Please correct." error.
        # So it is needed for pre-cut seq names here, I choose remove first 11 characters now.
        # Update: This strategy will import new problems in down-stream analysis. Change vsearch to v2.8.0 temperarily.
        """
    }
    for_mothur0.subscribe { it.each { it.copyTo(mothur_temp) } }
}


/*
 * STEP 2 - mothur in-box
 */
if (params.singleEnd) {
    process mothur_in_box_single {
        input:
        file fasta from for_group.collect()
        file reference
        file taxonomy

        // output:
        // file 'final.opti_mcc.0.03.biom' into 'final.biom'

        script:
        """
        #!/usr/bin/env mothur
        # To make group file for single-end data
        make.group(fasta=${fasta.join('-')}, groups=${fasta.collect { it.simpleName }.join('-')})
        system(cp groups $workflow.launchDir/${params.outdir}/mothur_temp/merge.groups)
        system(cat *.fasta > stability.trim.contigs.fasta)
        system(cp stability.trim.contigs.fasta $workflow.launchDir/${params.outdir}/mothur_temp)
        set.dir(input=$workflow.launchDir/${params.outdir}/mothur_temp, output=$workflow.launchDir/${params.outdir}/mothur_temp)
        summary.seqs(fasta=stability.trim.contigs.fasta)
        screen.seqs(fasta=stability.trim.contigs.fasta, group=merge.groups, summary=stability.trim.contigs.summary, optimize=maxlength, maxambig=0, criteria=97.5)
        summary.seqs()
        unique.seqs(fasta=stability.trim.contigs.good.fasta)
        count.seqs(name=stability.trim.contigs.good.names, group=merge.good.groups)
        summary.seqs(count=stability.trim.contigs.good.count_table)
        align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=$workflow.launchDir/${reference})
        summary.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table)
        screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, summary=stability.trim.contigs.good.unique.summary, optimize=start-end, criteria=97.5, maxhomop=8)
        summary.seqs(fasta=current, count=current)
        filter.seqs(fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=.)
        unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.good.count_table)
        pre.cluster(fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filter.count_table, diffs=3)
        chimera.vsearch(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)
        remove.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)
        summary.seqs(fasta=current, count=current)
        classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=$workflow.launchDir/${reference}, taxonomy=$workflow.launchDir/${taxonomy}, cutoff=80, iters=1000)
        remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.gg.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Eukaryota)
        summary.tax(taxonomy=current, count=current)
        cluster.split(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.gg.wang.pick.taxonomy, splitmethod=classify, taxlevel=4, cutoff=0.03)
        make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)
        classify.otu(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.gg.wang.pick.taxonomy, label=0.03)
        get.oturep(column=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, name=stability.trim.contigs.good.names, list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, fasta=stability.trim.contigs.good.unique.fasta)
        rename.file(taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy, shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared, prefix=final)
        make.biom(shared=final.opti_mcc.shared, constaxonomy=final.taxonomy)
        """
    }
} else {
    process mothur_in_box_paired {     
        input:
        // Files are not really used here. Only to ensure that this process runs after the previous one
        file done from done.collect()
        file reference
        file taxonomy

        """
        #!/usr/bin/env mothur
        make.file(inputdir=$workflow.launchDir/${params.outdir}/mothur_temp)
        make.contigs(file=stability.files, checkorient=t, processors=96)
        summary.seqs(fasta=stability.trim.contigs.fasta)
        screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, summary=stability.trim.contigs.summary, optimize=maxlength, maxambig=0, criteria=97.5)
        summary.seqs()
        unique.seqs(fasta=stability.trim.contigs.good.fasta)
        count.seqs(name=stability.trim.contigs.good.names, group=stability.contigs.good.groups)
        summary.seqs(count=stability.trim.contigs.good.count_table)
        align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=$workflow.launchDir/${reference})
        summary.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table)
        screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, summary=stability.trim.contigs.good.unique.summary, optimize=start-end, criteria=97.5, maxhomop=8)
        summary.seqs(fasta=current, count=current)
        filter.seqs(fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=.)
        unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.good.count_table)
        pre.cluster(fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filter.count_table, diffs=3)
        chimera.vsearch(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)
        remove.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)
        summary.seqs(fasta=current, count=current)
        classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=$workflow.launchDir/${reference}, taxonomy=$workflow.launchDir/${taxonomy}, cutoff=80, iters=1000)
        remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.gg.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Eukaryota)
        summary.tax(taxonomy=current, count=current)
        cluster.split(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.gg.wang.pick.taxonomy, splitmethod=classify, taxlevel=4, cutoff=0.03)
        make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)
        classify.otu(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.gg.wang.pick.taxonomy, label=0.03)
        get.oturep(column=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, name=stability.trim.contigs.good.names, list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, fasta=stability.trim.contigs.good.unique.fasta)
        rename.file(taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy, shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared, prefix=final)
        make.biom(shared=final.opti_mcc.shared, constaxonomy=final.taxonomy)
        """
    }
}

// Temperarily cancel this process due to dependencies error:
// /share/apps/anaconda3/lib/R/bin/exec/R: /lib64/libstdc++.so.6: version `CXXABI_1.3.8' not found (required by /share/apps/anaconda3/lib/R/bin/exec/../../lib/../../libicuuc.so.58)
// TODO: fix it in conda
/*
 * Step - Output Description HTML
 */
// process output_documentation {
//     publishDir "${params.outdir}/Documentation", mode: 'copy'

//     input:
//     file output_docs

//     output:
//     file "results_description.html"

//     script:
//     """
//     markdown_to_html.r $output_docs results_description.html
//     """
// }


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/metaamplicon] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/metaamplicon] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/metaamplicon] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/metaamplicon] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/metaamplicon] Pipeline Complete"

}
