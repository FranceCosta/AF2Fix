#!/usr/bin/env nextflow

params.pfamuser = "admin"     // user name for MySQL pfam database
params.pfampassword = "" // password for MySQL pfam database
params.pfamhost = "mysql-pfam-rel" // hostname for MySQL database
params.pfamport = "" // port for MySQL database
params.proteins_table = "${workflow.projectDir}/input/protein_table_solenoids.csv" // input protein table to run pipeline on customized proteins. Run with -entry improve
params.outdir  = "${workflow.projectDir}/output"
//params.afdatabase = "/hps/nobackup/agb/interpro/mblum/alphafold/pfam-alphafold.sqlite.old" // db with per-domain af distribution
params.afdatabase = "${workflow.projectDir}/assets/pfam-alphafold.sqlite"
params.domains_to_consider = 50 // Number of protein domains to include in the pipeline
params.proteins_per_class = 10 // Number of protein candidates to consider per class of AlphaFold confidence (<50, 50-70, >70 plddt). 10 more will be included and discarded by later script
params.colabfold_database =  "${workflow.projectDir}/assets/" //"/nfs/research/agb/research/francesco/software/localcolabfold/databases/"
params.colabfold_cache = "${workflow.projectDir}/tmp/"
params.test_msa = "${workflow.projectDir}/test"

// use module aliasing to call the same process multiple times
include { RANK as RANK_NO_TEMPLATE } from "${workflow.projectDir}/modules/af2rank" params(outdir: params.outdir, projectDir: "${workflow.projectDir}")
include { RANK as RANK_TEMPLATE_MSA } from "${workflow.projectDir}/modules/af2rank" params(outdir: params.outdir, projectDir: "${workflow.projectDir}")
include { RANK as RANK_TEMPLATE_SEQ } from "${workflow.projectDir}/modules/af2rank" params(outdir: params.outdir, projectDir: "${workflow.projectDir}")

include { XHPI as MOLPROB_NO_TEMPLATE } from "${workflow.projectDir}/modules/molprobity" params(outdir: params.outdir, projectDir: "${workflow.projectDir}")
include { XHPI as MOLPROB_TEMPLATE_MSA } from "${workflow.projectDir}/modules/molprobity" params(outdir: params.outdir, projectDir: "${workflow.projectDir}")
include { XHPI as MOLPROB_TEMPLATE_SEQ } from "${workflow.projectDir}/modules/molprobity" params(outdir: params.outdir, projectDir: "${workflow.projectDir}")

process DOWNLOAD_COLABFOLD_WEIGHTS {

    container 'docker://athbaltzis/colabfold_proteinfold:1.1.0'
    time '2h'
    tag 'download_colabfold_weights'
    publishDir "${params.colabfold_cache}/", mode: 'copy', overwrite: true
    errorStrategy 'retry' // tries a second time
        
    output:
    path 'params*', emit: weights
    
    """
    download_colabfold_weights.py
    """
    /*
    // bug prone code:
    // do not download if already existing
    script:
    if( file("${params.colabfold_cache}/params/download_finished.txt").exists() == true )
        """
        touch download_finished.txt
        """
    
    else if( file("${params.colabfold_cache}/params/download_finished.txt").exists() == false )
        """
        download_colabfold_weights.py
        """
    */
}

process QUERY_DB {
    label 'python38'
    memory '500MB'
    tag 'query-db'
    publishDir "$params.outdir/DB_query", mode: 'copy'
  
    output:
    path 'proteins_per_domain.csv', emit: proteins_per_domain
    
    """
    query-db_v3.py --database ${params.afdatabase} --number_to_consider ${params.domains_to_consider} --proteins_per_class ${params.proteins_per_class}
    """
}

process QUERY_DB_WITH_INPUT {
    label 'python38'
    memory '500MB'
    tag 'query-db'
    publishDir "$params.outdir/DB_query", mode: 'copy'
  
    output:
    path 'proteins_per_domain.csv', emit: proteins_per_domain
    
    """
    query-db_input.py --database ${params.afdatabase} --max_proteins ${params.proteins_per_class} --proteins_table ${params.proteins_table}
    """
}

process GET_SEQUENCES {
    
    memory '500MB'
    tag 'get_sequences'
    publishDir "$params.outdir/DB_query", mode: 'copy'

    input:
    path proteins_per_domain

    output:
    path 'target_proteins.tsv', emit: proteins_info
    
    """
    get_sequences.sh -i ${proteins_per_domain} -o target_proteins.tsv -n ${params.proteins_per_class} -u ${params.pfamuser} -p ${params.pfampassword} -h ${params.pfamhost} -P ${params.pfamport}
    """
}

process DOWNLOAD_DOMAIN_LENGTHS {

    time '2h'
    tag 'download_domain_lwnghts'
    publishDir "${params.outdir}/DB_query/", mode: 'copy'
    
    input:
    path proteins_info
    
    output: 
    path 'domain_lengths.csv'
    
    """
    domain_lenghts.sh -i ${proteins_info} -u ${params.pfamuser} -p ${params.pfampassword} -h ${params.pfamhost} -P ${params.pfamport}
    """

}

process GET_MSAS {
   
   container 'docker://athbaltzis/colabfold_proteinfold:1.1.0'
   memory '180GB' // normally 180
   tag 'get_msas'
   publishDir "$params.outdir/AF2/MSAs", mode: 'copy'
   
   input:
   path proteins_info
   
   output:
   path ("PF*"), emit: path_to_msa

   """
   generate_MSAs.sh -i ${proteins_info} -d ${params.colabfold_database}
   """

}

process GET_TEST_MSA {
    memory '300MB'
    //publishDir "$params.outdir/MSAs", mode: 'copy', pattern: 'PF*/*.a3m' // this does not work
    //publishDir "$params.outdir/MSAs", mode: 'copy' // this works
    publishDir "$params.outdir/AF2/MSAs", mode: 'copy', pattern: "PF*" //this works
    //publishDir "$params.outdir/MSAs", mode: 'copy', pattern: "PF*/" // this does not work
    //publishDir "$params.outdir/MSAs", mode: 'copy', pattern: "*/*.a3m" // this does not work


    output:
    path ('PF*'),       emit: path_to_msa
    path ('PF*/*.a3m'), emit: path_to_msa_file

    """
    cp -r ${params.test_msa}/* ./
    """
}

process PREDICT_NO_TEMPLATE {
   
    container 'docker://athbaltzis/colabfold_proteinfold:1.1.0'
    tag 'predict_no_template'
    memory '10G'
    clusterOptions '--gres=gpu:a100:1'
    containerOptions '--nv'
    publishDir "$params.outdir/AF2/no_template/${domain}", mode: 'copy', pattern: '*_relaxed_rank_001_*.pdb'
    publishDir "$params.outdir/AF2/no_template/${domain}", mode: 'copy', pattern: '*_scores_rank_001_*.json'
    publishDir "$params.outdir/AF2/no_template/${domain}", mode: 'copy', pattern: 'log.txt' 

    input:
    path path_to_msa // use val if path not working
    path weights
    
    output:
    path '*_relaxed_rank_001_*.pdb', emit: pdb_structure
    path '*_scores_rank_001_*.json'
    path 'log.txt'
    tuple path('*_relaxed_rank_001_*.pdb'), val(domain), emit: pdb_structures_domain
    tuple path("./"), path(path_to_msa), emit: output_and_msa

    script:
    domain = path_to_msa.baseName

    """
    colabfold_batch --data ./ --num-recycle 5 --use-gpu-relax --num-relax 1 --model-type auto\
                    ${path_to_msa} ./
    """

}

process GET_TEMPLATES {
    
    container 'docker://athbaltzis/colabfold_proteinfold:1.1.0'
    tag 'get_templates'
    publishDir "$params.outdir/AF2/templates/$domain/", mode: 'copy', pattern: '*.cif'
    publishDir "$params.outdir/AF2/templates/$domain/", mode: 'copy', pattern: '.mapper.json'

    input:
    tuple path(output_folder, stageAs: 'af2results'), path(path_to_msa)

    output:
    path '*.cif', optional: true
    path '.mapper.json', optional: true
    tuple path('*.cif'), path(path_to_msa), emit: templates_msa, optional: true

    script:
    domain = path_to_msa.baseName

    """
    PDBtoMMCIF.py --AFres_dir ./af2results --output_dir ./
    """
}

process PREDICT_TEMPLATE_MSA {

    container 'docker://athbaltzis/colabfold_proteinfold:1.1.0'
    tag 'predict_no_template'
    memory '10G'
    clusterOptions '--gres=gpu:a100:1'
    containerOptions '--nv'
    publishDir "$params.outdir/AF2/template_MSA/$domain", mode: 'copy'
    errorStrategy "ignore"

    input:
    tuple path(template), path(msa)
    path weights
    
    output:
    path '*_relaxed_rank_001_*.pdb', emit: pdb_structure, optional: true
    path '*_scores_rank_001_*.json', optional: true
    path 'log.txt', optional: true
    tuple path('*_relaxed_rank_001_*.pdb'), val(domain), emit: pdb_structures_domain, optional: true

    script:
    domain = msa.baseName

    """
    if [ "\$(ls -A *.cif)" ]; then
        colabfold_batch --data ./ --templates --custom-template-path ./\
                    --use-gpu-relax --num-relax 1  --num-recycle 5 --model-type auto ${msa} ./
	else
        echo "${msa} has no templates"
	fi
    """

}

process PREDICT_TEMPLATE_SEQ {

    container 'docker://athbaltzis/colabfold_proteinfold:1.1.0'
    tag 'predict_no_template'
    memory '15G'
    clusterOptions '--gres=gpu:a100:1'
    containerOptions '--nv'
    publishDir "$params.outdir/AF2/template_single_seq/$domain", mode: 'copy'
    errorStrategy "ignore"

    input:
    tuple path(templates_path), path(msa)
    path weights

    output:
    path '*_relaxed_rank_001_*.pdb', emit: pdb_structure, optional: true
    path '*_scores_rank_001_*.json', optional: true
    path 'log.txt', optional: true
    tuple path('*_relaxed_rank_001_*.pdb'), val(domain), emit: pdb_structures_domain, optional: true

    script:
    domain = msa.baseName
        
    """
    if [ "\$(ls -A *.cif)" ]; then
        
        # get sequences

        touch sequences.fa
        for msa in \$( ls ${msa}/*.a3m )
        do
            echo '>'\$(basename \$msa .a3m ) >> sequences.fa
            echo \$(head \$msa -n 2 | tail -n 1) >> sequences.fa
        done

        colabfold_batch --data ./ --templates --custom-template-path ./ --msa-mode single_sequence \
                    --use-gpu-relax --num-relax 1  --num-recycle 5 --model-type auto sequences.fa ./
	else
        echo "${msa} has no templates"
	fi
    """
        
}

process COLLECT_RESULTS {
   
    label 'python38'
    tag 'collect_results'
    publishDir "$params.outdir/summary_results/", mode: 'copy' // this finishes after pipeline has completed lol
    // saveAs

    input:
    val rank_mock
    val xhpi_mock
    val procheck_mock


    output:
    path "*.csv", emit: output_tables

    """
    # consider case in which params.outdir was given a relative path
    if [ ! -d "${params.outdir}" ]; then
        outdir="${workflow.projectDir}/${params.outdir}"
    else
        outdir="${params.outdir}"
    fi
    
    collect_results.py --AFres_dir \$outdir/AF2/ --AF2Rank_res_dir \$outdir/AF2Rank/\
                        --database_res_dir \$outdir/DB_query/ --xhpi_res_dir \$outdir/xhpi/\
                        --procheck_res_dir \$outdir/ramachandran
    """

}

workflow preparation {
    
    // this runs the classic pipeline without input
    main:
        DOWNLOAD_COLABFOLD_WEIGHTS()
        QUERY_DB()
        GET_SEQUENCES(QUERY_DB.out.proteins_per_domain)
        DOWNLOAD_DOMAIN_LENGTHS(GET_SEQUENCES.out.proteins_info)
    emit:
        weights       = DOWNLOAD_COLABFOLD_WEIGHTS.out.weights
        proteins_info = GET_SEQUENCES.out.proteins_info
}

workflow preparation_improve {
    
    // this takes a list of proteins and domains as input
    main:
        DOWNLOAD_COLABFOLD_WEIGHTS()
        QUERY_DB_WITH_INPUT()
        GET_SEQUENCES(QUERY_DB_WITH_INPUT.out.proteins_per_domain)
        DOWNLOAD_DOMAIN_LENGTHS(GET_SEQUENCES.out.proteins_info)
    emit:
        weights       = DOWNLOAD_COLABFOLD_WEIGHTS.out.weights
        proteins_info = GET_SEQUENCES.out.proteins_info
}

workflow predictions {
    
    take:
        weights
        proteins_info
    
    main:
        GET_MSAS(proteins_info)
        // 
        PREDICT_NO_TEMPLATE(GET_MSAS.out.path_to_msa.flatten(), weights)

        // create a set of unique paths to domain folders containing no template MSA AlphaFold results
        //PREDICT_NO_TEMPLATE.out.pdb_structure.collect().flatMap{ file -> "${file.toString().split('/')[-2]}" }.unique().set{ no_template_domains }
        GET_TEMPLATES(PREDICT_NO_TEMPLATE.out.output_and_msa)
        //GET_TEMPLATES.out.templates_msa.view { "$it" }
        
        PREDICT_TEMPLATE_MSA(GET_TEMPLATES.out.templates_msa, weights)
        PREDICT_TEMPLATE_SEQ(GET_TEMPLATES.out.templates_msa, weights)
        // AF2Rank
        // use .transpose() to have as input one .pdb at time
        RANK_NO_TEMPLATE(PREDICT_NO_TEMPLATE.out.pdb_structures_domain.transpose(), "no_template")
        RANK_TEMPLATE_MSA(PREDICT_TEMPLATE_MSA.out.pdb_structures_domain.transpose(), "template_MSA")
        RANK_TEMPLATE_SEQ(PREDICT_TEMPLATE_SEQ.out.pdb_structures_domain.transpose(), "template_single_seq")
        // Molprobity
        MOLPROB_NO_TEMPLATE(PREDICT_NO_TEMPLATE.out.pdb_structures_domain, "no_template")
        MOLPROB_TEMPLATE_MSA(PREDICT_TEMPLATE_MSA.out.pdb_structures_domain, "template_MSA")
        MOLPROB_TEMPLATE_SEQ(PREDICT_TEMPLATE_SEQ.out.pdb_structures_domain, "template_single_seq")

        // collect results after completion of previous processes
        // not very robust as files may not be in the publishDir by time 
        // COLLECT_RESULTS is run
        // alternatevly one could use symlinks to a scratch directory and
        // delete il later
        
        RANK_NO_TEMPLATE.out.csv_table.flatten()\
            .mix(RANK_TEMPLATE_MSA.out.csv_table.flatten(), RANK_TEMPLATE_SEQ.out.csv_table.flatten())\
            .collectFile(skip: 1, keepHeader: true)\
            .last()\
            .set{ rank_out }

        XHPI_NO_TEMPLATE.out.summary.flatten()\
            .mix(XHPI_TEMPLATE_MSA.out.summary.flatten(), XHPI_TEMPLATE_SEQ.out.summary.flatten())\
            .collectFile(skip: 1, keepHeader: true)\
            .last()\
            .set{ xhpi_out }

        RAMACHANDRAN_NO_TEMPLATE.out.summary.flatten()\
            .mix(RAMACHANDRAN_TEMPLATE_MSA.out.summary.flatten(), RAMACHANDRAN_TEMPLATE_SEQ.out.summary.flatten())\
            .collectFile(skip: 1, keepHeader: true)\
            .last()\
            .set{ procheck_out }
        
    emit:
        rank_out = rank_out
        xhpi_out = xhpi_out
        procheck_out = procheck_out

} 

workflow collect_results {
    take:
        rank_out
        xhpi_out
        procheck_out

    main:
        COLLECT_RESULTS(rank_out, xhpi_out, procheck_out)
}

workflow {
    
    // classic workflow to run AF2Fix on a set of PFAM domain proteins
    preparation()
    predictions(preparation.out.weights, preparation.out.proteins_info)
    // maybe this way collect_result() will run after files have been moved to 
    // outputdir by predictions()
    collect_results(predictions.out.rank_out, predictions.out.xhpi_out, predictions.out.procheck_out)
}

workflow improve {
    
    // workflow to improve specific proteins using PFAM googd qulaity structures as templates
    preparation_improve()
    predictions(preparation_improve.out.weights, preparation_improve.out.proteins_info)
    collect_results(predictions.out.rank_out, predictions.out.xhpi_out, predictions.out.procheck_out)

}

workflow test {
    // test if nextflow is copying things to publishDir
    //GET_TEST_MSA()
    //GET_TEST_MSA.out.path_to_msa_file\
    //    .flatten()
    //    .view { "MSA generated: $it"}

    // test correct tuple unpacking
    Channel.of( [['alpha', 'beta', 'delta'], 1], [['gamma', 'epsilon'], 2] )\
        .transpose()\
        .view { it }
}
