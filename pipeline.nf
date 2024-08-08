#!/usr/bin/env nextflow

params.pfamuser = "admin"     // user name for MySQL pfam database
params.pfampassword = "" // password for MySQL pfam database
params.pfamhost = "mysql-pfam-rel" // hostname for MySQL database
params.pfamport = "" // port for MySQL database
params.proteins_table = "${workflow.projectDir}/example/protein_table.csv" // input protein table to run pipeline on customized proteins. Run with -entry improve
params.outdir  = "${workflow.projectDir}/output"
params.afdatabase = "/hps/nobackup/agb/interpro/mblum/alphafold/pfam-alphafold.sqlite.old" // db with per-domain af distribution
//params.afdatabase = "${workflow.projectDir}/assets/pfam-alphafold.sqlite"
params.domains_to_consider = 50 // Number of protein domains to include in the pipeline
params.proteins_per_class = 10 // Number of protein candidates to consider per class of AlphaFold confidence (<50, 50-70, >70 plddt). 10 more will be included and discarded by later script
params.colabfold_database = "/hps/nobackup/agb/research/francesco/data/colabfold/20240226_databases"
//params.colabfold_database =  """${workflow.projectDir}/assets/" //
params.colabfold_cache = "${workflow.projectDir}/tmp/"
params.test_msa = "${workflow.projectDir}/test"

// use module aliasing to call the same process multiple times
include { RANK as RANK_NO_TEMPLATE } from "${workflow.projectDir}/modules/af2rank" params(outdir: params.outdir, projectDir: "${workflow.projectDir}")
include { RANK as RANK_TEMPLATE_MSA } from "${workflow.projectDir}/modules/af2rank" params(outdir: params.outdir, projectDir: "${workflow.projectDir}")
include { RANK as RANK_TEMPLATE_SEQ } from "${workflow.projectDir}/modules/af2rank" params(outdir: params.outdir, projectDir: "${workflow.projectDir}")

include { MOLPROBITY as MOLPROB_NO_TEMPLATE } from "${workflow.projectDir}/modules/molprobity" params(outdir: params.outdir, projectDir: "${workflow.projectDir}")
include { MOLPROBITY as MOLPROB_TEMPLATE_MSA } from "${workflow.projectDir}/modules/molprobity" params(outdir: params.outdir, projectDir: "${workflow.projectDir}")
include { MOLPROBITY as MOLPROB_TEMPLATE_SEQ } from "${workflow.projectDir}/modules/molprobity" params(outdir: params.outdir, projectDir: "${workflow.projectDir}")

process DOWNLOAD_COLABFOLD_WEIGHTS {

    container 'docker://ghcr.io/sokrypton/colabfold:1.5.3-cuda12.2.2'
    time '20m'
    memory '200MB'
    tag 'download_colabfold_weights'
    publishDir "${params.colabfold_cache}/", mode: 'copy', overwrite: true
    errorStrategy 'retry' // tries a second time
        
    output:
    path 'params*', emit: weights
    
    """
    download_colabfold_weights.py
    """
}

process QUERY_DB {
    label 'python38'
    memory '500MB'
    tag 'query-db'
    time '2h'
    publishDir "$params.outdir/DB_query", mode: 'copy'
  
    output:
    path 'proteins_per_domain.csv', emit: proteins_per_domain
    
    beforeScript 'module load mysql'

    """
    query-db_v3.py --database ${params.afdatabase} --number_to_consider ${params.domains_to_consider} --proteins_per_class ${params.proteins_per_class}
    """
}

process QUERY_DB_TEST {
    label 'python38'
    memory '500MB'
    tag 'query-db'
    time '2h'
    publishDir "$params.outdir/DB_query", mode: 'copy'
  
    output:
    path 'proteins_per_domain.csv', emit: proteins_per_domain
    
    """
    echo "PF11361 A0A368WH52_9FIRM,A0A315REV6_9FIRM,A0A1E5HXZ1_9FIRM,A0A1I0AKN1_9FIRM,A0A1N6XNT3_9FIRM,A0A4R6SCZ4_9FIRM,M5ECE6_9FIRM,A0A4R6LDU5_9FIRM,E3DM68_HALPG,A0A2T5RRP8_9FIRM,A0A4R7Z519_9FIRM,E4RK56_HALHG,A0A7X6TSH1_CLOSP,A0A143YLW6_9LACT,A0A1W1IE60_9LACT,A0A847D4P1_9LACT,A0A0R2PT75_9MICO,A0A5C5ECL3_9LACT,A0A3P7S2P6_9FIRM,A0A1H6ETX2_9ACTN" > proteins_per_domain.csv
    echo "PF15176 A0A3B5R9I5_XIPMA,A0A3Q7THB7_VULVU,A0A663F048_AQUCH,K7FNB5_PELSI,A0A811ZG69_NYCPR,A0A3P9IZ15_ORYLA,A0A383ZP48_BALAS,A0A3B3WEP9_9TELE,A0A6P9BDQ4_PANGU,G1NBL5_MELGA,A0A2I0M5B1_COLLI,A0A3Q1M5I3_BOVIN,K4G0F4_CALMI,A0A1U7TKA7_CARSF,A0A6I9IN50_VICPA,A0A5E4CM12_MARMO,A0A5C6NUQ8_9TELE,A0A1S3FIE4_DIPOR,A0A6P5PGU9_MUSCR,A0A5N4EAV5_CAMDR" >> proteins_per_domain.csv
    """
}

process QUERY_DB_WITH_INPUT {
    label 'python38'
    memory '500MB'
    tag 'query-db'
    time '2h'
    publishDir "$params.outdir/DB_query", mode: 'copy'
  
    output:
    path 'proteins_per_domain.csv', emit: proteins_per_domain
    
    beforeScript 'module load mysql'
    """
    query-db_input.py --database ${params.afdatabase} --max_proteins ${params.proteins_per_class} --proteins_table ${params.proteins_table}
    """
}

process GET_SEQUENCES {
    
    memory '500MB'
    tag 'get_sequences'
    time '30m'
    publishDir "$params.outdir/DB_query", mode: 'copy'

    input:
    path proteins_per_domain

    output:
    path 'target_proteins.tsv', emit: proteins_info
    
    beforeScript 'module load mysql'
    
    script:
    """
    get_sequences.sh -i ${proteins_per_domain} -o target_proteins.tsv -n ${params.proteins_per_class} -u ${params.pfamuser} -p ${params.pfampassword} -h ${params.pfamhost} -P ${params.pfamport}
    """
}

process DOWNLOAD_DOMAIN_LENGTHS {

    time '2h'
    tag 'download_domain_lenghts'
    memory '500MB'
    publishDir "${params.outdir}/DB_query/", mode: 'copy'
    
    input:
    path proteins_info
    
    output: 
    path 'domain_lengths.csv'
    
    beforeScript 'module load mysql'
    
    script:
    """
    domain_lenghts.sh -i ${proteins_info} -u ${params.pfamuser} -p ${params.pfampassword} -h ${params.pfamhost} -P ${params.pfamport}
    """

}

process GET_MSAS {
   
    //container 'docker://athbaltzis/colabfold_proteinfold:1.1.0' // this is failing with current db
    //container 'docker://230218818/colabfold:1.5.5' // not having colabfold_search command
    container 'docker://ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2' // has a bug and does not save MSAs under protein name
    //container 'docker://ghcr.io/sokrypton/colabfold:1.5.3-cuda12.2.2' same
    //container 'docker://ghcr.io/sokrypton/colabfold@sha256:ed7c0bb3280379cf1c5831d8229517b36bfb40cffe27e665adee1d0b6e8d8355' same
    //container 'docker://nfcore/proteinfold_colabfold:1.0.0' // no colabfold_search command
    //container 'docker://ghcr.io/sokrypton/colabfold@sha256:9576b632fccdc80f499e79ef5190e86ec7698e30b2ca012748e8eb176201b9bc' same MSA name bug
    //container 'docker://ghcr.io/labdao/colabfold:nightly' no colabfold_search command
    //container 'docker://ghcr.io/sokrypton/colabfold@sha256:40fec143e49896c06761cc02086d44d8303715e07b219d8c1ddc21cb7b41ebef' same MSA name bug
    //container 'docker://lmansouri/colabfold:1.1.0' not compatible with db

    
    memory '180GB' // normally 180
    tag 'get_msas'
    time '48h'
    publishDir "$params.outdir/AF2/MSAs", mode: 'copy'
    
    input:
    path proteins_info
    
    output:
    path ("PF*"), emit: path_to_msa

    """
    #generate_MSAs.sh -i ${proteins_info} -d ${params.colabfold_database}
    # get sequences from file
    tail -n +2 ${proteins_info} | cut -d ' ' -f 1,2,11 | awk -F ' ' '{print ">"\$1"_"\$2"\\n"\$3""}' > sequences.fa
    # generate MSAs
    colabfold_search.py --use-env 0 --use-templates 0 --threads 96 sequences.fa ${params.colabfold_database} ./
    # split MSAs into domain folders
    for MSA in \$( ls *.a3m )
    do
        domain=\$( head \$MSA -n 1 | tr ">" " " | tr "_" " " | awk '{ print \$1 }' )
        protein=\$( head \$MSA -n 1 | tr ">" " " | tr "_" " " | awk '{ print \$2 }' )
        mkdir -p \$domain
        cp \$MSA \$domain/\$protein.a3m
    done
    """

}

process GET_TEST_MSA {
    memory '300MB'
    time '10m'
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
   
    container 'docker://ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2'
    tag 'predict_no_template'
    memory '20G'
    time '10h'
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
    
    container 'docker://ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2'

    tag 'get_templates'
    time '10m'
    memory '200MB'
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

    container 'docker://ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2'

    tag 'predict_no_template'
    memory '20G'
    time '10h'
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

    container 'docker://ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2'

    tag 'predict_no_template'
    memory '20G'
    time '10h'
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
    time '2h'
    tag 'collect_results'
    memory '4GB'
    publishDir "$params.outdir/summary_results/", mode: 'copy' // this finishes after pipeline has completed lol
    // saveAs

    input:
    val rank_mock
    val molprob_mock


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
                        --database_res_dir \$outdir/DB_query/ --molprobity_res_dir \$outdir/molprobity/
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

        MOLPROB_NO_TEMPLATE.out.summary.flatten()\
            .mix(MOLPROB_TEMPLATE_MSA.out.summary.flatten(), MOLPROB_TEMPLATE_SEQ.out.summary.flatten())\
            .collectFile(skip: 1, keepHeader: true)\
            .last()\
            .set{ molprob_out }
        
    emit:
        rank_out = rank_out
        molprob_out = molprob_out

} 

workflow collect_results {
    take:
        rank_out
        molprob_out

    main:
        COLLECT_RESULTS(rank_out, molprob_out)
}

workflow {
    
    // classic workflow to run AF2Fix on a set of PFAM domain proteins
    preparation()
    predictions(preparation.out.weights, preparation.out.proteins_info)
    // maybe this way collect_result() will run after files have been moved to 
    // outputdir by predictions()
    collect_results(predictions.out.rank_out, predictions.out.molprob_out)
}

workflow improve {
    
    // workflow to improve specific proteins using PFAM googd qulaity structures as templates
    preparation_improve()
    predictions(preparation_improve.out.weights, preparation_improve.out.proteins_info)
    collect_results(predictions.out.rank_out, predictions.out.xhpi_out, predictions.out.procheck_out)

}

workflow test {

    // Test the pipeline on a small set
    DOWNLOAD_COLABFOLD_WEIGHTS()
    QUERY_DB_TEST()
    GET_SEQUENCES(QUERY_DB_TEST.out.proteins_per_domain)
    DOWNLOAD_DOMAIN_LENGTHS(GET_SEQUENCES.out.proteins_info)

    GET_MSAS(GET_SEQUENCES.out.proteins_info)
    // 
    PREDICT_NO_TEMPLATE(GET_MSAS.out.path_to_msa.flatten(), DOWNLOAD_COLABFOLD_WEIGHTS.out.weights)

    // create a set of unique paths to domain folders containing no template MSA AlphaFold results
    //PREDICT_NO_TEMPLATE.out.pdb_structure.collect().flatMap{ file -> "${file.toString().split('/')[-2]}" }.unique().set{ no_template_domains }
    GET_TEMPLATES(PREDICT_NO_TEMPLATE.out.output_and_msa)
    //GET_TEMPLATES.out.templates_msa.view { "$it" }
    
    PREDICT_TEMPLATE_MSA(GET_TEMPLATES.out.templates_msa, DOWNLOAD_COLABFOLD_WEIGHTS.out.weights)
    PREDICT_TEMPLATE_SEQ(GET_TEMPLATES.out.templates_msa, DOWNLOAD_COLABFOLD_WEIGHTS.out.weights)
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

    MOLPROB_NO_TEMPLATE.out.summary.flatten()\
        .mix(MOLPROB_TEMPLATE_MSA.out.summary.flatten(), MOLPROB_TEMPLATE_SEQ.out.summary.flatten())\
        .collectFile(skip: 1, keepHeader: true)\
        .last()\
        .set{ molprob_out }

    COLLECT_RESULTS(rank_out, molprob_out)
    
}
