process PROCHECK {
    container "docker://francecosta/ccp4:v0.0.5"
    tag 'procheck'
    queue 'short'
    memory '1GB'
    publishDir "$params.outdir/ramachandran/${condition}/${domain}", mode: 'copy'
    //cache false
    errorStrategy 'ignore'

    input:
    tuple path(structure), val(domain)
    val condition

    output:
    path '*.csv', emit: summary

    //  mv ${condition}/*.csv ${condition}/${domain}/${structure_path.toString().split("/")[-1].split("_relaxed_rank_001")[0]}.csv
    """
    for structure in \$( ls *0.pdb ); do
        cp \$structure \$( echo \$structure | cut -d _ -f 1 ).pdb
        structure=\$( echo \$structure | cut -d _ -f 1 ).pdb
        procheck.py --PDB_file \$structure --output ./
    done
    """

}
