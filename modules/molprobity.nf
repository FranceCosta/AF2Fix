process MOLPROBITY {
    container "docker://francecosta/molprobity:v0.0.1"
    tag 'molprobity'
    memory '200MB'
    publishDir "$params.outdir/molprobity/${condition}/${domain}", mode: 'copy'
    //cache false
    errorStrategy 'ignore'

    input:
    tuple path(structure), val(domain)
    val condition

    output:
    path '*.csv', emit: summary

    """
    for structure in \$( ls *0.pdb ); do
        cp \$structure \$( echo \$structure | cut -d _ -f 1 ).pdb
        structure=\$( echo \$structure | cut -d _ -f 1 ).pdb
        molprobity.py --PDB_file \$structure
    done
    """
}