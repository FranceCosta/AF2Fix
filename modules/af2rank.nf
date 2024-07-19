process RANK {
    //conda '/nfs/research/agb/research/francesco/anaconda3/envs/alphafold_v2' // this works
    //conda "$params.projectDir/assets/alphafold_environment_old.yml" // not working, module jax not found
    container "docker://neoformit/alphafold:v2.3.1_2" // this with export XLA_PYTHON_CLIENT_PREALLOCATE=false
    // does not use GPU. Same wo export XLA_PYTHON_CLIENT_PREALLOCATE=false

    module 'cuda-11.1.1-gcc-9.3.0-oqr2b7d'
    tag 'rank'
    memory '20G'
    clusterOptions '--gres=gpu:a100:1'
    publishDir "$params.outdir/AF2Rank/$condition/$domain", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple path(structure), val(domain)
    val condition

    output:
    path '*.pdb', emit: pdb_structure
    path '*.csv', emit: csv_table

    // AF2Rank.py --input_pdb ${structure} --mask_sidechains_add_cb --output_dir ./ ${structure.toString().split("/")[-1].split("_relaxed_rank_001")[0]}
    """
    export XLA_PYTHON_CLIENT_PREALLOCATE=false
    unset TF_FORCE_UNIFIED_MEMORY
    for structure in \$( ls *.pdb ); do
        AF2Rank.py --input_pdb \$structure --mask_sidechains_add_cb --output_dir ./ \$( echo \$structure | cut -d _ -f 1 )
    done
    """

}
