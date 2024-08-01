#!/usr/bin/env bash
#SBATCH -o .log
#SBATCH -e .err
#SBATCH --mem 20GB
#SBATCH --job-name "AF2Fix"

set -euo pipefail
IFS=$'\n\t'

[ -d work/tmp ] || mkdir -p work/tmp

module load openjdk-16.0.2-gcc-9.3.0-xyn6nf5 || echo "Could not load module openjdk-16.0.2-gcc-9.3.0-xyn6nf5"
module load nextflow-22.10.1-gcc-11.2.0-ju5saqw || echo "Could not load module nextflow-22.10.1-gcc-11.2.0-ju5saqw"

export NXF_SINGULARITY_CACHEDIR='assets/'
export NXF_OPTS='-Dnxf.pool.type=sync -Dnxf.pool.maxThreads=10000'
export HADOOP_USER_CLASSPATH_FIRST=true

# in case of building errors
# export SINGULARITY_DISABLE_CACHE=yes

nextflow run -c .config pipeline.nf -with-trace -resume
