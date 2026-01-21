#!/bin/bash
#SBATCH -J sarek
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-24:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=5G                         # Memory total in MiB (for all cores)
#SBATCH -o %x_%j.out                 # File to which STDOUT will be written, including job ID (%j)

set -eu
module load java

set -x

ts=`date +%Y%m%d-%H%M%S`

export TMPDIR=/n/scratch/users/c/ch305/tmp
export _JAVA_OPTIONS="-Djava.io.tmpdir=$TMPDIR -Dtmp.dir=$TMPDIR"
export SINGULARITY_WORKDIR=$TMPDIR

nextflow run nf-core/sarek -r 3.5.1 -profile singularity \
  -params-file nf-params.json -c ~/o2.config \
  -with-trace -with-report report-${ts}.html -with-timeline timeline-${ts}.html \
  "$@"

# whole exome sequencing

