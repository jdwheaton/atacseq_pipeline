#!/bin/bash
#SBATCH --mem=2G
#SBATCH -o snakemake.out
#SBATCH -e snakemake.err
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mail-type=END

# This script reads from "config.yml" and executes a specific snakemake pipeline
# depending on the value of the "paired" variable.

# Borrowed this YAML parser from Stefan Farestam on Stackoverflow
function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

eval $(parse_yaml config.yml)

if [ "$paired" == true ]; then
  # Execute sbatch command for paired end
  echo "Executing paired-end workflow"
  snakemake -s atacseq_pe_pipeline.py -p -j 100 --latency-wait 120 \
  --use-singularity --singularity-args "-H $PWD -B /work/mc394/" \
  --cluster-config cluster.json \
  --cluster "sbatch -n {threads} --mem={cluster.mem} -t {cluster.time} \
  -o {cluster.output} -e {cluster.error} -p {cluster.partition} --mail-type=FAIL"
elif [ "$paired" == false ]; then
  # Execute sbatch command for single end
  echo "Executing single-end workflow"
  snakemake -s atacseq_se_pipeline.py -p -j 100 --latency-wait 120 \
  --use-singularity --singularity-args "-H $PWD -B /work/mc394/" \
  --cluster-config cluster.json \
  --cluster "sbatch -n {threads} --mem={cluster.mem} -t {cluster.time} \
  -o {cluster.output} -e {cluster.error} -p {cluster.partition} --mail-type=FAIL"
else
  echo "You did not specify a value for 'paired' in the config file!"
fi
