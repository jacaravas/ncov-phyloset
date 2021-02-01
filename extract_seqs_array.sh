#! /bin/bash -l

read -a PARAM <<< $(/bin/sed -n ${SGE_TASK_ID}p $1)

#$ -N EXTRACT_SEQS
#$ -cwd
#$ -q all.q
NOW=$( date '+%F_%H:%M:%S' )
echo $NOW
#$ -o logs/$JOB_NAME_$JOB_ID.$TASK_ID.out
#$ -e logs/$JOB_NAME_$JOB_ID.$TASK_ID.err
#$ -pe smp 1

module load Python3/default
module load mafft

module list

# Parameters passed by job control file
reference=${PARAM[0]}
query=${PARAM[1]}
output=${PARAM[2]}

echo "Variables read from job control file"
echo "$1"
echo "$reference"
echo "$query"
echo "$output"

python extract_sequences.py --reference "$reference" --query "$query" --name "$output" --genome

module unload Python3/default
module unload mafft
