#!/bin/bash
# bash virtual_pcr.sh -f CCCAGTACAGCCCATTCACC -r AAAGCGGATAGGGCTCCTGTA -i FAT05083_pass_a4a97209_0.fastq.gz -m 7500 -l 9000


# ============================================================================================================
# The first part of the program extracts sequences between primers. 
# The implementation was extracted from this project in github https://github.com/Nucleomics-VIB/InSilico_PCR

version="1.1.3; 2022-02-08"

source /home/bioinfo/miniconda3/etc/profile.d/conda.sh
conda activate InSilico_PCR || \
  ( echo "# the conda environment 'InSilico_PCR' was not found on this machine" ;
    echo "# please read the top part of the script!" \
    && exit 1 )

########## Variables Initialization ##############

# REM: the workflow can be restarted by deleting all log files after a given step

# speed-up
thr=8
pigt=2
mem="1g"

# extract reads corresponding to the 16S PCR
forwardl="Forward"
reversel="Reverse"

# be stringent to avoid noisy reads
cut=1.0

barcodes_len=5

# split in 500k read bins and zip
lines=2000000

# quality phred scale of your data (for demo data it is 33)
qual=33

# expected amplicon size limits, set to 10 and 10000 by default
# adjust if you notice that the sequence extraction has unwanted tail(s)
readminlen=10
readmaxlen=100000



# whether to include the primer matches or to clip them (t/f)
primincl="t"

######## end of user editable region ############

# Arguments

Help()
{
   # Display Help
   echo "Optional Arguments:"
   echo "   -h      Print this help message"
   echo "   -t      Number of threads. [DEFAULT: 8]"
   echo "   -q      Quality phred scale of your data. [DEFAULT: 33]"
   echo "   -m      Expected amplicon size limits, minimum [DEFAULT: 10]"
   echo "   -l      Expected amplicon size limits, limit [DEFAULT: 100000]"
   echo "   -c      Cutoff value (identity/similarity) [0-1] [DEFAULT: 0.8]"
   echo "   -b      Barcodes length. [DEFAULT: 5]"
   echo 
   echo "Required Arguments:"
   echo "   -n      Primer Pair Name"
   echo "   -f      Forward Primer"
   echo "   -r      Reverse Primer"
   echo "   -i      Path to input file containing the sequences"
}


while getopts hc:b:t:q:m:M:n:f:r:i: flag
do
    case "${flag}" in
        h) # display Help
            Help
            exit 1;;
        c) cut=${OPTARG};;
        b) barcodes_len=${OPTARG};;
        t) thr=${OPTARG};;
        q) qual=${OPTARG};;
        m) readminlen=${OPTARG};;
        M) readmaxlen=${OPTARG};;
        n) primername=${OPTARG};;
        f) forwardp=${OPTARG};;
        r) reversep=${OPTARG};;
        i) infile=${OPTARG};;
    esac
done

shift $(( OPTIND - 1 ))

if [ -z "$primername" ]; then
        echo 'Missing -n (Primer Pair Name)'
        exit 1
elif [ -z "$forwardp" ]; then
        echo 'Missing -f (Forward Primer)'
        exit 1
elif [ -z "$reversep" ]; then
        echo 'Missing -r (Reverse Primer)'
        exit 1
elif [ -z "$infile" ]; then
        echo 'Missing -i (Sequences input file)'
        exit 1      

fi

# extract string for output names
name=$(basename ${infile%\.fq\.gz})

##################
# prepare folders

# working default to local folder
data="$(pwd)"


# split folder
split="split_data_${name}"
mkdir -p ${split}

# run logs
logs="run_logs_${name}"
mkdir -p ${logs}

# tmp output folder
tmpout="bbmap_out"
mkdir -p ${tmpout}

# parallel folders
WORKDIR=$PWD
mkdir -p $PWD/tmp
export TMPDIR=$PWD/tmp

# keep track of all
runlog=${logs}/runlog.txt
exec &> >(tee -i ${runlog})


###########################################
# split large file into chunks for parallel

if [[ ! -f ${logs}/done.splitting ]]; then
  echo "# splitting the data in multiple smaller files and compressing (may take some time!)"
  zcat ${infile} | \
    split -a 3 -d -l ${lines} --filter='pigz -p '${pigt}' \
    > $FILE.fq.gz' - ${split}/${name}_ && \
  touch  ${logs}/done.splitting
else
  echo "# splitting already done"
fi

#################################
# run in BBMap msa.sh in parallel
#################################

# compute job number & threads/job
jobs=$(ls ${split}|wc -l)
# limit to thr
if [[ "${jobs}" -gt "${thr}" ]]; then
	jobs=${thr}
fi
jobt=$((${thr}/${jobs}))

######################################################
# find forward primer using a fraction of the threads


if [[ ! -f ${logs}/done.searching.${name}_${forwardl} ]]; then
  echo "# searching for forward primer sequence: ${forwardp} in all files"
  find ${split} -type f -name "${name}_???.fq.gz" -printf '%P\n' |\
    sort -n |\
    parallel --workdir ${WORKDIR} --tmpdir ${TMPDIR} -j ${jobs} msa.sh -Xmx${mem} threads=${jobt} \
      qin=${qual} \
      in=${split}/{} \
      out=${tmpout}/forward.sam \
      literal="${forwardp}" \
      rcomp=t \
      addr=t \
      replicate=t \
      cutoff="${cut}" && \
    touch ${logs}/done.searching.${name}_${forwardl}
else
  echo "# forward search already done"
fi



######################################################
# find reverse primer using a fraction of the threads

if [[ ! -f ${logs}/done.searching.${name}_${reversel} ]]; then
  echo "# searching for reverse primer sequence: ${reversep} in all files"
  find ${split} -type f -name "${name}_???.fq.gz" -printf '%P\n' |\
    sort -n |\
    parallel --workdir ${WORKDIR} --tmpdir ${TMPDIR} -j ${jobs} msa.sh -Xmx${mem} threads=${jobt} \
      qin=${qual} \
      in=${split}/{} \
      out=${tmpout}/reverse.sam \
      literal="${reversep}" \
      rcomp=t \
      addr=t \
      replicate=t \
      cutoff="${cut}" && \
  touch ${logs}/done.searching.${name}_${reversel}
else
  echo "# reverse search already done"
fi


python stretch_length_for_barcodes.py ${tmpout}/forward.sam ${tmpout}/reverse.sam ${barcodes_len}

##########################################
# extract regions with BBMap cutprimers.sh

if [[ ! -f ${logs}/done.cutprimer.${name}_${forwardl}_${reversel} ]]; then
  echo "# extracting template sequences between primer matches"
  find ${split} -type f -name "${name}_???.fq.gz" -printf '%P\n' |\
    sort -n |\
    parallel --workdir ${WORKDIR} --tmpdir ${TMPDIR} -j ${thr} cutprimers.sh -Xmx${mem} \
      qin=${qual} \
      in=${split}/{} \
      out=${tmpout}/{}_16s.fq \
      sam1=${tmpout}/forward.sam \
      sam2=${tmpout}/reverse.sam \
      include=${primincl} \
      fixjunk=t && \
  touch ${logs}/done.cutprimer.${name}_${forwardl}_${reversel}
fi

###########################################################
# merge results and keep only reads in amplicon size range

if [[ ! -f ${logs}/done.merging.${name}_${forwardl}_${reversel} ]]; then
  # clear existing results
  # final="output/${name}_filtered_${forwardl}_${reversel}.fa"
  file_name="${name/.fastq.gz/""}"  
  # final="output/${file_name}_${primername}_extracted_regions.fa"
  final="output/${primername}.fq"
  final2="output/${primername}.fa"

  # cat /dev/null > ${final}
  echo "# filtering results at min:${readminlen} and max:${readmaxlen} and merging to ${final}"
  find ${tmpout} -type f -name "${name}_???.fq.gz_16s.fq" | \
    sort -n |\
    xargs cat |\
    seqkit seq -m "${readminlen}" -M "${readmaxlen}" > $final && \
    seqkit fq2fa $final -o $final2 && \
  touch ${logs}/done.merging.${name}_${forwardl}_${reversel}
else
  echo "# merging and filtering already done"
  echo "# force redo by deleting ${logs}/done.merging.${name}_${forwardl}_${reversel}"
fi


# return to normal
conda deactivate

# ===============================================================================


python sep_samples.py barcodes/barcodes.tsv ${primername}

rm -r $split
rm -r $logs
rm -r $tmpout
rm -r $runlog
rm -r $PWD/tmp

