#!/bin/bash
# Pathogen-Host Protein Interaction Analysis Pipeline
# ----------------------------------------------------------
# This script processes pathogen and host protein sequences to identify interactions.
# Tools used: SignalP6.0, OrthoFinder, CD-HIT, AlphaFold3
# Arguments that have to be specified by users:
# -a: Account name for Slurm
# -e: Email for notifications
# -p: Pathogen FASTA directory
# -i: Host FASTA directory
# -l: The directory for AlphaFold3 model parameters, which can be obtained from https://docs.google.com/forms/d/e/1FAIpQLSfWZAgo1aYk0O4MuAXZj8xRQ8DafeFJnldNOnh_13qAx2ceZw/viewform
# Output:
# Results are saved in specified output directory, default './AF3_out'.


#-------------------------------------------------
# Receive arguments from user inputs
#-------------------------------------------------
# Display usage instructions
function usage() {
    echo "Usage: $0 -a <account> 
    -e <email> 
    -c <cpu, default=4> 
    -m <memory in GB, default=62> 
    -d <days, default=1> 
    -p <pathogen fasta dir> 
    -i <host fasta dir> 
    -l <model param dir> 
    -f <cd-hit cutoff, range from 0.4 to 1, default=0.5> 
    -o <output dir, default './AF3_out'>"
    exit 1
}

# Initialize default values
cpu=4
mem=62
days=1
cd_c=0.5
param_dir='/blue/jhuguet/z.bao/AF3_model_parameters/af3.bin'
AFOUT=$(pwd)/AF3_out

# Parse command-line arguments using getopts
while getopts ":a:e:c:m:d:p:i:l:f:o:h" opt; do
    case $opt in
        a) acnt="$OPTARG" ;;            # Account name
        e) email="$OPTARG" ;;           # Email address
        c) cpu="$OPTARG" ;;             # Number of CPUs
        m) mem="$OPTARG" ;;             # Memory (in GB)
        d) days="$OPTARG" ;;            # Time (in days)
        p) input_fa_dir_p="$OPTARG" ;;  # Pathogen fasta directory
        i) input_fa_dir_h="$OPTARG" ;;  # Host fasta directory
        l) param_dir="$OPTARG" ;;       # Model parameter directory
        f) cd_c="$OPTARG" ;;            # CD-HIT cutoff value
        o) AFOUT="$OPTARG" ;;           # Output directory
        h) usage ;;                     # Show usage instructions
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# Ensure mandatory arguments are provided
if [[ -z "$acnt" || -z "$email" || -z "$input_fa_dir_p" || -z "$input_fa_dir_h" || -z "$param_dir" ]]; then
    echo "Missing required arguments." >&2
    usage
fi

# Validate memory value
if ! [[ "$mem" =~ ^[0-9]+$ ]]; then
    echo "Error: Memory value (-m) must be a positive integer." >&2
    exit 1
fi

# Validate CD-HIT cutoff
if (( $(echo "$cd_c < 0.4" | bc -l) || $(echo "$cd_c > 1" | bc -l) )); then
    echo "Error: CD-HIT cutoff (-f) must be between 0.4 and 1." >&2
    exit 1
fi

# Ensure necessary directories are absolute paths
input_fa_dir_p=$(realpath "$input_fa_dir_p")
input_fa_dir_h=$(realpath "$input_fa_dir_h")
param_dir=$(realpath "$param_dir")
AFOUT=$(realpath "$AFOUT")

# Check if required directories exist and are writable
function validate_directory() {
    local dir="$1"
    if [[ ! -d "$dir" || ! -r "$dir" ]]; then
        echo "Error: Directory $dir does not exist or is not readable." >&2
        exit 1
    fi
}
validate_directory "$input_fa_dir_p"
validate_directory "$input_fa_dir_h"
validate_directory "$param_dir"

# set a date time for all steps whenever possible
datetime=$(date +%Y%m%d_%H%M%S) # One time stamp for all steps

#-------------------------------------------------
# Generate a Slurm batch script dynamically
#-------------------------------------------------
cat <<EOF > slurm_job_${datetime}.sh # accept parameters from previous section
#!/bin/sh
#SBATCH --account=$acnt
#SBATCH --qos=$acnt
#SBATCH --partition=gpu
#SBATCH --gpus=a100:1
#SBATCH --job-name=pathogen-host_protein_interactions
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=$email
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$cpu
#SBATCH --mem=${mem}GB
#SBATCH --time=${days}-00:00:00
#SBATCH --output=log/pathogen-host_protein_interactions_%j.out

hostname;pwd
date
mkdir -p log

# parameters from user input arguments
input_fa_dir_p=$input_fa_dir_p
input_fa_dir_h=$input_fa_dir_h
cd_c=$cd_c
param_dir=$param_dir
AFOUT=$AFOUT

datetime=$datetime

EOF

cat <<'EOF' >> slurm_job_${datetime}.sh # quote EOF to avoid shell interpretation. '>>' Append.


# STEP 1 PATHOGEN ANALYSIS
#1) SignalP prediction
#2) OrthoFinder processing
#3) CD-HIT clustering
##----------------------------------------------------

echo "STEP 1 PATHOGEN ANALYSIS"
# step 1.1 signalp to predict secret protein
echo "step 1.1 pathogen signalP6.0 "
ml signalp/6.0h
ml seqtk

for fa_file in $(ls $input_fa_dir_p/*.fa* | awk -F '/' '{print $NF}');  # -F'/' field separator $NF refers to the last field in the record. # in case of .faa, fasta
do
    echo $fa_file
    out_sp6_dir_p=$(pwd)/temp/p1_SignalP6_out/sp${datetime}
    mkdir -p $out_sp6_dir_p/$fa_file
    mkdir -p $out_sp6_dir_p/trimedSPheaders_forOrtho
    
    # 1.1.1 signalP6
    signalp6 --fastafile $input_fa_dir_p/$fa_file --organism other --output_dir /$out_sp6_dir_p/$fa_file # using abs path for --output_dir is a must!
    
    # 1.1.2 process signalP6 results for OrthoFinder
    # get seq names that are marked with "signal_peptide"
    for seq_nm in $(grep -w "signal_peptide" $out_sp6_dir_p/$fa_file/output.gff3 | awk '{print $1}') # cut with default delimter may not be good here.
    do
        # store the single seq name for seqtk
        echo $seq_nm > $out_sp6_dir_p/temp_protein_nm.txt
        
        # extract the single protein
        seqtk subseq $input_fa_dir_p/$fa_file $out_sp6_dir_p/temp_protein_nm.txt > $out_sp6_dir_p/temp_protein.fa
        
        # get the end positon
        end_pos=$(grep $seq_nm $out_sp6_dir_p/$fa_file/output.gff3 | cut -f5 ) # awk is not good here since the name column size varies
        
        # trim first end_pos protein bases for each single protein
        seqtk trimfq -b $end_pos $out_sp6_dir_p/temp_protein.fa >> $out_sp6_dir_p/trimedSPheaders_forOrtho/$fa_file'_SP_Trimed_SPHeaders_ForOrtho.fasta'
    done
done


# Step 1.2 OrthoFinder
ml orthofinder

# Check the number of fa files in $input_fa_dir_p
num_files_p=$(find "$input_fa_dir_p" -type f -name "*.fa*" | wc -l)

if (( num_files_p > 1 )); then
    echo "Number of pathogen FASTA files: $num_files_p"
    echo "step 1.2 Running OrthoFinder for pathogen files..."
    
    # Run OrthoFinder for Pathogen
    out_of_dir_p=$(pwd)/temp/p2_OrthoF_out
    mkdir -p $out_of_dir_p
    orthofinder -f $out_sp6_dir_p/trimedSPheaders_forOrtho/ -a 6 -o $out_of_dir_p/of${datetime} # DO NOT mkdir the output folder, orelse ERROR: directory already exists
    
    # specify input fasta directory for CD-HIT
    cd_input_dir=$(ls -td $out_of_dir_p/of${datetime}/*/ | head -1)Orthogroup_Sequences # * example Results_Apr25. $() returned with '/'
else
    echo "Number of pathogen FASTA files is $num_files_p. Skipping step 1.2 OrthoFinder for pathogen."
    cd_input_dir=$input_fa_dir_p
fi


# step 1.3 CD-HIT
echo "step 1.3 pathogen cdhit"

ml cdhit
out_cd_dir_p=$(pwd)/temp/p3_cdhit_out/"cdhit_"$datetime
mkdir -p $out_cd_dir_p

# 1.3.1 determine n based on c. -n 5 for thresholds 0.7~1.0 {4:0.6~0.7, 3:0.5~0.6, 2:0.4~0.5}
case $cd_c in
    0.[7-9]*|1.[0]*) cd_n=5 ;;   # thresholds 0.7~1.0
    0.[6]*) cd_n=4 ;;            # thresholds 0.6~0.7
    0.[5]*) cd_n=3 ;;            # thresholds 0.5~0.6
    0.[4]*) cd_n=2 ;;            # thresholds 0.4~0.5
    *) echo "Invalid cutoff value"; exit 1 ;;
esac

# 1.3.2 run CD-HIT
for fa in $( ls $cd_input_dir | awk -F '/' '{print $NF}')
do
    cd-hit -i $cd_input_dir/$fa -o $out_cd_dir_p/'representative_protein_from_'$fa'.fa' -c $cd_c -n $cd_n #
done

# 1.3.3 merge fas to one fa
cat $out_cd_dir_p/*.fa | awk 'NF' > $(pwd)/temp/p3_cdhit_out/$datetime'_final_pathogen_proteins.fa' # awk 'NF' remove empty lines

# 1.3.4 provide final path for pathogen proteins
fa1_path=$(pwd)/temp/p3_cdhit_out/$datetime'_final_pathogen_proteins.fa'




# STEP 2 HOST ANALYSIS
#1) OrthoFinder processing
#2) CD-HIT clustering
##----------------------------------------------------

echo "STEP 2 HOST ANALYSIS"
# Step 2.1 OrthoFinder
# Check the number of files in $input_fa_dir_h
num_files_h=$(find "$input_fa_dir_h" -type f -name "*.fa*" | wc -l)

if (( num_files_h > 1 )); then
    echo "Number of host FASTA files: $num_files_h"
    echo "step 2.1 Running OrthoFinder for host files..."
    
    # Run OrthoFinder for Host
    out_of_dir_h=$(pwd)/temp/h2_OrthoF_out
    mkdir -p "$out_of_dir_h"
    orthofinder -f $input_fa_dir_h -a 6 -o $out_of_dir_h/of${datetime} # DO NOT mkdir the output folder, orelse ERROR: directory already exists
    
    # specify input fasta directory for CD-HIT
    cd_input_dir=$(ls -td $out_of_dir_h/of${datetime}/*/ | head -1)Orthogroup_Sequences # * example Results_Apr25. $() returned with '/'
else
    echo "Number of host FASTA files is $num_files_h. Skipping step 2.1 OrthoFinder for host."
    cd_input_dir=$input_fa_dir_h
fi


# step 2.2 CD-HIT
echo "step 2.2 Host CD-HIT"
out_cd_dir_h=$(pwd)/temp/h3_cdhit_out/"cdhit_"$datetime
mkdir -p $out_cd_dir_h

# 2.2.1 run CD-HIT
for fa in $( ls $cd_input_dir | awk -F '/' '{print $NF}')
do
    cd-hit -i $cd_input_dir/$fa -o $out_cd_dir_h/'representative_protein_from_'$fa'.fa' -c $cd_c -n $cd_n # -n 5 for thresholds 0.7~1.0 {4:0.6~0.7, 3:0.5~0.6, 2:0.4~0.5}
done

# 2.2.2 merge fas (a fa has >=1 seqs) to one fa
cat $out_cd_dir_h/*.fa | awk 'NF' > $(pwd)/temp/h3_cdhit_out/$datetime'_final_host_proteins.fa' # awk 'NF' remove empty lines

# 2.2.3 provide final path for host proteins
fa2_path=$(pwd)/temp/h3_cdhit_out/$datetime'_final_host_proteins.fa'




# STEP 3 ALPHAFOLD3
#1) prepare josn files
#2) run AlphaFold3
#3) summarize AF3 results
##----------------------------------------------------

echo "STEP 3 ALPHAFOLD3"
# step 3.1 prepare json files from fa files
echo "step 3.1 prepare json files from fa files"
js_out_dir=$(pwd)/temp/AF3_input
mkdir -p $js_out_dir

python prepare_json_from_fa.py --fa1_path=$fa1_path --fa2_path=$fa2_path --out_dir=$js_out_dir --num 100 --today=$datetime  # no "," between arguments
echo "prepare json done!"


# step 3.2 run AF3 on HPG
echo "step 3.2 run AF3 on HPG"
ml singularity
# js_dir=$(ls -td $js_out_dir/*/ | head -1) # another layer of folder within js_out_dir, get the latest one. "/" is here.
js_dir=$(ls -d $js_out_dir/${datetime}*/) # another layer of folder within js_out_dir, get the one with dtatetime. "/" is here.

# The AF3 results from one json will be put into one res_folder
for stem in $(ls $js_dir | cut -d "_" -f3- | cut -d "." -f1 ); # if too many jsons at a time, add more pipes such as | head -3 or head -6 | tail -3
do
    echo $stem
    json=$(ls $js_dir*$stem*) # input json file
    echo $json
    
    afout=$AFOUT/$stem # output AF3 directory. Results from one json will be put into one directory
    mkdir -p $afout
    
    singularity exec \
         --nv \
         --bind $js_dir:/root/af_input \
         --bind $afout:/root/af_output \
         --bind $param_dir:/root/models \
         --bind /data/reference/alphafold/3.0.0:/root/public_databases \
         /apps/alphafold/3.0.0/alphafold3.sif \
         python /app/alphafold/run_alphafold.py \
         --db_dir=/data/reference/alphafold/3.0.0 \
         --json_path=$json \
         --model_dir=$param_dir \
         --output_dir=$afout
done

# step 3.3 post analysis the AF3 results
echo "step 3.3 post analysis the AF3 results"
ml yolo/v8 # to provide pandas

python post_analysis.py --af3_out_dir=$AFOUT --summary_path=$AFOUT/af3_results_summary_${datetime}.csv

EOF


#-------------------------------------------------
# Submit the dynamically created job
#-------------------------------------------------
job_id=$(sbatch slurm_job_${datetime}.sh | awk '{print $4}')
if [[ -z "$job_id" ]]; then
    echo "Error: Failed to submit Slurm job." >&2
    exit 1
fi
echo "Submitted batch job ID: $job_id"
