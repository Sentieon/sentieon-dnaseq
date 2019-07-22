#!/bin/sh
export TMPDIR=$PWD/tmp
mkdir -p $TMPDIR
export SENTIEON_TMPDIR=$TMPDIR
release=$HOME/ssd/release/sentieon-genomics-201808.05

# Reference files
hg38_dir=$HOME/ssd/reference
fasta=$hg38_dir/Homo_sapiens_assembly38.fasta
known_Mills=$hg38_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
known_1000GIndels=$hg38_dir/1000G_phase1.snps.high_confidence.hg38.vcf.gz
known_dbsnp=$hg38_dir/dbsnp_146.hg38.vcf.gz
calling_intervals_list=$hg38_dir/wgs_calling_regions.hg38.interval_list
# Fastq files
fastq_dir=$HOME/ssd/fastq
FASTQ1=${fastq_dir}/140127_D00360_0011_AHGV6ADXX_R1.fastq.gz
FASTQ2=${fastq_dir}/140127_D00360_0011_AHGV6ADXX_R2.fastq.gz
# Log file
logfile=$PWD/Sentieon_match_4.1_$(date +"%Y%m%d%H%M").log
exec >$logfile 2>&1

# Sample ID
sample_id=NA12878
lib=LIBRARY
ID_SM=${sample_id}
PU=${lib}_${sample_id}
LB=${sample_id}
SM=${sample_id}
# Sequencing platform
PL=ILLUMINA #platform
GROUP_ID="@RG\tID:$ID_SM\tSM:$SM\tPU:$PU\tPL:$PL\tLB:$LB"
prefix=./${sample_id}_
cpu_cores=$(nproc)

#######################################################
check_error() 
{
    if [ "$1" = "0" ]
    then
        echo "$2 completed successfully."
        return
    fi
    echo Error peforming $2.
    exit $1
}

delta_time()
{
    start=$1
    end=$2
    dt=$((end-start))
    dh=$(echo "$dt/3600" | bc)
    dt2=$(echo "$dt-3600*$dh" | bc)
    dm=$(echo "$dt2/60" | bc)
    ds=$(echo "$dt2-60*$dm" | bc)
    echo "$dh:$dm:$ds"
}

run()
{
    cmd=$1
    what=$2
    echo "Starting to run $what: $cmd"
    start=`date +"%D %T"`
    start_s=`date +%s`
    echo "$what start time: $start"
    eval "$cmd"
    check_error $? "$what"
    end=`date +"%D %T"`
    end_s=`date +%s`
    echo "$what end time: $end"
    runtime=$(delta_time $start_s $end_s)
    echo "$what runtime: $runtime"
}

echo "############"
echo "# Sentieon #"
echo "############"

mkdir -p sentieon_match_4.1
cd sentieon_match_4.1
output_prefix=${prefix}sentieon_match_4.1_
# speed up memory allocation malloc in bwa
export LD_PRELOAD=$release/lib/libjemalloc.so
export MALLOC_CONF=lg_dirty_mult:-1
mem_kb=$(cat /proc/meminfo | grep "MemTotal" | awk '{print $2}')
export bwt_max_mem="$((mem_kb / 1024 / 1024 - 2))g"

# ******************************************
#  Pipeline settings and function definitions
# ******************************************
dedup_param="--traverse_param 1000000/10000"
cram_option="--cram_write_options version=3.0,compression_level=1"
bam_option="--bam_compression 1"
ext="bam"

# ******************************************
# 1. Mapping reads with BWA-MEM. Output coordinate-sorted BAM file.
# ******************************************
bwa_opts="-K 10000000 -M -Y"
stage="Alignment"
cmd="$release/bin/sentieon bwa mem -R '$GROUP_ID' $bwa_opts -t $cpu_cores $fasta $FASTQ1 $FASTQ2"
command="($cmd || echo -n 'BWA Error' ) | $release/bin/sentieon util sort $cram_option $bam_option -i - -r $fasta -t $cpu_cores -o ${output_prefix}sorted.$ext --sam2bam"
run "$command" "$stage"

# ******************************************
# 2. Remove Duplicate Reads
# ******************************************
command="$release/bin/sentieon driver -r $fasta -t $cpu_cores -i ${output_prefix}sorted.$ext $dedup_param --algo LocusCollector --fun score_info ${output_prefix}score.txt.gz"
run "$command" "Dedup_prep"
command="$release/bin/sentieon driver -r $fasta -t $cpu_cores -i ${output_prefix}sorted.$ext $dedup_param --algo Dedup --rmdup $cram_option $bam_option --score_info ${output_prefix}score.txt.gz --metrics ${output_prefix}dedup_metrics.txt ${output_prefix}deduped.$ext" 
run "$command" "Dedup"

# ******************************************
# 3. BQSR
# ******************************************

command="$release/bin/sentieon driver -r $fasta -t $cpu_cores -i ${output_prefix}deduped.$ext --algo QualCal -k $known_dbsnp -k $known_1000GIndels -k $known_Mills ${output_prefix}recal_data.table"
run "$command" "BQSR"

# ******************************************
# 4. HC Variant caller
# ******************************************
command="$release/bin/sentieon driver -r $fasta -t $cpu_cores -i ${output_prefix}deduped.$ext -q ${output_prefix}recal_data.table --interval $calling_intervals_list --algo Haplotyper --emit_mode gvcf ${output_prefix}output.g.vcf.gz"
run "$command" "Haplotype GVCF caller"

command="$release/bin/sentieon driver -r $fasta -t $cpu_cores --interval $calling_intervals_list --algo GVCFtyper --genotype_model multinomial -d $known_dbsnp -v ${output_prefix}output.g.vcf.gz ${output_prefix}output.vcf.gz"
run "$command" "GVCFtyper"
