#!/bin/sh
export TMPDIR=$PWD/tmp
mkdir -p $TMPDIR
gatk_Launcher=$HOME/ssd/release/gatk-4.0.2.1/gatk
samtools_release=$HOME/ssd/release/samtools-1.9/bin
bwa_release=$HOME/ssd/release/bwa.kit

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
logfile=$PWD/GATK4.0_vs_Sentieon_$(date +"%Y%m%d%H%M").log
exec >$logfile 2>&1

# Sample ID
sample_id=NA12878
lib=LIBRARY
ID_SM=${sample_id}
PU=${lib}_${sample_id}
LB=${sample_id}
SM=${sample_id}
# Sequencing platform.
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
#######################################################
echo "###########"
echo "# GATK4.0 #"
echo "###########"

mkdir -p gatk4.0
cd gatk4.0

output_prefix=${prefix}gatk4.0_
# Split the whole genome for BQSR scatter-gather
python - $hg38_dir/Homo_sapiens_assembly38.dict $TMPDIR << EOF 
import sys
with open(sys.argv[1], "r") as ref_dict_file:
    sequence_tuple_list = []
    longest_sequence = 0
    for line in ref_dict_file:
        if line.startswith("@SQ"):
            line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
            sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
    longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
# We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
# the last element after a :, so we add this as a sacrificial element.
hg38_protection_tag = ":1+"
# initialize the tsv string with the first sequence
tsv_string = "-L " + sequence_tuple_list[0][0] + hg38_protection_tag
temp_size = sequence_tuple_list[0][1]
for sequence_tuple in sequence_tuple_list[1:]:
    if temp_size + sequence_tuple[1] <= longest_sequence:
        temp_size += sequence_tuple[1]
        tsv_string += " " + "-L " + sequence_tuple[0] + hg38_protection_tag
    else:
        tsv_string += "\n" + "-L " + sequence_tuple[0] + hg38_protection_tag
        temp_size = sequence_tuple[1]
# add the unmapped sequences as a separate line to ensure that they are recalibrated as well
with open("%s/sequence_grouping.txt"%sys.argv[2],"w") as tsv_file:
  tsv_file.write(tsv_string)
  tsv_file.write("\n")
  tsv_file.close()
tsv_string += '\n' + "-L unmapped\n"
with open("%s/sequence_grouping_with_unmapped.txt"%sys.argv[2],"w") as tsv_file_with_unmapped:
  tsv_file_with_unmapped.write(tsv_string)
  tsv_file_with_unmapped.close()
  
EOF

# Split wgs_calling_regions list into 9 sub-regions for HaplotypeCaller scatter-gather
awk 'BEGIN {
  prev_total = 0
  frag = 1
  container = ""
} 
{ if ( $1 !~ /^@/ ) 
  {
    len = ($3 - $2 + 1)
    if ( prev_total + len < 324860607 ) {
      prev_total += len
      container = container sprintf("-L %s:%d-%d ", $1, $2, $3)
    }
    else {
      a1 = prev_total + len - 324860607
      a2 = 324860607 - prev_total
      if ( a1 > a2 ) { print container; container = sprintf("-L %s:%d-%d ", $1, $2, $3); prev_total = len}
      else { container = container sprintf("-L %s:%d-%d ", $1, $2, $3); print container; container = ""; prev_total = 0}
      frag += 1
    }
  }
}
END {
  if ( container ) { print container }
}' $calling_intervals_list > $TMPDIR/ScatterIntervalList.txt

stage="Alignment"
start_time_mod=$(date +%s)
    # Alignment and sort
    bwa_opts="-K 10000000 -M -Y"
    cmd1="$bwa_release/bwa mem $bwa_opts -t $cpu_cores -R '$GROUP_ID' $fasta $FASTQ1 $FASTQ2 \
        | $samtools_release/samtools sort -@ $cpu_cores -o ${output_prefix}sorted.bam -"
    cmd2="$samtools_release/samtools index ${output_prefix}sorted.bam"
    run "$cmd1 && $cmd2" "bwa/sort"
end_time_mod=$(date +%s)
if [[ "$OSTYPE" == "darwin"* ]]; then start_date=$(date -j -f "%s" $start_time_mod); else start_date=$(date -d @$start_time_mod); fi
if [[ "$OSTYPE" == "darwin"* ]]; then end_date=$(date -j -f "%s" $end_time_mod); else end_date=$(date -d @$end_time_mod); fi
echo "Stage "$stage" Started: "$start_date"; Ended: "$end_date"; Elapsed time: "$(($end_time_mod - $start_time_mod))" sec"

stage="Dedup"
start_time_mod=$(date +%s)
    # Remove PCR duplicates
    cmd="$gatk_Launcher --java-options '-Xms4000m -Djava.io.tmpdir=$TMPDIR -Dsamjdk.compression_level=2' MarkDuplicates \
        -I ${output_prefix}sorted.bam \
        -O ${output_prefix}deduplicated.bam \
        -M ${output_prefix}duplication.metrics \
        --REMOVE_DUPLICATES true \
        --CREATE_INDEX true"
    run "$cmd" "Dedup"
end_time_mod=$(date +%s)
if [[ "$OSTYPE" == "darwin"* ]]; then start_date=$(date -j -f "%s" $start_time_mod); else start_date=$(date -d @$start_time_mod); fi
if [[ "$OSTYPE" == "darwin"* ]]; then end_date=$(date -j -f "%s" $end_time_mod); else end_date=$(date -d @$end_time_mod); fi
echo "Stage "$stage" Started: "$start_date"; Ended: "$end_date"; Elapsed time: "$(($end_time_mod - $start_time_mod))" sec"

stage="BQSR"
start_time_mod=$(date +%s)
    # BQSR scatter
    bqsr_block=1
    input_bqsr_reports=""
    while read -r line
    do      
      cmd="$gatk_Launcher --java-options '-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m -Djava.io.tmpdir=$TMPDIR' \
        BaseRecalibrator \
        -I ${output_prefix}deduplicated.bam \
        -R $fasta \
        --known-sites $known_Mills \
        --known-sites $known_1000GIndels \
        --known-sites $known_dbsnp \
        $line \
        -O $TMPDIR/${sample_id}_bqsr_${bqsr_block}.grp"
      run "$cmd" "Bqsr_table_${bqsr_block}" &
      input_bqsr_reports="$input_bqsr_reports -I $TMPDIR/${sample_id}_bqsr_${bqsr_block}.grp"
      bqsr_block=$((bqsr_block+1))
    done < $TMPDIR/sequence_grouping.txt
    wait

    # BQSR gather
    cmd="$gatk_Launcher --java-options '-Xms3000m -Djava.io.tmpdir=$TMPDIR' \
        GatherBQSRReports \
        $input_bqsr_reports \
        -O ${output_prefix}bqsr.grp"
    run "$cmd" "GatherBQSRReports"

    # ApplyBQSR scatter
    bqsr_block=1
    input_recalibrated_bams=""
    while read -r line
    do      
      cmd="$gatk_Launcher --java-options '-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m -Djava.io.tmpdir=$TMPDIR' \
        ApplyBQSR \
        -R $fasta \
        -I ${output_prefix}deduplicated.bam \
        --bqsr-recal-file ${output_prefix}bqsr.grp \
        $line \
        -O $TMPDIR/${sample_id}_${bqsr_block}_recaled.bam"
      run "$cmd" "ApplyBQSR_${bqsr_block}" &
      input_recalibrated_bams="$input_recalibrated_bams -I $TMPDIR/${sample_id}_${bqsr_block}_recaled.bam"
      bqsr_block=$((bqsr_block+1))
    done < $TMPDIR/sequence_grouping_with_unmapped.txt
    wait

    # ApplyBQSR gather
    cmd="$gatk_Launcher --java-options '-Xms3000m -Djava.io.tmpdir=$TMPDIR -Dsamjdk.compression_level=2' \
        GatherBamFiles \
        $input_recalibrated_bams \
        --CREATE_INDEX true \
        -O ${output_prefix}recalibrated.bam"
    run "$cmd" "GatherBamFiles"      
end_time_mod=$(date +%s)
if [[ "$OSTYPE" == "darwin"* ]]; then start_date=$(date -j -f "%s" $start_time_mod); else start_date=$(date -d @$start_time_mod); fi
if [[ "$OSTYPE" == "darwin"* ]]; then end_date=$(date -j -f "%s" $end_time_mod); else end_date=$(date -d @$end_time_mod); fi
echo "Stage "$stage" Started: "$start_date"; Ended: "$end_date"; Elapsed time: "$(($end_time_mod - $start_time_mod))" sec"

stage="HaplotypeCaller"
start_time_mod=$(date +%s)
    # HaplotypeCaller scatter
    hc_block=1
    variant_file=""
    while read -r line
    do 
      cmd="$gatk_Launcher --java-options '-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=$TMPDIR' \
          HaplotypeCaller \
          -R $fasta \
          -I ${output_prefix}recalibrated.bam \
          $line \
          --dbsnp $known_dbsnp \
          -O $TMPDIR/${sample_id}_${hc_block}.g.vcf.gz \
          -ERC GVCF"
      run "$cmd" HC_${hc_block} &
      variant_file="$variant_file -I $TMPDIR/${sample_id}_${hc_block}.g.vcf.gz"
      hc_block=$((hc_block+1))
    done <$TMPDIR/ScatterIntervalList.txt 
    wait

    cmd="$gatk_Launcher --java-options '-Xms2000m -Djava.io.tmpdir=$TMPDIR' \
        MergeVcfs \
        $variant_file \
        -O ${output_prefix}output.g.vcf.gz"
    run "$cmd" "MergeVcfs"
end_time_mod=$(date +%s)
if [[ "$OSTYPE" == "darwin"* ]]; then start_date=$(date -j -f "%s" $start_time_mod); else start_date=$(date -d @$start_time_mod); fi
if [[ "$OSTYPE" == "darwin"* ]]; then end_date=$(date -j -f "%s" $end_time_mod); else end_date=$(date -d @$end_time_mod); fi
echo "Stage "$stage" Started: "$start_date"; Ended: "$end_date"; Elapsed time: "$(($end_time_mod - $start_time_mod))" sec"

stage="GenotypeGVCFs"
    start_time_mod=$(date +%s)
    cmd="$gatk_Launcher --java-options '-Xms4000m -Djava.io.tmpdir=$TMPDIR' \
        GenotypeGVCFs \
        -R $fasta \
        --dbsnp $known_dbsnp \
        -V ${output_prefix}output.g.vcf.gz \
        -O ${output_prefix}output.vcf.gz"
    run "$cmd" "GenotypeGVCFs"
    end_time_mod=$(date +%s)
if [[ "$OSTYPE" == "darwin"* ]]; then start_date=$(date -j -f "%s" $start_time_mod); else start_date=$(date -d @$start_time_mod); fi
if [[ "$OSTYPE" == "darwin"* ]]; then end_date=$(date -j -f "%s" $end_time_mod); else end_date=$(date -d @$end_time_mod); fi
echo "Stage "$stage" Started: "$start_date"; Ended: "$end_date"; Elapsed time: "$(($end_time_mod - $start_time_mod))" sec"

rm -r $TMPDIR
