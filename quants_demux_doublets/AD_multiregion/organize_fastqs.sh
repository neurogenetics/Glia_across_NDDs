FASTQ_LIST=$(cat /data/ADRD/AD_multiregion/fastq_processing/cleaned_AD_fastqs.txt)

for FASTQ in ${FASTQ_LIST}
do
mkdir -p /data/ADRD/AD_multiregion/fastq_processing/all_fastqs/$FASTQ
mv /data/ADRD/AD_multiregion/fastqs/'Sequencing Round 1'/$FASTQ* /data/ADRD/AD_multiregion/fastq_processing/all_fastqs/$FASTQ/
mv /data/ADRD/AD_multiregion/fastqs/'Sequencing Round 2'/$FASTQ* /data/ADRD/AD_multiregion/fastq_processing/all_fastqs/$FASTQ/
done

# for samples which were sequenced more than once, need to separate
for dir in */
do
base_dir=$(basename "$dir")
S_nums=($(ls "$dir" | grep -oP '(?<=_S)\d+' | sort -n | uniq))
if [ ${#S_nums[@]} -gt 1 ]; then
  highest_S=${S_nums[-1]}
  new_dir="${base_dir}-2"
  mkdir -p "$new_dir"
  for file in "$dir"/*_S${highest_S}*
  do
  mv "$file" "$new_dir/"
  done
fi
done
