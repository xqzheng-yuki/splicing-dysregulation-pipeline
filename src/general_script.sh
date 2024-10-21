wget -O /mnt/gtklab01/xiaoqing/salmon/decoy/index/salmon_partial_sa_index.tgz http://refgenomes.databio.org/v3/assets/archive/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1/salmon_partial_sa_index

# common used directory
cd /mnt/gtklab01/xiaoqing
cd /mnt/gtklab01/xiaoqing/salmon/full_output
cd /mnt/gtklab01/linglab/tdp43/fastq

cat > sample_name.txt
nohup bash salmon.sh > output.log 2> error.log &
# mv /path/to/my_folder /path/to/new_location/new_name

more unmapped_names.txt | grep "d" > decoy_mapping.txt

READS_FILE="/mnt/gtklab01/xiaoqing/sample_name.txt"
RESULT_DIR="/mnt/gtklab01/xiaoqing/salmon/full_output"
NAME_SUF="aux_info/unmapped_names.txt"
while read SAMPLE_NAME; do
    cut -d" " -f 2 $RESULT_DIR/$SAMPLE_NAME/$NAME_SUF | sort | uniq -c
done < $READS_FILE