### Usage example: sh concat_bigscape.sh /Volumes/TFL210426/bigscape/bigscape_output/ 2022-01-11_12-53-12_hybrids_glocal/

cd $1
rm $1bigscape_all_c030.txt
find $1network_files/$2 -type f -name *_c0.30.network | cat > $1network_files.txt
while read line; do cat $line >> $1pre_bigscape_all_c030.txt; done < $1network_files.txt
awk '/Clustername/&&c++>0 {next} 1' $1pre_bigscape_all_c030.txt >> $1bigscape_all_c030.txt