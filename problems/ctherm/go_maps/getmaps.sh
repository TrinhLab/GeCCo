file="Ctherm_Pangenome.base_genome-DSM1313.core_pangenome.FeatureSet.tsv"
# Files must be sorted for join
sed -nE 's/.*(Cthe_[0-9]+).*(WP_[0-9.]+).*/\2,\1/p' $file | sort > /tmp/wp2cthe.csv
sed -nE 's/(CLO1313_RS[0-9]+).*(WP_[0-9.]+).*/\2,\1/p' $file | sort > /tmp/wp2clo.csv
join -t , -j 1  /tmp/wp2clo.csv /tmp/wp2cthe.csv > idmap.csv

#info
echo "Cthe matches: $(cat /tmp/wp2cthe.csv | wc -l)"
echo "CLO matches: $(cat /tmp/wp2clo.csv | wc -l)"
echo "Overlapping matches: $(cat idmap.csv | wc -l)"
