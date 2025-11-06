stat_file=$1
out_file=$2
awk '{print $1}' $stat_file | sed 's/*//g'|sed '/^[[:blank:]]*$/d' | head -25 > $out_file
