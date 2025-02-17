count="1"
file=$1
while read -r line; do

echo count $count

input_filename=$(echo $line | awk '{print $1}')

echo $input_filename

sh SubmitToCondor.sh output_${1}_${count} ${input_filename}

count=$[$count+1]

done < $file

exit 0

