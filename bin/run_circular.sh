set -e 

function usage(){ 
		echo "usage : bash run_circular.sh <input fasta> <output fasta>" 
}

if [[ "$#" -ne 2 ]]; then 
	usage
	exit 1 
fi 

tmp=`mktemp -d -p .`
BIN=/databis/hilpert/plasmidome_project/bin


awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' $1 > $tmp/oneline.fasta
/data/chochart/SCRIPT/CIRCULAR/detect.circular.seq.pl -f $tmp/oneline.fasta -k 50 > $tmp/circular.fasta	
grep "circ" $tmp/circular.fasta | cut -f 1 -d " " | sed 's/>//g' > $tmp/circular.id
python3 $BIN/seq_from_list.py --input_fasta $1 --keep $tmp/circular.id --output_fasta $2 

rm -r $tmp 
