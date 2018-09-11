set -e 

function usage(){
	echo "usage : bash draw_devenn.sh <fasta assembly> <file 1> <file 2> <prefix> <name 1> <name 2> <color 1> <color 2> <graph output> <stats output> <mode>
	
	<mode> : length|number : length will compute with length of contigs and number with number of contigs"
}	

if [[ "$#" -ne 11 ]]; then 
	usage
	exit 1 
fi 

assembly=$1
f1=$2
f2=$3
pref=$4
n1=$5
n2=$6
c1=$7
c2=$8
out=$9
outstat=${10}
mode=${11}

intersect_pref=$pref.$n1.$n2
intersect_common=results/intersect/$intersect_pref.common.txt

BIN=/databis/hilpert/plasmidome_project/bin
BIN2=/databis/hilpert/plasmidome_realdata2/bin

python3 $BIN2/intersect_files.py $f1 $f2 results/intersect $intersect_pref $n1 $n2 	

if [[ $mode == "length" ]]; then 
	echo -e "$n1\t$n2\t$n1.$n2" > $outstat
	echo "# Compute $f1 total length" 
	a1=$(python3 $BIN/total_length_contig_list.py $assembly $f1) 
	echo "# Compute $f2 total length" 
	a2=$(python3 $BIN/total_length_contig_list.py $assembly $f2)
	echo "# Compute $f1 & $f2 common length" 
	a3=$(python3 $BIN/total_length_contig_list.py $assembly $intersect_common)
	
	echo -e "$a1\t$a2\t$a3" >> $outstat  
	
	
elif [[ $mode == "number" ]]; then 
	echo "# Compute $f1 number contigs" 
	a1=$(wc -l $f1 | cut -f 1 -d " ") 
	echo "# Compute $f2 number contigs" 
	a2=$(wc -l $f2 | cut -f 1 -d " ")
	echo "# Compute $f1 & $f2 common number contigs" 
	a3=$(wc -l $intersect_common | cut -f 1 -d " ")
fi 	

echo "# Draw graph" 
Rscript --vanilla $BIN2/draw_devenn.R $a1 $a2 $a3 $n1 $n2 $c1 $c2 $out $pref 2> log/$pref.R
