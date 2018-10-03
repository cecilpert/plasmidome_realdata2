set -e 

function usage(){
	echo "usage : bash draw_devenn.sh <fasta assembly> <file 1> <file 2> <prefix> <name 1> <name 2> <color 1> <color 2> <output prefix>"
}	

if [[ "$#" -ne 9 ]]; then 
	usage
	exit 1 
fi 

assembly=$1
f1=$2
f2=$3
pref=$4
name1=$5
name2=$6
c1=$7
c2=$8
out=$9

outstat=$out.stats
outlength=$out.length.pdf
outnumber=$out.number.pdf

intersect_pref=$pref.$name1.$name2
intersect_common=results/intersect/$intersect_pref.common.txt

BIN=/databis/hilpert/plasmidome_project/bin
BIN2=/databis/hilpert/plasmidome_realdata2/bin

python3 $BIN2/intersect_files.py $f1 $f2 results/intersect $intersect_pref $name1 $name2 	

echo -e "Part\tNumber.contigs\tLength" > $outstat


echo "# Compute $f1 total length" 
a1=$(python3 $BIN/total_length_contig_list.py $assembly $f1) 
echo "# Compute $f2 total length" 
a2=$(python3 $BIN/total_length_contig_list.py $assembly $f2)
echo "# Compute $f1 & $f2 common length" 
a3=$(python3 $BIN/total_length_contig_list.py $assembly $intersect_common) 
	
	
echo "# Compute $f1 number contigs" 
n1=$(wc -l $f1 | cut -f 1 -d " ") 
echo "# Compute $f2 number contigs" 
n2=$(wc -l $f2 | cut -f 1 -d " ")
echo "# Compute $f1 & $f2 common number contigs" 
n3=$(wc -l $intersect_common | cut -f 1 -d " ")
 	
echo -e "$name1\t$n1\t$a1" >> $outstat
echo -e "$name2\t$n2\t$a2" >> $outstat
echo -e "$name1.$name2.common\t$n3\t$a3" >> $outstat

perc1=$(echo $a1 $a3 | awk '{print $2/$1}')  
perc2=$(echo $a2 $a3 | awk '{print $2/$1}')
percn1=$(echo $n1 $n3 | awk '{print $2/$1}') 
percn2=$(echo $n2 $n3 | awk '{print $2/$1}') 

echo -e "Common/$name1\t$percn1\t$perc1" >> $outstat
echo -e "Common/$name2\t$percn2\t$perc2" >> $outstat 
 
echo "# Draw graph" 
Rscript --vanilla $BIN2/draw_devenn.R $a1 $a2 $a3 $name1 $name2 $c1 $c2 $outlength $pref 2> log/$pref.R
Rscript --vanilla $BIN2/draw_devenn.R $n1 $n2 $n3 $name1 $name2 $c1 $c2 $outnumber $pref 2> log/$pref.R
