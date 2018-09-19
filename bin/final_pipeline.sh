set -e 

function usage(){
	echo 'usage : all_analysis.sh -a assembly.fasta -o outdir 
	[OPTIONS]
	-h : print help
	#Output  
	--prefix <prefix> : prefix for output files (default : assembly fasta name) 
	#Databases 
	--chrm_db <fasta> : fasta file with bacterial chromosomes sequences you want to use (default : databases/all_chrm.fasta) 
	--rna_db <fasta> : fasta file with rRNA sequences (must be back transcribed) you want to use (default : databases/rRNA/SILVA_132_LSUParc_SSUParc_tax_silva_trunc.T.fasta) 
	--phylo_db <hmm> : hmm profile(s) with phylogenetic markers (default : databases/phylogenetic_markers/wu2013/bacteria_and_archaea_dir/BA.hmm)
	' 
}

function treat_args(){
	if [[ ! $assembly ]]; then 
		echo "No assembly.fasta given. Use -a option." 
		exit 1 
	fi 
	
	if [[ ! $outdir ]]; then 
		echo "No output directory given. Use -o option" 
		exit 1 
	fi	
}	

function set_default(){
	if [[ ! $PREFIX ]]; then
		prefix=$(echo $assembly | rev | cut -f 1 -d "/" | cut -f 2- -d "." | rev) 
	fi
	if [[ ! $CHRM_DB ]]; then 
		chrm_db=databases/all_chrm.fasta
	fi
	if [[ ! $RNA_DB ]]; then 
		rna_db=databases/rRNA/SILVA_132_LSUParc_SSUParc_tax_silva_trunc.T.fasta
	fi
	if [[ ! $PHYLO_DB ]]; then 
		phylo_db=databases/phylogenetic_markers/wu2013/bacteria_and_archaea_dir/BA.hmm
	fi	
}	

function verif_file(){
	f=$1 
	message=$2
	message2=$3
	if [[ ! -f $f ]]; then 
		echo $message 
		exit 1
	else 
		echo $message2	
	fi
}	

function verif_result(){
	if [[ -f $1 && ! $FORCE ]]; then
		file_exist=1
	else 
		file_exist=0	
	fi 	
}

function define_paths(){ 
	all_contigs=$assembly 
	contigs_nochrm=$outdir/intermediate_contigs/$prefix.nochrm.fasta
	contigs_norna=$outdir/intermediate_contigs/$prefix.nochrm.norna.fasta
	contigs_nophylomark=$outdir/intermediate_contigs/$prefix.nochrm.norna.nophylomark.fasta
	predicted_proteins=$outdir/protein_prediction/$prefix.predicted_proteins.faa 
	contigs_learning=$outdir/intermediate_contigs/$prefix.nochrm.norna.nophylomark.learning.fasta 
	contigs_circular=$outdir/intermediate_contigs/$prefix.circular.fasta
		
}	

function complete_resume(){
	name=$1
	nb_contigs=$2
	length_contigs=$3	
	path=$4
	echo -e "$name\t$nb_contigs\t$length_contigs\t$path" >> $resume 
}	

function verif_blast_db(){
	file=$1 
	type=$2
	message_create=$3
	message_found=$4
	
	if [[ $type == "prot" ]]; then 
		local pref="phr" 	
	elif [[	$type == "nucl" ]]; then 
		local pref="nhr" 	
	fi
	
	if [ ! -f $file.$pref ];then 
		echo $message_create
		makeblastdb -in $file -dbtype $type &> $outdir/log/rna_search.blast_db.log
	else	
		echo $message_found	
	fi
}

function treat_blast_rna(){
	blast=$1
	awk -F "\t" '{if ($3 >= 80 && $4 >= 300) print}' $blast.tsv > $blast.id80.length300.tsv
	awk -F "\t" '{if ($4 >= 1200 && $15 >= 80) print}' $blast.id80.length300.tsv > $blast.id80.length1200.cov80.tsv
	awk -F "\t" '{if (($7 == 1 || $7 == $14 || $8 == 1 || $8 == $14) && ($9 == 1 || $9 == $16 || $10 == 1 || $10 == $16)) print}' $blast.id80.length300.tsv > $blast.id80.length300.ends.tsv
	cut -f 2 $blast.id80.length1200.cov80.tsv | sort -u > $blast.id80.length1200.cov80.contigs.id 
	cut -f 2 $blast.id80.length300.ends.tsv | sort -u > $blast.id80.length300.ends.contigs.id
	cat $blast.id80.length1200.cov80.contigs.id $blast.id80.length300.ends.contigs.id | sort -u > $blast.all.contigs.id
	
}	

function treat_hmm(){
	hmm=$1 
	grep -v "^#" $hmm.tsv | awk '{print $1}' | sort -u > $hmm.predicted_proteins.id 
	grep -v "^#" $hmm.tsv | awk '{print $3}' | sort -u > $hmm.proteins_family.id
	grep "^k" $hmm.predicted_proteins.id | rev | cut -f 2- -d "_" | rev > $hmm.contigs.id 
}

function run_align_chrm(){  
	mkdir -p $outdir/align_chrm
	out=$outdir/align_chrm/$prefix.align_chrm
	verif_file $chrm_db "ERROR Chromosomes database doesn't exist. Use --chrm_db to specify it." "[align_chrm] Chromosomes databases found in $chrm_db"
	echo "[align_chrm] Use $current_assembly" 
	echo "[align_chrm] Count chromosomes..." 
	nb_chrm=$(grep "^>" -c $chrm_db) 
	echo "[align_chrm] Run minimap2..." 
	source activate plasmidome 
	minimap2 -x asm5 -N $nb_chrm $chrm_db $assembly 1> $out.paf 2> $outdir/log/align_chrm.log
	source deactivate plasmidome
	echo "[align_chrm] Treat minimap2..." 
	awk '{if ($10/$2 >= 0.8) print}' $out.paf > $out.conserve.paf 
	cut -f 1 $out.conserve.paf | sort -u > $out.conserve.contigs.id 
	echo "[align_chrm] Delete chromosomic contigs..."
	python3 $BIN/delete_seq_from_file.py $current_assembly $out.conserve.contigs.id $contigs_nochrm normal 
}	

function run_rna_search(){
	mkdir -p $outdir/rna_search
	out=$outdir/rna_search/$prefix.rna
	verif_file $rna_db "ERROR rRNA database doesn't exist. Use --rna_db to specify it." "[rna_search] rRNA databases found in $rna_db"
	echo "[rna_search] Use $current_assembly" 
	verif_blast_db $current_assembly nucl "[rna_search] Make blast db for $current_assembly..." "[rna_search] Blast db found for $current_assembly"
	echo "[rna_search] Run Blast..." 
	bash $BIN2/parallelize_blast.sh $rna_db $current_assembly $out.blast.tsv 32 blastn
	echo "[rna_search] Treat Blast..." 
	treat_blast_rna $out.blast 
	echo "[rna_search] Delete rRNA contigs..."
	python3 $BIN/delete_seq_from_file.py $current_assembly $out.blast.all.contigs.id $contigs_norna normal
}	

function run_phylo_markers_search(){
	mkdir -p $outdir/phylogenetic_markers_search 
	out=$outdir/phylogenetic_markers_search/$prefix.phylo_markers.hmm
	verif_file $phylo_db "ERROR Phylogenetic markers HMM profile doesn't exist. Use --phylo_db to specify it" "[phylo_search] HMM profiles found in $phylo_db" 
	echo "[phylo_search] Use $predicted_proteins" 
	echo "[phylo_search] Run Hmmer..." 
	hmmsearch --tblout $out.tsv -E 1e-5 $phylo_db $predicted_proteins > $out
	treat_hmm $out
	echo "[phylo_search] Delete phylogenetic markers contigs..." 
	python3 $BIN/delete_seq_from_file.py $current_assembly $out.contigs.id $contigs_nophylomark normal
	
}	

function run_proteins_prediction(){
	mkdir -p $outdir/protein_prediction 
	out=$outdir/protein_prediction/$prefix.predicted_proteins
	echo "[protein_prediction] Use $current_assembly" 
	echo "[protein_prediction] Run Prodigal..."
	prodigal -i $current_assembly -c -m -p meta -f gff -a $out.faa -o $out.gff -q
}	

function run_learning(){
	mkdir -p $outdir/learning 
	echo "[learnig] Use $current_assembly" 
	echo "[learning] Run PlasFlow..." 
	bash $BIN/all_decontamination.sh -f $current_assembly -o $outdir/learning --plasflow 70 --prefix $prefix > $outdir/log/plasmid_prediction.log 
	mv $outdir/learning/plasflow/$prefix.plasflow0.7.plasmids.fasta $contigs_learning
}

function run_circular(){
	mkdir -p $outdir/circular 
	echo "[circular] Use $current_assembly" 
	out=$outdir/circular/$prefix.circular
	echo "[circular] Run circularization script..." 
	bash $BIN2/run_circular.sh $current_assembly $out.fasta
	mv $out.fasta $contigs_circular
	
}		

TEMP=$(getopt -o h,a:,o: -l prefix:,chrm_db:,rna_db:,phylo_db:,force  -- "$@")
eval set -- "$TEMP" 
while true ; do 
	case "$1" in 
		-a) 
			assembly=$2
			shift 2;; 
		-o) 
			outdir=$2
			shift 2;; 
		--prefix) 
			PREFIX=1
			prefix=$2 
			shift 2;; 
		--chrm_db) 
			CHRM_DB=1
			chrm_db=$2
			shift 2;;	
		--rna_db)
			RNA_DB=1
			rna_db=$2
			shift 2;;
		--phylo_db)
			PHYLO_DB=1 
			phylo_db=$2
			shift 2;; 
		--force)
			FORCE=1
			shift ;; 
		-h) 
			usage 
			shift ;;
		--)  
			shift ; break ;; 					
	esac 
done

BIN=/databis/hilpert/plasmidome_project/bin
BIN2=/databis/hilpert/plasmidome_realdata2/bin

treat_args
set_default
define_paths

mkdir -p $outdir 
mkdir -p $outdir/log 
mkdir -p $outdir/intermediate_contigs
resume=$outdir/resume.tsv

echo -e "step\tcontigs_number\tcontigs_length\tpath" > $resume 
current_assembly=$all_contigs 
complete_resume "All contigs" $(grep "^>" -c $current_assembly) $(python3 $BIN/total_length_fasta.py $current_assembly) $current_assembly

echo "==== CLEANING PART ====" 
echo "## STEP 1 : CHROMOSOMES ALIGNMENT" 
start=$(date +%s) 
verif_result $contigs_nochrm 
if [[ $file_exist == 1 ]]; then 
	echo "Chromosomes alignment results already exists. Use --force to overwrite" 
else 
	run_align_chrm 
fi
complete_resume "Contigs without chromosomes" $(grep "^>" -c $contigs_nochrm) $(python3 $BIN/total_length_fasta.py $contigs_nochrm) $(readlink -f $contigs_nochrm)
current_assembly=$contigs_nochrm
end=$(date +%s)
echo "Time elapsed : $((end-start)) s" 

echo "## STEP 2 : rRNA SEARCH" 
start=$(date +%s)
verif_result $contigs_norna
if [[ $file_exist == 1 ]]; then 
	echo "rRNA search results already exists. Use --force to overwrite" 
else 
	run_rna_search
fi
complete_resume "Contigs without chromosomes and RNA" $(grep "^>" -c $contigs_norna) $(python3 $BIN/total_length_fasta.py $contigs_norna) $(readlink -f $contigs_norna)
current_assembly=$contigs_norna
end=$(date +%s)
echo "Time elapsed : $((end-start)) s" 

echo "## STEP 3 : PHYLOGENETIC MARKERS SEARCH" 
start=$(date +%s)
echo "# STEP 3.1 : Proteins prediction" 
verif_result $predicted_proteins 
if [[ $file_exist == 1 ]]; then 
	echo "Proteins prediction results already exists. Use --force to overwrite" 
else 
	run_proteins_prediction 
fi
echo "# STEP 3.2 : Phylogenetic markers search" 
verif_result $contigs_nophylomark 
if [[ $file_exist == 1 ]]; then 
	echo "Phylogenetic markers search results already exists. Use --force to overwrite" 
else 
	run_phylo_markers_search 
fi
complete_resume "Contigs without chromosomes, RNA and phylogenetic markers" $(grep "^>" -c $contigs_nophylomark) $(python3 $BIN/total_length_fasta.py $contigs_nophylomark) $(readlink -f $contigs_nophylomark)
current_assembly=$contigs_nophylomark
end=$(date +%s)
echo "Time elapsed : $((end-start)) s" 

echo ""
echo "==== PLASMID PREDICTION PART ====" 

echo "## STEP 1 : PLASMID PREDICTION FROM CLEANED CONTIGS WITH LEARNING" 
start=$(date +%s)
verif_result $contigs_learning 
if [[ $file_exist == 1 ]]; then 
	echo "Plasmid prediction results already exists. Use --force to overwrite" 
else 
	run_learning
fi
complete_resume "Predict plasmids (from cleaned contigs)" $(grep "^>" -c $contigs_learning) $(python3 $BIN/total_length_fasta.py $contigs_learning) $(readlink -f $contigs_learning)
end=$(date +%s)
echo "Time elapsed : $((end-start)) s" 

echo "## STEP 2 : CIRCULAR CONTIGS"
start=$(date +%s)
current_assembly=$all_contigs
verif_result $contigs_circular  
if [[ $file_exist == 1 ]]; then 
	echo "Circular contigs results already exists. Use --force to overwrite" 
else 
	run_circular
fi
complete_resume "Circular contigs" $(grep "^>" -c $contigs_circular) $(python3 $BIN/total_length_fasta.py $contigs_circular) $(readlink -f $contigs_circular)
end=$(date +%s)
echo "Time elapsed : $((end-start)) s" 

