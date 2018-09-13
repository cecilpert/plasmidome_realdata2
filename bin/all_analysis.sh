function usage(){
	echo 'usage : all_analysis.sh <prefix> 
	[options] : 
	-p : prefix 
	--plasmid_align
	--learning : launch learning step to predict plasmids
	--learning_option <stringent|normal|both> stringent launch plasflow with a high decision probability so you will be relatively sure of the prediction but you will loose a lot of plasmids too. Normal launch a combinaison of cBar and PlasFlow which try to optimize plasmids discovery and plasmids loss. Both launch the two options, with two results (default : normal) 
	--rna_search : launch rna_search step 
	--rna_search_option <blast|barrnap|both> :  Barrnap use hmmer databases (https://github.com/tseemann/barrnap) and Blast will launch blast alignment against SSU and LSU SILVA databases. Blast is longer to execute. Both launch the two options. (default : barrnap)
	--rna_db : path to rna db directory if using blast rna search (default : databases/rRNA) 
	--chrm_align
	--plasmid_db <plasmid db> (default : databases/all_plasmids.fasta)
	--chrm_db <chrm db> (default : databases/all_chrm.fasta) 
	--graphs : launch graphical representation of different metrics
	--graph_rna : launch graphical representation of comparison between plasmid predict contigs and contigs with rRNA
	--search_plasmids_markers : search plasmidic proteins (replicon, mob, mpf, orit) in assembly 
	-h : print help
	' 
}

TEMP=$(getopt -o h,p: -l plasmid_align,plasmid_db:,chrm_align,chrm_db:,learning,learning_option:,graphs,rna_search,rna_db:,rna_search_option:,graph_rna,search_plasmids_markers -- "$@")
eval set -- "$TEMP" 
while true ; do 
	case "$1" in 
		-p) 
			prefix=$2
			shift 2;; 
		--plasmid_align) 
			PLASMID_ALIGN=1
			shift;; 
		--chrm_align)
			CHRM_ALIGN=1 
			shift ;; 	
		--learning)
			LEARNING=1
			shift ;; 
		--learning_option) 
			learning_option=$2 
			shift 2;;
		--plasmid_db) 
			plasmid_db=$2
			shift 2;;
		--chrm_db) 
			chrm_db=$2
			shift 2;;	
		--graphs)
			GRAPHS=1
			shift;;
		--rna_search) 
			RNA=1 
			shift;;
		--rna_search_option)
			rna_option=$2
			shift 2;; 
		--rna_db) 
			rna_db=$2
			shift 2;; 
		--graph_rna)
			GRAPH_RNA=1
			shift;;	
		--search_plasmids_markers) 
			PLASMIDS_MARK=1 
			shift;; 	
		-h) 
			usage 
			shift ;;
		--)  
			shift ; break ;; 					
	esac 
done

BIN2=/databis/hilpert/plasmidome_realdata2/bin
circular=results/circular/megahit.$prefix.1kb.circular.id 
assembly=results/$prefix\_assembly/megahit.$prefix.1kb.fasta
align_plasmids=results/align_plasmids/$prefix.contig80.id
align_chrm=results/align_chrm/$prefix.contig80.id
learning_plasmids=results/learning/plasflow/$prefix.plasflow0.7.plasmids.id 
learning_plasmids_nochrm=results/intersect/$prefix.learning_plasmids.align_chrm.learning_plasmids.specific.txt
all_markers=results/plasmids_markers/$prefix.all_markers.contigs.id 
rep=results/plasmids_markers/rep/$prefix.rep.blast.id80.cov80.contigs.id 
orit=results/plasmids_markers/orit/$prefix.orit.blast.id90.cov90.contigs.id
mob=results/plasmids_markers/mob/$prefix.mob.predicted_proteins.blast.id80.cov80.contigs.id
mpf=results/plasmids_markers/mpf/$prefix.mpf.predicted_proteins.blast.id80.cov80.contigs.id
rna=results/rna_search/blast/new_sorting/$prefix.allRNA.blast.all.contigs.id

echo "oooo"
mkdir -p results/graphs_learning 

echo "[graph] Learning_plasmids_nochrm & RNA" 
graph=results/graphs_learning/$prefix.learning_plasmids_nochrm.rna
bash $BIN2/draw_devenn.sh $assembly $learning_plasmids_nochrm $rna $prefix learning_plasmids_nochrm contigs_with_rna blue red $graph

exit

echo "[graph] Learning_plasmids & Circular" 
graph=results/graphs_learning/$prefix.learning_plasmids.circular
bash $BIN2/draw_devenn.sh $assembly $learning_plasmids $circular $prefix learning_plasmids circular blue green $graph

echo "[graph] Align_plasmids & Align_chrm" 
graph=results/graphs_learning/$prefix.align_plasmids.align_chrm
bash $BIN2/draw_devenn.sh $assembly $align_plasmids $align_chrm $prefix align_plasmids align_chrm blue red $graph

echo "[graph] Learning_plasmids & Align_chrm" 
graph=results/graphs_learning/$prefix.learning_plasmids.align_chrm
bash $BIN2/draw_devenn.sh $assembly $learning_plasmids $align_chrm $prefix learning_plasmids align_chrm blue red $graph

echo "[graph] Learning_plasmids & Align_plasmids" 
graph=results/graphs_learning/$prefix.learning_plasmids.align_plasmids
bash $BIN2/draw_devenn.sh $assembly $learning_plasmids $align_plasmids $prefix learning_plasmids align_plasmids blue green $graph

echo "[graph] Learning_plasmids_nochrm & All_markers" 
graph=results/graphs_learning/$prefix.learning_plasmids_nochrm.all_plasmids_markers
bash $BIN2/draw_devenn.sh $assembly $learning_plasmids_nochrm $all_markers $prefix learning_plasmids_no_chrm all_plasmids_markers blue green $graph

echo "[graph] Learning_plasmids_nochrm & REP" 
graph=results/graphs_learning/$prefix.learning_plasmids_nochrm.rep
bash $BIN2/draw_devenn.sh $assembly $learning_plasmids_nochrm $rep $prefix learning_plasmids_no_chrm contigs_with_rep blue green $graph

echo "[graph] Learning_plasmids_nochrm & mob" 
graph=results/graphs_learning/$prefix.learning_plasmids_nochrm.mob
bash $BIN2/draw_devenn.sh $assembly $learning_plasmids_nochrm $mob $prefix learning_plasmids_no_chrm contigs_with_mob blue green $graph

echo "[graph] Learning_plasmids_nochrm & mpf" 
graph=results/graphs_learning/$prefix.learning_plasmids_nochrm.mpf
bash $BIN2/draw_devenn.sh $assembly $learning_plasmids_nochrm $mpf $prefix learning_plasmids_no_chrm contigs_with_mpf blue green $graph

echo "[graph] Learning_plasmids_nochrm & orit" 
graph=results/graphs_learning/$prefix.learning_plasmids_nochrm.orit
bash $BIN2/draw_devenn.sh $assembly $learning_plasmids_nochrm $orit $prefix learning_plasmids_no_chrm contigs_with_orit blue green $graph
