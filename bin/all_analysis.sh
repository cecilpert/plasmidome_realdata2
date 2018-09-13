function usage(){
	echo 'usage : all_analysis.sh -p prefix
	[options] : 
	--plasmid_align
	--learning : launch learning step to predict plasmids
	--learning_option <cbar_plaslow|plasflow70|plasflow80> (default : plasflow70)
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

function treat_args(){
	if [[ ! $prefix ]]; then
		usage
		echo "You must give prefix"
		exit 1 
	fi	
	if [[ $learning_option != "cbar_plasflow" && $learning_option != "plasflow70" && $learning_option != "plasflow80" ]]; then 
		usage 
		echo "--learning_option must be cbar_plaslow, plasflow70 or plasflow80" 
		exit 1 
	fi	
	if [[ $rna_option != "barrnap" && $rna_option != "blast" && $rna_option != "both" ]]; then 
		usage 
		echo "--rna_search_option must be barrnap,blast or both" 
		exit 1 
	fi
}

function define_paths(){
	BIN2=/databis/hilpert/plasmidome_realdata2/bin
	BIN=/databis/hilpert/plasmidome_project/bin
	circular=results/circular/megahit.$prefix.1kb.circular.id 
	assembly=results/$prefix\_assembly/megahit.$prefix.1kb.fasta
	assembly_nochrm=results/$prefix\_assembly/megahit.$prefix.1kb.nochrm.fasta
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
}	

function verif_file(){
	f=$1 
	message=$2
	if [[ ! -f $f ]]; then 
		echo $message 
		exit 1
	fi
}	

function verif_result(){
	if [[ -f $1 && ! $FORCE ]]; then
		file_exist=1
	else 
		file_exist=0	
	fi 	
}	

function draw_devenn(){
	graph=$1
	f1=$2
	f2=$3
	name1=$4 
	name2=$5
	col1=$6
	col2=$7
	verif_result $graph.pdf 
	if [[ $file_exist == 1 ]]; then 
		echo "$graph already exists, use --force to overwrite"
	else 	
		bash $BIN2/draw_devenn.sh $assembly $f1 $f2 $prefix $name1 $name2 $col1 $col2 $graph	
	fi		
}	

learning_option="plasflow80"
rna_option="blast" 

TEMP=$(getopt -o h,p: -l plasmid_align,plasmid_db:,chrm_align,chrm_db:,learning,learning_option:,graphs,rna_search,rna_db:,rna_search_option:,graph_rna,search_plasmids_markers,force -- "$@")
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
		--force)
			FORCE=1	
			shift;;
		-h) 
			usage 
			shift ;;
		--)  
			shift ; break ;; 					
	esac 
done

treat_args
define_paths

mkdir -p results 

if [[ $CHRM_ALIGN ]]; then 
	echo "## ALIGN CHROMOSOMES STEP" 
	
	align_chrm_dir=results/align_chrm 
	mkdir -p $align_chrm_dir 
	
	bash $BIN/delete_seq_parallel.sh $assembly $align_chrm $assembly_nochrm normal 6
fi 


if [[ $LEARNING ]]; then 

	echo "## LEARNING STEP ($learning_option)" 
	verif_file $assembly_nochrm "Assembly file with no chromosomes doesn't exists. Launch --align_chrm before" 

	learning_dir=results/learning_nochrm
	mkdir -p $learning_dir
	
	if [ $learning_option == "cbar_plasflow" ]
	then 
		bash $BIN/all_decontamination.sh -f $assembly -o $learning_dir --cbar_plasflow 70 --prefix $prefix 
		
	elif [ $learning_option == "plasflow70" ]
	then
		bash $BIN/all_decontamination.sh -f $assembly -o $learning_dir --plasflow 70 --prefix $prefix
	elif [ $learning_option == "plasflow80" ] 
	then 
		bash $BIN/all_decontamination.sh -f $assembly -o $learning_dir --plasflow 80 --prefix $prefix
	fi
	
	
fi 
	

if [[ $GRAPHS ]]; then 

	graph_dir=results/new_graphs 
	mkdir -p $graph_dir
	
	graph_learning_dir=$graph_dir/$learning_option
	mkdir -p $graph_learning_dir
	
	echo "[graph] Align_plasmids & Align_chrm" 
	graph=$graph_dir/$prefix.align_plasmids.align_chrm
	draw_devenn $graph $align_plasmids $align_chrm Align.plasmids Align.chromosomes blue red

	echo "[graph] Learning_plasmids_nochrm & RNA" 
	graph=results/graphs_learning/$prefix.learning_plasmids_nochrm.rna
	bash $BIN2/draw_devenn.sh $assembly $learning_plasmids_nochrm $rna $prefix learning_plasmids_nochrm contigs_with_rna blue red $graph

	echo "[graph] Learning_plasmids & Circular" 
	graph=results/graphs_learning/$prefix.learning_plasmids.circular
	bash $BIN2/draw_devenn.sh $assembly $learning_plasmids $circular $prefix learning_plasmids circular blue green $graph

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

fi 
