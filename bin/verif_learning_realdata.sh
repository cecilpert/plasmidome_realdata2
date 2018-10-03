set -e 

function usage(){ 
	echo "bash verif_learning_realdata.sh -f <assembly.fasta> -o <outdir> --phylo_db <hmm profile> --rna_db <fasta> --chrm_db <fasta> --markers_db <dir> 
	Options : 
	--prefix <prefix> : prefix for results files, default : name of the fasta assembly
	--force : overwrite results if already exists
	--log <directory> : directory to put log files. Default : outdir
	--phylo_db <hmm profile> : phylogenetic markers hmm profile 
	--rna_db <fasta> : fasta file with 16S and 23S RNA 
	--chrm_db <fasta> : fasta file with chromosomes 
	--markers_db <dir> : directory with plasmids markers databases (must be named mob.proteins.faa, mpf.proteins.faa, orit.fas, rep.dna.fas) 
	--learning_thr <int (between 0 and 100)> : threshold for plasflow learning to detect plasmids (default : 80)
	--predicted_proteins <fasta> : predicted proteins for assembly. If not specified, will launch Prodigal to predict proteins  
	" 
}

function treat_args(){
	if [[ ! $assembly ]];then
		quit=1
		echo "You must give fasta assembly. Use -f option"
	fi		
	if [[ ! $outdir ]]; then 
		quit=1 
		echo "You must give output directory. Use -o option" 
	fi
	if [[ ! $phylo_db ]]; then 
		quit=1
		echo "You must give phylogenetic markers profiles. Use --phylo_db option" 
	fi
	if [[ ! $rna_db ]]; then 
		quit=1
		echo "You must give RNA database. Use --rna_db option" 
	fi
	if [[ ! $chrm_db ]]; then 
		quit=1
		echo "You must give chromosomes database. Use --chrm_db option" 
	fi
	if [[ ! $markers_db ]]; then 
		quit=1
		echo "You must give plasmids markers databases. Use --markers_db option" 
	fi
	if [[ $quit ]]; then 
		exit 1
	fi 
}	

function verif_args(){
	mob_db=$markers_db/mob.proteins.faa
	mpf_db=$markers_db/mpf.proteins.faa
	rep_db=$markers_db/rep.dna.fas
	orit_db=$markers_db/orit.fas
	verif_file $assembly "[verif_learning] Assembly file doesn't found in $assembly" "[verif_learning] Assembly file found in $assembly"
	verif_file $phylo_db  "[verif_learning] Phylogenetic markers database doesn't found in $phylo_db" "[verif_learning] Phylogenetic markers database found in $phylo_db"
	verif_file $rna_db  "[verif_learning] RNA database doesn't found in $rna_db" "[verif_learning] RNA database found in $rna_db"
	verif_file $chrm_db  "[verif_learning] Chromosomes database doesn't found in $chrm_db" "[verif_learning] Chromosomes database found in $chrm_db"
	verif_file $mob_db  "[verif_learning] MOB database doesn't found in $mob_db" "[verif_learning] MOB database found in $mob_db"
	verif_file $mpf_db  "[verif_learning] MPF database doesn't found in $mpf_db" "[verif_learning] MPF database found in $mpf_db"
	verif_file $rep_db  "[verif_learning] REP database doesn't found in $rep_db" "[verif_learning] REP database found in $rep_db"
	verif_file $orit_db  "[verif_learning] OriT database doesn't found in $orit_db" "[verif_learning] OriT database found in $orit_db"
	mkdir -p $outdir 
	if [[ ! $prefix ]]; then 
		prefix=$(echo $assembly | rev | cut -f 1 -d "/" | cut -f 2- -d "." | rev)
	fi

	if [[ ! $log ]]; then 
		log=$outdir/log	
	fi
	mkdir -p $log 	
}

learning_thr=70

TEMP=$(getopt -o h,f:,o: -l prefix:,force,log:,rna_db:,chrm_db:,phylo_db:,markers_db:,proteins:  -- "$@")
eval set -- "$TEMP" 
while true ; do 
	echo $1
	case "$1" in 
		-f)
			assembly=$2
			shift 2;; 
		-o)
			outdir=$2
			shift 2;; 
		--chrm_db) 
			chrm_db=$2
			shift 2;; 
		--phylo_db) 
			phylo_db=$2
			shift 2;; 	
		--markers_db) 
			markers_db=$2
			shift 2;; 	
		--rna_db) 
			rna_db=$2
			shift 2;; 	
		--prefix)
			prefix=$2 
			shift 2;; 
		--force) 
			FORCE=1
			shift;; 
		--log)
			log=$2
			shift 2;; 
		--learning_thr)
			learning_thr=$2 
			shift 2;; 
		--proteins)
			predicted_proteins=$2; 
			shift 2;; 
		-h) 
			usage 
			shift ;;
		--)  
			shift ; break ;; 					
	esac 
done	



BIN=/databis/hilpert/plasmidome_project/bin
BIN2=/databis/hilpert/plasmidome_realdata2/bin

source $BIN/common_functions.sh 

treat_args
verif_args

if [[ ! $predicted_proteins ]]; then 
	echo "PRELIMINARY STEP : PROTEINS PREDICTION" 
	dir=$outdir/protein_prediction
	mkdir -p $dir
	verif_result $dir/$prefix.predicted_proteins.faa 
	if [[ $file_exist == 1 ]]; then 
		echo "Predicted proteins file already exists. Use --force to overwrite."
	else 
		run_proteins_prediction $dir $prefix $assembly 
	fi
	predicted_proteins=$dir/$prefix.predicted_proteins.faa		
fi

echo "## STEP 1 : LEARNING" 
echo "[learning] Threshold : $learning_thr" 
dir=$outdir/learning
mkdir -p $dir 
thr=$(echo $learning_thr | awk '{print $1/100}')
verif_result $dir/$prefix.plasflow$thr.plasmids.fasta
if [[ $file_exist == 1 ]]; then
	echo "Learning results already exists. Use --force to overwrite" 
else 		
	echo "[learning] Run PlasFlow..." 
	bash $BIN/run_plasflow.sh $prefix $assembly -o $dir --thres $thr > $log/learning$learning_thr.log 
fi

echo "## STEP 2 : RNA SEARCH" 
dir=$outdir/rna_search 
mkdir -p $dir 
if [[ $FORCE ]]; then 
	bash $BIN2/pipeline_scripts/run_rna_search.sh -f $assembly -o $dir -d $rna_db --force
else 
	bash $BIN2/pipeline_scripts/run_rna_search.sh -f $assembly -o $dir -d $rna_db
fi	

echo "## STEP 3 : CHROMOSOMES SEARCH" 
dir=$outdir/chrm_search 
mkdir -p $dir 
if [[ $FORCE ]]; then 
	bash $BIN2/pipeline_scripts/eliminate_chrm.sh -f $assembly -o $dir -d $chrm_db --force
else 
	bash $BIN2/pipeline_scripts/eliminate_chrm.sh -f $assembly -o $dir -d $chrm_db
fi	

echo "## STEP 4 : PHYLOGENETIC MARKERS SEARCH"
dir=$outdir/phylogenetic_markers
mkdir -p $dir
if [[ $FORCE ]]; then 
	bash $BIN2/pipeline_scripts/search_phylogenetic_markers.sh -f $assembly -o $dir -d $phylo_db --proteins $predicted_proteins --force
else 
	bash $BIN2/pipeline_scripts/search_phylogenetic_markers.sh -f $assembly -o $dir -d $phylo_db --proteins $predicted_proteins
fi

echo "## STEP 5 : SEARCH PLASMIDS MARKERS" 
dir=$outdir/plasmids_markers 
mkdir -p $dir 
if [[ $FORCE ]]; then 
	bash $BIN2/pipeline_scripts/search_plasmids_markers.sh -f $assembly -o $dir -d $markers_db --force --proteins $predicted_proteins
else 
	bash $BIN2/pipeline_scripts/search_plasmids_markers.sh -f $assembly -o $dir -d $markers_db --proteins $predicted_proteins
fi	

echo "## STEP 6 : SEARCH CIRCULAR SEQUENCES" 
dir=$outdir/circular 
mkdir -p $dir 
out=$dir/$prefix.circular
verif_result $out.fasta 
if [[ $file_exists == 1 ]]; then 
	echo "Circular sequences results already exists. Use --force to overwrite" 
else 
	bash $BIN2/run_circular.sh $assembly $out.fasta
	grep "^>" $out.fasta | cut -f 1 -d " " | cut -f 2 -d ">" > $out.id 
fi	

