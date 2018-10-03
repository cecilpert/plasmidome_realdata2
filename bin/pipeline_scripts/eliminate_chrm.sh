set -e 

function usage(){ 
	echo "bash eliminate_chrm.sh -f <assembly.fasta> -o <outdir> -d <Chromosomes database.fasta>
	Options : 
	--prefix <prefix> : prefix for results files, default : name of the fasta assembly
	--force : overwrite results if already exists
	--log <directory> : directory to put log files. Default : outdir" 
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
	if [[ ! $chrm_db ]]; then 
		quit=1
		echo "You must give chromosomes database. Use -d option" 
	fi
	if [[ $quit ]]; then 
		exit 1
	fi 
}	

function verif_args(){
	verif_file $assembly "[chrm_search] Assembly file doesn't found in $assembly" "[chrm_search] Assembly file found in $assembly"
	verif_file $chrm_db  "[chrm_search] Chromosomes database doesn't found in $chrm_db" "[chrm_search] Chromosomes database found in $chrm_db"
	mkdir -p $outdir 
	if [[ ! $prefix ]]; then 
		prefix=$(echo $assembly | rev | cut -f 1 -d "/" | cut -f 2- -d "." | rev)
	fi

	if [[ ! $log ]]; then 
		log=$outdir/log	
	fi
	mkdir -p $log 	
	if [[ ! $cov ]]; then 
		cov=0.8
	fi	
}

TEMP=$(getopt -o h,f:,o:,d: -l prefix:,force,log:,cov:  -- "$@")
eval set -- "$TEMP" 
while true ; do 
	case "$1" in 
		-f)
			assembly=$2
			shift 2;; 
		-o)
			outdir=$2
			shift 2;; 
		-d) 
			chrm_db=$2
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
		--cov)
			cov=$2
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

echo "## CHROMOSOMES SEARCH"
echo "Coverage : $cov" 
out=$outdir/$prefix
verif_result $out.paf 
if [[ $file_exist == 1 ]]; then 
	echo "Chromosomes search results already exists. Use --force to overwrite" 
else 
	nb_chrm=$(grep "^>" -c $chrm_db) 
	echo "[chrm_search] Count chromosomes..." 
	echo "[chrm_search] Run minimap2..." 
	source activate plasmidome 
	minimap2 -x asm5 -N $nb_chrm $chrm_db $assembly 1> $out.paf 2> $outdir/log/align_chrm.log
	source deactivate plasmidome
	echo "[align_chrm] Treat minimap2..." 
	awk '{if ($10/$2 >= '$cov') print}' $out.paf > $out.$cov.paf 
	cut -f 1 $out.$cov.paf | sort -u > $out.$cov.chrm.id 
	echo "[align_chrm] Delete chromosomic contigs..."
	python3 $BIN/delete_seq_from_file.py $assembly $out.$cov.chrm.id $out.$cov.nochrm.fasta normal 
	grep "^>" $out.$cov.nochrm.fasta | cut -f 1 -d " " | cut -f 2 -d ">" > $out.$cov.nochrm.id
fi	


