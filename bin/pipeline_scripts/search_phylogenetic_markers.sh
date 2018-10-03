set -e 

function usage(){ 
	echo "bash search_phylogenetic_markers.sh -f <assembly.fasta> -o <outdir> -d <Phylogenetic markers hmm profile> --proteins <predicted proteins fasta>
	Options : 
	--prefix <prefix> : prefix for results files, default : name of the fasta assembly
	--force : overwrite results if already exists"
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
		echo "You must give phylogenetic markers profile. Use -d option" 
	fi
	if [[ $quit ]]; then 
		exit 1
	fi 
}	

function verif_args(){
	verif_file $assembly "[phylo_markers_search] Assembly file doesn't found in $assembly" "[rna_search] Assembly file found in $assembly"
	verif_file $phylo_db "[phylo_markers_search] Phylogenetic markers profile doesn't found in $phylo_db" "[phylo_markers_search] Phylogenetic markers found in $phylo_db"
	mkdir -p $outdir 
	if [[ ! $prefix ]]; then 
		prefix=$(echo $assembly | rev | cut -f 1 -d "/" | cut -f 2- -d "." | rev)
	fi
}

function treat_hmm(){
	hmm=$1 
	grep -v "^#" $hmm.tsv | awk '{print $1}' | sort -u > $hmm.predicted_proteins.id 
	grep -v "^#" $hmm.tsv | awk '{print $3}' | sort -u > $hmm.proteins_family.id
	cat $hmm.predicted_proteins.id | rev | cut -f 2- -d "_" | rev > $hmm.contigs.id 
}	

TEMP=$(getopt -o h,f:,o:,d: -l prefix:,force,proteins:  -- "$@")
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
			phylo_db=$2
			shift 2;; 
		--prefix)
			prefix=$2 
			shift 2;; 
		--force) 
			FORCE=1
			shift;; 
		--proteins) 
			predicted_proteins=$2
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

echo "[phylo_search] Run Hmmer..." 
out=$outdir/$prefix.phylo_markers.hmm
verif_result $out.tsv 
if [[ $file_exist == 1 ]]; then 
	echo "Phylogenetic markers search results already exists. Use --force to overwrite" 
else 
	hmmsearch --tblout $out.tsv -E 1e-5 $phylo_db $predicted_proteins > $out
	treat_hmm $out
fi
