#!/bin/bash

# Commands - A script to produce an AUTOMATED HMMSEARCH AND ANALYSIS TOOL

# Caputure pseudo-unique prefix
RIGHT_NOW=`date "+%y%m%d%H%M"`

WORKING_DIR=$(pwd) # pwd=print working directory
OUTPUT_DIR=$WORKING_DIR/Project_Results/"$RIGHT_NOW"

# Check if the results directory exists
if [ ! -d "$(pwd)/Project_Results" ]; then
	mkdir $(pwd)/Project_Results
fi

# Check if within the same minute another execution had been started.
# If so set the prefix with addition of the seconds.
if [ -d "$(pwd)/Project_Results/$RIGHT_NOW" ]; then
    PREFIX="$RIGHT_NOW"_`date +%s`
else
	PREFIX="$RIGHT_NOW"
fi

#==============================================================================
# Parsing Input
#==============================================================================
. ../AHATool/AHATool_Resources/shflags

# Define strings:
DEFINE_string prefix $PREFIX "The prefix the tool will use for produced files." p
DEFINE_string input 'sequences.fasta' "the input file (fasta, aln or hmm)." i
#DEFINE_string start 'build' "start of execution (search or build)." s
DEFINE_string output $OUTPUT_DIR "The path of the workspace" o
DEFINE_string database 'nr.fa' "database options: 1. nr_db; 2. custom_db" d
DEFINE_string update 'yes' "database update if possible? yes/no?" u
DEFINE_string cladogram 'yes' "Prepare tree file for cladogram? yes/no?" c	
# Define floats:
DEFINE_float evalue 0.0000000001 "e-value (recommended: 1e-10)." e
DEFINE_integer threads 2 "processor options: 1, 2, 4" t

# Other flags:
FLAGS "$@" || exit $?
eval set -- "${FLAGS_ARGV}"

#==============================================================================
# Logging
#==============================================================================
exec &> >(tee -a "$PREFIX"_log_file.txt)
LOGFILE="$PREFIX"_log_file.txt

#==============================================================================
# Preparatory: Checking nr.gz
#==============================================================================
# Checking if nr.gz is older than 30 minutes. If an update of a downloaded 
# 	nr.gz just happened there is no need to execute the code below.
# 1800 seconds amounts to 30 minutes (180000s = 50h).
if [ -f "nr.fa" ]; then
	if [ "${FLAGS_update}" = "yes" ] &&[ "${FLAGS_database}" = "nr.fa" ] && [ `stat --format=%Y nr.fa` -le $(( `date +%s` - 180000 )) ]; then 
		# Checking if an update of nr.gz is possible
		printf "Checking whether there is a new update of nr.gz available.\n"
		perl ../AHATool/AHATool_Resources/update_FASTAdb.pl nr;
		if [ $? -eq 0 ]; then
			echo "${`tput setaf 2`} No update was necessary. ${`tput sgr0`}"
		else
			echo "${`tput setaf 2`} Update done. ${`tput sgr0`}"
			exit
		fi
	fi
fi

#==============================================================================
# About/ Intro
#==============================================================================
##### Constants
TITLE="HMM Tool for $HOSTNAME"
TIME_STAMP="Executed on $(date +"%x %r") by $USER"
COUNTER=0
HITS_COUNT=0


PREFIX=${FLAGS_prefix}
DATABASE=${FLAGS_database}
THREADS=${FLAGS_threads}
INPUT_FILE=${FLAGS_input}
OUTPUT_DIR=${FLAGS_output}

START="build" #${FLAGS_start} #of execution (hmmbuild or hmmsearch)
E_VALUE=${FLAGS_evalue}
INPUT_DIR=${FLAGS_ARGV:1:-1}

# Replace link to non-existent input directory with link to working directory
if [ ! -d "$INPUT_DIR" ] 
	then INPUT_DIR=$WORKING_DIR
fi

if [ ! -f "$INPUT_FILE" ] 
	then echo $INPUT_FILE " is missing."
	exit
fi

# Make subdirectory in input directory
mkdir $OUTPUT_DIR

cp $INPUT_FILE $OUTPUT_DIR

# Colors:
RED=`tput setaf 1`
GREEN=`tput setaf 2`
ORANGE=`tput setaf 130`
ICE=`tput setaf 123`
RESET=`tput sgr0`

F0="$PREFIX"_"$INPUT_FILE" 
F1="$PREFIX"_temporaryfile1.txt 
F2="$PREFIX"_temporaryfile2.txt 
F3="$PREFIX"_temporaryfile3.txt 
F4="$PREFIX"_temporaryfile4.txt 
F5="$PREFIX"_temporaryfile5.txt 

##### Functions
check_for_packages ()
{
#$1 is the path to the file containing the list of software packages 
for s_ware in $(cat $1)
do
	if ! exists=$(apt-cache policy $s_ware)
		then echo "${RED} $s_ware is missing. ${RESET}"; COUNTER=$(($COUNTER+1))
	fi
done
if (($COUNTER==0)); then 
	echo "${GREEN} All needed packages are installed. ${RESET}"
fi
}   # end of check_for_packages

split_file ()
{
# Cut everything after column 18 so that we are left without 'description of target'
awk -v f=1 -v t=18 '{for(i=f;i<=t;i++) printf("%s%s",$i,(i==t)?"\n":OFS)}' $1 > $3
# Enter a tab after each empty space ' '
sed 's/  \?/\t/g' $3 > $2
# delete any double tab '\t\t' within a line
sed 's|\t\t*|\t|g' $2 > $3
# Stripping the "\tab" at the end of each line 
sed 's/\t$//' $3 > $2
# Cut first 18 columns so that we are left with 'description of target'
awk '{print substr($0, index($0,$19))}' $1 > $3
#also works:
#awk '{ s = ""; for (i = 19; i <= NF; i++) s = s $i " "; print s }' $1 > $3
}   # end of split_file

create_header ()
{
# Insert header
if [ "${FLAGS_database}" = "nr.fa" ] || [ "${FLAGS_database}" = "vr.fa" ]; then 
	sed -i '1s/^/target name\taccession\tquery name\taccession\te-value\tscore\tbias\te-value\tscore\tbias\texp\treg\tclu\tov\tenv\tdom\trep\tinc\tdescription of target\tsize (aa)\tSignalP prediction\tCS pos.\taccesssion number\tspecies\ttaxonomy\tsequence\n/' $1
else #custom
	sed -i '1s/^/target name\taccession\tquery name\taccession\te-value\tscore\tbias\te-value\tscore\tbias\texp\treg\tclu\tov\tenv\tdom\trep\tinc\tSignalP prediction\tCS pos.\tsize (aa)\tstart\tend\tdirection\n/' $1 #\tsize (aa)
fi
}   # end of create_header

check_success ()
{
if [ $? -eq 0 ]; then
	echo "${GREEN} Done. ${RESET}"
else
	echo "${RED} $1 failed. ${RESET}"
	exit
fi
}  # end of check_success

embed_signalP ()
{
printf "Starting signal peptide prediction...\n"
{
	signalp6 -fasta "$PREFIX"_intermediate_hits.fa -org gram- -format short -prefix "$PREFIX"_gram-_short -mature
	signalp6 -fasta "$PREFIX"_intermediate_hits.fa -org gram+ -format short -prefix "$PREFIX"_gram+_short -mature
	signalp6 -fasta "$PREFIX"_intermediate_hits.fa -org arch -format short -prefix "$PREFIX"_arch_short -mature
} 1>/dev/null 2>&1

# Getting the "Prediction":
sed -e '1,2d' < "$PREFIX"_gram-_short_summary.signalp6 > $3
cut -f2 $3 > $2
cut -c-15 $2  > $1
# Stripping the "\tab" at the end of each line 
sed 's/\t$//' $1> $2 

# Getting the "CS Position":
cut -f7 $3 > $1
sed 's/[.].*$//' $1> $3
sed -r 's/.*[:] //' $3 > $1

#patch'em
paste $2 $1 > $3
check_success "Output files"

find $WORKING_DIR/ -maxdepth 1 -name "$PREFIX"_gram-_short_mature.fasta -exec mv {} $OUTPUT_DIR \;
find $WORKING_DIR/ -maxdepth 1 -name "$PREFIX"_gram+_short_mature.fasta -exec mv {} $OUTPUT_DIR \;
find $WORKING_DIR/ -maxdepth 1 -name "$PREFIX"_arch_short_mature.fasta -exec mv {} $OUTPUT_DIR \;
find $WORKING_DIR/ -maxdepth 1 -name "$PREFIX"_gram-_short_summary.signalp6 -exec mv {} $OUTPUT_DIR \;
find $WORKING_DIR/ -maxdepth 1 -name "$PREFIX"_gram+_short_summary.signalp6 -exec mv {} $OUTPUT_DIR \;
find $WORKING_DIR/ -maxdepth 1 -name "$PREFIX"_arch_short_summary.signalp6 -exec mv {} $OUTPUT_DIR \;
}  # end of embed_signalP

reformat_fasta ()
{
# Paste sequence at the end of SummaryTable.txt
sed 's/>.*/>/' $1 > $2
echo ">" >> $2
# Concatenate lines until appearance of '>'
sed -e '
    1d
    :loop
	$!N
	s/\n>//
    Tloop
    y/\n/,/
' $2 > $3
}  # end of reformat_fasta

efetch_file ()
{
for i in `cat $1`; do 
    printf ${i}"\t"; \
	efetch -db protein -id ${i} -format xml \
    | xtract -pattern Seq-entry -element Org-ref_taxname, OrgName_lineage, PubMedId ; \ #, PubMedId
done > $2

# ...nucleotide sequence
for i in `cat $1`; do 
    printf ${i}"\t"; \
	efetch -db protein -format fasta_cds_na -id ${i};\
    printf "\n"; \
done > $PREFIX"_coding_sequence.fasta" #$PREFIX"_protein_nuccore_mrna.fasta"

cat "SubSource_name\t \t Prot-ref_name_E \n" > $PREFIX"_additional_information.txt"
for i in `cat $1`; do 
    printf ${i}"\t"; \
#extract even more? - does not work quiet yet...
	efetch -db protein -id ${i} -format xml \
	    | xtract -pattern Seq-entry -element SubSource_name, Prot-ref_name_E; \ #PubMedId,  
	#esearch -db protein -query -id ${i} | elink -target nucleotide -name protein_nuccore | efetch -format xml \
	#| xtract -pattern Seq-entry -element PubMedId, Prot-ref_name_E, SubSource_name; \
    printf "\n"; \
done > $PREFIX"_additional_information.txt" #$PREFIX"_protein_nuccore_mrna.fasta"
} # end of efetch_file


check_for_taxonomy ()
{
# Format for PDB IDs:
#reg_ex_0=[0-9][a-zA-Z0-9_]{3} 
reg_ex_1=pdb_ID
# Format for GenBank and RefSeq accession numbers:
reg_ex_2=[a-zA-Z]{2}[a-zA-Z_][0-9]{5}
reg_ex_3=[a-zA-Z]{2}[a-zA-Z_][0-9]{7}
reg_ex_4=[a-zA-Z]{2}[a-zA-Z_][0-9]{9}
# Format for UniProt accession numbers:
reg_ex_5=[OPQ][0-9][A-Z0-9]{3}[0-9]
reg_ex_6=[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]

#LOOP THAT ITERATES OVER THE POSSIBLE COMBINATIONS:
if (($(cat $2 | wc -l)==0)); then
	if [ "${FLAGS_database}" = "nr.fa" ]; then 
		echo "${RED} Gathering taxonomy information failed within the $2 quarter of the data. ${RESET}"
	else
		echo "${RED} Gathering taxonomy information failed for some of the sequences. ${RESET}"
	fi	
	sed -i -e 's/$/   \t - \t - /' "$1" > "$2"
	sed -i -e 's/^[0-9][a-zA-Z0-9_][a-zA-Z0-9_][a-zA-Z0-9_]$/pdb_ID/' $2
elif (($(cat $1 | wc -l)!=$(cat $2 | wc -l))); then
	printf "${ORANGE}An error occurred while gathering some tax. information within the "$2" 
quarter of the data. Please check for missing taxa in the Summary Table.\n${RESET}"
	sed -i -e 's/pdb_ID\t/pdb_ID   \t - \t -\n /' $2
	for  var in `seq 1 6`; 
	do
		if (($var==1)); then i=$reg_ex_1; elif (($var==2)); then i=$reg_ex_2; elif (($var==3)); then i=$reg_ex_3; elif (($var==4)); then i=$reg_ex_4; elif (($var==5)); then i=$reg_ex_5; elif (($var==6)); then i=$reg_ex_6; else i=$reg_ex_7; fi
		for  var_2 in `seq 1 6`;
		do
			if (($var_2==1)); then j=$reg_ex_1; elif (($var_2==2)); then j=$reg_ex_2; elif (($var_2==3)); then j=$reg_ex_3; elif (($var_2==4)); then j=$reg_ex_4; elif (($var_2==5)); then j=$reg_ex_5; elif (($var_2==6)); then j=$reg_ex_6; else j=$reg_ex_7; fi

			number=$(awk '/'$i'\t'$j'/{ print NR; exit }' $2)
			while [ $number>0 ] ; do
				sed -i -e "$number s/\t/   \t -\t -\n/" $2
				number=$(awk '/'$i'\t'$j'/{ print NR; exit }' $2)

			done
			number=$(awk '/'$j'\t'$i'/{ print NR; exit }' $2)
			while [ $number>0 ] ; do
				sed -i -e "$number s/\t/   \t -\t -\n/" $2
				number=$(awk '/'$j'\t'$i'/{ print NR; exit }' $2)

			done
			number=$(awk '/'$i'\n/{ print NR; exit }' $2)
			while [ $number>0 ] ; do
				sed -i -e "$number s/\n/   \t -\t -\n/" $2
				number=$(awk '/'$i'\n/{ print NR; exit }' $2)

			done
			number=$(awk '/'$i'\t\n/{ print NR; exit }' $2)
			while [ $number>0 ] ; do
				sed -i -e "$number s/\t\n/   \t -\t -\n/" $2
				number=$(awk '/'$i'\n/{ print NR; exit }' $2)

			done
		done
	done
	number=$(awk 'length<13{ print NR; exit }' $2)
	if [ $number>0 ] ; then
		sed -i -e "$number s/\t/   \t - \t - \n/" $2
	fi
fi
while (($(cat $1 | wc -l)>$(cat $2 | wc -l))); do
	printf "\t -\t -\n" >> "$2"
done
echo "Done with: "$2" quarter. "
} # end of check_for_taxonomy

format_time () {
  ((h=$1/3600))
  ((m=($1%3600)/60))
  ((s=$1%60))
  printf "%02d:%02d:%02d\n" $h $m $s
 }  # end of format_time

check_online_connection () {
wget -q --spider https://www.ncbi.nlm.nih.gov
if [ $? -eq 0 ]; then
    echo "Checked for online connection to NCBI."
else
    echo "${RED} As no connection to NCBI could be established the execution of AHATool is stopped.${RESET}"
	echo "This could be caused by a lack of online connectivity. Please check network options in order to establish an online connection."
	exit
fi
 }  # end of check_online_connection

##### Main
printf "
=======================================================================
Welcome to AHATool: an Automatic HMM and Analysis Tool.
V.2 [28.Jul 2021]
Microbiology and Biotechnology - Streit's lab
University of Hamburg (D)
Developed by Nele Schulte (and P.Pérez-García)
=======================================================================
$TIME_STAMP
AHATool will be executed with the following parameters: 
Directory to work in: $INPUT_DIR
Database: $DATABASE
Number of threads: $THREADS
Input file: $INPUT_FILE
E-value: $E_VALUE
Output prefix: $PREFIX
RAxML tree: ${FLAGS_cladogram}
"
# Start of execution: $START profile HMM #outdated

#==============================================================================
# Initial Checking
#==============================================================================
check_online_connection
# software_needed=all software that should be in installed
printf "=======================================================================
Checking for needed software:
---------------------------------\n"
software_needed=$WORKING_DIR/../AHATool/AHATool_Resources/SOFTWAREneeded.txt
check_for_packages $software_needed

if (($COUNTER>0)); then 
	echo "${RED}The tool can not proceed due to the lack of $COUNTER software packages (listed above).
Please approve installation. Do the installation via 'sudo apt-get install [name of software package]'. ${RESET}"
	exit
fi

# Resetting the count
COUNTER=0
SIGNALP_V=6.0
# Check for SignalP
if ! signalp6 --version | grep -q 'SignalP 6.0'; then
	echo "${RED} Check for SignalP failed. The tool can not proceed due to the lack of SignalP.
Please approve installation. Do the installation as described on the official website: http://www.cbs.dtu.dk/services/SignalP/portable.php. ${RESET}"
	exit
fi

# Check for epost (needed for eulilities like epost)
if ! epost -version | grep -q '[0-9.]'; then
	echo "${RED} Check for epost failed. The tool can not proceed due to the lack of EDirect.
Please approve installation. Do the installation as described on the official website: https://www.ncbi.nlm.nih.gov/books/NBK179288/. ${RESET}"
	exit
fi

# Check for blastp
if ! blastp -version | grep -q 'blastp:'; then
	echo "${RED} Check for blastp failed. The tool can not proceed due to the lack of ncbi-blast.
Please approve installation. Do the installation as described on the official website: https://www.ncbi.nlm.nih.gov/books/NBK279690/. ${RESET}"
	exit
fi

printf "=======================================================================
Checking the input folder...
---------------------------------\n"
if [ ! -d "$INPUT_DIR" ]; then
	echo "$INPUT_DIR"
    echo "${RED}The input directory was not found."
    echo "Kindly check and input the correct file path. ${RESET}"
    exit
fi

# Use wc to do a count of the number of lines (-l) in the output of ls -1.
printf "Files found: "
ls -1 | wc -l
#looking for relevant kinds of FASTA files: .fasta, .faa, .fa:
printf "FASTA files found: "
COUNTER=$(($COUNTER + $(find $INPUT_DIR -maxdepth 1 -name "*.fasta" | wc -l)))
COUNTER=$(($COUNTER + $(find $INPUT_DIR -maxdepth 1 -name "*.faa" | wc -l)))
COUNTER=$(($COUNTER + $(find $INPUT_DIR -maxdepth 1 -name "*.fa" | wc -l)))
printf "$COUNTER \n"

#If only one file with matching format is found use that file
if (($COUNTER==1)); then 
	INPUT_FILE=$(find $INPUT_DIR -type f \( -iname \*.fasta -o -iname \*.faa -o -iname \*.fa -o -iname \*.hmm \))
fi

if (($COUNTER==0)); then 
	echo "${RED}The tool can not proceed due to the lack of an INPUT files 
(any input file with the format .hmm, 1	a.fasta, .faa or .fa).${RESET}"
	exit
fi

printf "======================================================================= 
Checking for needed files:
---------------------------------\n"
# Check in the input folder (given as parameter) for database

if [ -e "$INPUT_DIR/$DATABASE" ]
			then echo "${GREEN} Database ($DATABASE) exists. ${RESET}" # in $2"
			else echo "${RED} Database ($DATABASE) is missing in $INPUT_DIR. ${RESET}" 
				echo "${RED} The tool can not proceed due to the lack of a database. ${RESET}"
				exit
fi

# Checking for existence of an index file (namely $DATABASE.ssi)
if [ -e "$INPUT_DIR/$DATABASE.ssi" ]
then 
	echo "${GREEN} Index file ($DATABASE.ssi) exists. ${RESET}"
else 
	echo "Index file ($DATABASE.ssi) is missing. An index file is being generated..."
	# Creating an index file; First thing: create an "SSI index" for that file:
	esl-sfetch --index $DATABASE
	check_success "Creation of index file"
fi

printf "=======================================================================\n"
printf "${ICE}The files created within this run will be identifiable by the prefix 
\"$PREFIX\".\n${RESET}"
printf "${ICE}They will be saved in the folder /Project_Results 
and the subfolder \"$RIGHT_NOW\".\n${RESET}"

if [[ $INPUT_FILE == *.hmm ]]; then
	START="search"
fi
#==============================================================================
# Multiple Sequence Alignment (MSA)
#==============================================================================

if [ "$START" == "build" ] || (($COUNTER==0)); then 

#checking for file *".aln" in the input directory
COUNTER=$(($COUNTER + $(find $INPUT_DIR -maxdepth 1 -name "*.aln" | wc -l)))

# This file decides where the tool starts (hmm -build oder -search)
if (($COUNTER==0)); then
    echo "No matching .aln file found."
else 
    echo "${GREEN} A file with an .aln format was found. ${RESET}"
fi
	printf "=======================================================================
Initiating MSA using t-coffee...\n"
	{
	t_coffee "$INPUT_FILE" -run_name "$F0".aln
	} 1>/dev/null 2>&1
	# -mode accurate
		
	printf "Reformating Clustal to FASTA...\n"
	# Clustal to Fasta reformatting # -in "$F0"
	t_coffee -other_pg seq_reformat -in "$F0".aln -output fasta_aln > "$F0".fasta
	#t_coffee -other_pg seq_reformat -in sequences.aln -output fasta_aln > "$PREFIX"_"$INPUT_FILE"_"aln.fasta" 
	check_success "MSA"
#else 	
	##touch $WORKING_DIR/"$F0".aln
	#cat $(find $INPUT_DIR -maxdepth 1 -name "*"$INPUT_FILE".aln") >"$F0".aln
#==============================================================================
# Build an HMM
#==============================================================================
# Resetting the count
COUNTER=0 

#checking for file *".hmm" in the input directory
COUNTER=$(($COUNTER + $(find $INPUT_DIR -maxdepth 1 -name "*.hmm" | wc -l)))

# This file decides where the tool starts (hmm -build oder -search)
if (($COUNTER==0)); then
    echo "$No .hmm file found. AHATool will start with building a HMM profile."
else 
    echo "${GREEN} A file with an .hmm format was found. ${RESET}"
	if [ "$START" = "search" ]; then 
		echo "Therefore tool will start with a search right away."
	fi
fi

#if [ "$START" = "build" ] || (($COUNTER==0)); then 
	printf "=======================================================================
Constructing HMM profile from MSA...\n"
	start=`date +%s`
	{
	hmmbuild "$PREFIX"_"$INPUT_FILE".hmm "$F0".fasta
	} 1>/dev/null 2>&1
	end=`date +%s`
	printf "${ICE}... completed in $(format_time `expr $end - $start`)\n${RESET}"
	check_success "HMM profile"
#else 	
#	cat $(find $INPUT_DIR -maxdepth 1 -name "*"$INPUT_FILE".hmm") > "$PREFIX"_"$INPUT_FILE".hmm
#	cat $(find $INPUT_DIR -maxdepth 1 -name "*"$INPUT_FILE".html") > "$PREFIX"_"$INPUT_FILE".html
#	cat $(find $INPUT_DIR -maxdepth 1 -name "*"$INPUT_FILE".fasta") > "$PREFIX"_"$INPUT_FILE".fasta
#	cat $(find $INPUT_DIR -maxdepth 1 -name "*"$INPUT_FILE".dnd") > "$PREFIX"_"$INPUT_FILE".dnd
elif  [ "$START" = "search" ]; then 
	cp $INPUT_FILE "$PREFIX"_"$INPUT_FILE"
fi
# HMM search: searching the constructed protein profile HMM against a protein sequence database.
printf "=======================================================================
Searching profile HMM against given database...\n"
if [ "$START" = "search" ]; then 
	{
	hmmsearch -o "$PREFIX"_"$INPUT_FILE"_hmm.out -A "$PREFIX"_"$INPUT_FILE"_hmm.aln --tblout "$PREFIX"_"$INPUT_FILE"_hmm.tbl --noali -E $E_VALUE --cpu $THREADS "$PREFIX"_"$INPUT_FILE" $DATABASE
	} 1>/dev/null 2>&1
else
	{
	hmmsearch -o "$PREFIX"_"$INPUT_FILE"_hmm.out -A "$PREFIX"_"$INPUT_FILE"_hmm.aln --tblout "$PREFIX"_"$INPUT_FILE"_hmm.tbl --noali -E $E_VALUE --cpu $THREADS "$PREFIX"_"$INPUT_FILE".hmm $DATABASE
	} 1>/dev/null 2>&1
fi
check_success "Profile HMM search"

# Counting the hits:
HITS_COUNT=$(grep "^[^#;]" "$PREFIX"_"$INPUT_FILE"_hmm.tbl |wc -l)
printf "${GREEN} $HITS_COUNT sequences found with e-value ≤ $E_VALUE ${RESET}"

# Moving files that are not needed for further use
if [ -f "$F0".fasta ]; then
	find $WORKING_DIR/ -maxdepth 1 -name "$F0".fasta -exec mv {} $OUTPUT_DIR \;
fi
if [ -f  "$F0".aln ]; then
	find $WORKING_DIR/ -maxdepth 1 -name "$F0".aln -exec mv {} $OUTPUT_DIR \;
fi
if ! [ "$START" = "search" ]; then
	find $WORKING_DIR/ -maxdepth 1 -name "$PREFIX"_"$INPUT_FILE".hmm -exec mv {} $OUTPUT_DIR \;
else
	find $WORKING_DIR/ -maxdepth 1 -name "$PREFIX"_"$INPUT_FILE" -exec mv {} $OUTPUT_DIR \;
fi
find $WORKING_DIR/ -maxdepth 1 -name "$PREFIX"_"$INPUT_FILE".html -exec mv {} $OUTPUT_DIR \;
find $WORKING_DIR/ -maxdepth 1 -name "$PREFIX"_"$INPUT_FILE".fasta -exec mv {} $OUTPUT_DIR \;
find $WORKING_DIR/ -maxdepth 1 -name "$PREFIX"_"$INPUT_FILE".dnd -exec mv {} $OUTPUT_DIR \;
find $WORKING_DIR/ -maxdepth 1 -name "$PREFIX"_"$INPUT_FILE"_hmm.out -exec mv {} $OUTPUT_DIR \;
find $WORKING_DIR/ -maxdepth 1 -name "$PREFIX"_"$INPUT_FILE"_hmm.aln -exec mv {} $OUTPUT_DIR \;

#==============================================================================
# Summary Table
#==============================================================================
printf "\n=======================================================================
Creating a summary table...\n"

# Extract 1st column with IDs
grep "^[^#;]" "$PREFIX"_"$INPUT_FILE"_hmm.tbl | cut -d" " -f1 > "$F1"
{
esl-sfetch -f $DATABASE "$F1" > "$PREFIX"_all_hits.fa
} 1>/dev/null 2>&1
check_success "Creation of hit list"
sed 's/\].*/].../' "$PREFIX"_all_hits.fa > "$PREFIX"_intermediate_hits.fa

#==============================================================================
# Taxonomy
#==============================================================================
if [ "${FLAGS_database}" = "nr.fa" ]; then 

	printf "======================================================================= 
Gathering taxonomy information...\n"
	start=`date +%s`

	# Stripping the ".1" at the end of each line in and writing the output in a new file.
	# This way we are left with only the accession number.
	sed 's/[[:punct:]].$//' "$F1">"$F2"

	# Resetting the count
	COUNTER=0 

	# 1. get number of lines
	COUNTER=$(cat $F2 | wc -l)
	
	# Create file to make sure they can be read even when NCBI does not respond
	touch taxonomy_aa
	touch taxonomy_ab
	touch taxonomy_ac
	touch taxonomy_ad

	# Format PDB IDs:
	sed -i -e 's/^[0-9][a-zA-Z0-9_][a-zA-Z0-9_][a-zA-Z0-9_]$/pdb_ID/' $F2

	# 2. split the file to work with into four equal sized files
	lines_count=$((($COUNTER+3)/4)) #$COUNTER+3 in order to round up
	split -l "$lines_count" $F2 taxonomy_

	if [[ "$THREADS" -eq 4 ]]; then

		# 3. run efetch for all four files simultaneously (4 threads)
		{
		efetch_file taxonomy_aa first &
		efetch_file taxonomy_ab second &
		efetch_file taxonomy_ac third &
		efetch_file taxonomy_ad fourth
		wait
		} 1>/dev/null 2>&1

		# 4. Check for success of gethering taxonomy information
		check_for_taxonomy taxonomy_aa first 
		check_for_taxonomy taxonomy_ab second 
		check_for_taxonomy taxonomy_ac third 
		check_for_taxonomy taxonomy_ad fourth

	elif [[ $THREADS -eq 2 ]]; then

		# 3. run efetch for two files simultaneously (2 threads)
		{
		efetch_file taxonomy_aa first &
		efetch_file taxonomy_ab second
		efetch_file taxonomy_ac third &
		efetch_file taxonomy_ad fourth
		wait
		} 1>/dev/null 2>&1

		# 4. Check for success of gethering taxonomy information
		check_for_taxonomy taxonomy_aa first 
		check_for_taxonomy taxonomy_ab second
		check_for_taxonomy taxonomy_ac third 
		check_for_taxonomy taxonomy_ad fourth

	else #if (($THREADS=1)); then

		# 3. run efetch for all four files seperately (1 thread)
		{
		efetch_file taxonomy_aa first
		efetch_file taxonomy_ab second
		efetch_file taxonomy_ac third
		efetch_file taxonomy_ad fourth
		wait
		} 1>/dev/null 2>&1

		# 4. Check for success of gethering taxonomy information
		check_for_taxonomy taxonomy_aa first
		check_for_taxonomy taxonomy_ab second
		check_for_taxonomy taxonomy_ac third
		check_for_taxonomy taxonomy_ad fourth

	fi
	rm taxonomy_aa
	rm taxonomy_ab
	rm taxonomy_ac
	rm taxonomy_ad

	# 5. patch taxonomy_1 though taxonomy_4 back together
	sed -i -e '$a\' first
	cat second >> first 
	sed -i -e '$a\' first
	cat third >> first 
	sed -i -e '$a\' first
	cat fourth >> first 
	wait

	end=`date +%s`
	printf "${ICE}... completed in $(format_time `expr $end - $start`)\n${RESET}"
	rm second
	rm third
	rm fourth
	printf "Adding lineage information to the summary table...\n"

	# Remove the lines containing the string '#'
	sed '/#/d' "$PREFIX"_"$INPUT_FILE"_hmm.tbl>"$F3"
	find $WORKING_DIR/ -maxdepth 1 -name "$PREFIX"_"$INPUT_FILE"_hmm.tbl -exec mv {} $OUTPUT_DIR \;
	split_file "$F3" "$F1" "$F2" #split_file_144

	# Paste description after tbl
	paste "$F1" "$F2" > "$F3"
	reformat_fasta "$PREFIX"_intermediate_hits.fa "$F4" "$F1"

	# Remove all ","
	sed 's/,//g' "$F1" > "$F4"

	# Cut so that we are left without 'PubMedId'
	sed -i -e 's/\t[0-9]*\t[0-9]*$//' first

	paste first "$F4" > "$F2"
	rm first

	# Calculate size(aa) 
	awk '{ print length }' "$F4" > "$F1"

	# Paste size(aa) after description
	paste "$F3" "$F1" > "$F4"
	rm -f "$F3"

	# Run SignalP
	start=`date +%s`
	embed_signalP "$F5" "$F1" "$PREFIX"_temporaryfile
	end=`date +%s`
	printf "${ICE}... completed in $(format_time `expr $end - $start`)\n${RESET}"
	rm "$PREFIX"_intermediate_hits.fa

	# Paste SignalP at the beginnig of intermediate file
	paste "$PREFIX"_temporaryfile "$F2" > "$F5"
	sed -i '1s/^/ \n/' "$F5"

	# Insert header
	create_header "$F4"

	# Final paste
	paste "$F4" "$F5" > "$PREFIX"_Summary.tsv 
	check_success "Creation of summary table"
	cp "$PREFIX"_Summary.tsv "$PREFIX"_Summary.xls

else
	split_file "$PREFIX"_"$INPUT_FILE"_hmm.tbl "$F4" "$F1" 
	find $WORKING_DIR/ -maxdepth 1 -name "$PREFIX"_"$INPUT_FILE"_hmm.tbl -exec mv {} $OUTPUT_DIR \; 
	sed '/#/d' "$F4">"$F3" 

	# Enter a tab instead of each ';'
	sed 's/[;]/ /g' "$F1" > "$F2"
	sed 's/[#]/ /g' "$F2" > "$F1"

	# Save copy of first line for header
	head -4 "$F1" > "$PREFIX"_temporaryfile1
	sed 's/[=]/= /g' "$F1" > "$F2"

	# Delete all occurences of any word containing '='
	sed 's/[a-zA-Z0-9_!]*[=][a-zA-Z0-9_!]*//g' "$F2" > "$F1"

	# Stripping the "\tab" at the end of each line 
	sed 's/\t$//' "$F1"> "$F2" 

	# Removing the first three lines
	sed -e '1,3d' < "$F2" > "$F1"

	# Enter a tab after each empty space ' '
	sed 's/  \?/\t/g' "$F1" > "$F2"
	sed 's|\t\t*|\t|g' "$F2" >"$F1"

	# Removing the first three lines
	sed -e '1,3d' < "$PREFIX"_temporaryfile1 > "$F2"
	sed 's/[=]/ =/g' "$F2" > "$PREFIX"_temporaryfile1 

	# Delete all occurences of any word containing '='
	sed 's/[a-zA-Z0-9_!]*[=][a-zA-Z0-9_!]*//g' "$PREFIX"_temporaryfile1 > "$F2"
	sed 's/[.0-9-]*[.0-9-][a-zA-Z0-9-]*//g' "$F2" > "$PREFIX"_temporaryfile1 

	# Enter a tab after each empty space ' '
	sed 's/  \?/\t/g' "$PREFIX"_temporaryfile1 > "$F2"
	sed 's|\t\t*|\t|g' "$F2"> "$PREFIX"_temporaryfile1 

	# Write header for "$F1"
	sed -i '1s/^/\n/' "$F1"
	paste "$F1" "$PREFIX"_temporaryfile1 > "$F2"
	sed 's|\t\t*|\t|g' "$F2"> "$PREFIX"_temporaryfile1

	# Delete all lines file starting from after a matching line
	sed -n '/hmmsearch/q;p' "$PREFIX"_temporaryfile1 > "$PREFIX"_temporaryfile

	# Sum up start and end of sequence to get size(aa)
	awk '{$3=$2-$1;} {print $3}' "$PREFIX"_temporaryfile > "$PREFIX"_temporaryfile1 
	paste "$PREFIX"_temporaryfile1 "$PREFIX"_temporaryfile > "$F2"
	sed 's|\t\t*|\t|g' "$F2"> "$PREFIX"_temporaryfile1

	start=`date +%s`
	embed_signalP "$F5" "$F1" "$F4"
	end=`date +%s`
	printf "${ICE}... completed in $(format_time `expr $end - $start`)\n${RESET}"

	# Paste SignalP at the beginnig of SummaryTable.txt
	paste "$F3" "$F4" >"$F5"

	# Insert header
	create_header "$F5"
	
	# Stripping the "\tab" at the end of each line 
	sed 's/\t$//' "$PREFIX"_temporaryfile1 > "$F2"

	# Delete first and last line
	sed 's/^0\tID\t/ID\t/' "$F2"> "$PREFIX"_temporaryfile1 
	head -n -1 "$PREFIX"_temporaryfile1 > "$F2"

	# Delete leading tab of each line
	sed 's/^\t//' "$F2"> "$PREFIX"_temporaryfile1  

	# Delete leading whitespace of each line
	sed 's/^ *//' "$PREFIX"_temporaryfile1 > "$F2" 

	# Pacht'em
	paste "$F5" "$F2" > "$PREFIX"_Summary.tsv 
	reformat_fasta "$PREFIX"_intermediate_hits.fa "$F4" "$F1" 

	# Remove all ","
	sed 's/,//g' "$F1" > "$F4"
	echo "now"
	sed -i '1s/^/sequence \n/' "$F4"

	paste  "$PREFIX"_Summary.tsv "$F4" > "$F1"

	printf "blasting hits against non-redunduant database...\n" 
	#printf "${ORANGE} For the $HITS_COUNT sequences found this will take approx. between
	# 15 and 60 minutes depending on the time of day (exhaustion of the BLASTing servers).\n${RESET}"
	printf "${ORANGE} A BLAST run takes on average between 15 and 25 minutes depending
 on the amount of sequences to BLAST and the time of day (exhaustion of the BLASTing servers).\n${RESET}"
	start=`date +%s`
	if [[ "$THREADS" -eq 4 ]] || [[ $THREADS -eq 2 ]]; then

		# run on 2 threads
		{
		blastp -task blastp-fast -db nr -query "$PREFIX"_all_hits.fa -out "$PREFIX"_summary_blast.txt -remote -max_target_seqs 1 -outfmt '10 sacc qacc evalue bitscore qcovs pident stitle staxids' &
		blastp -task blastp-fast -db nr -task blastp-fast -query "$PREFIX"_all_hits.fa -out "$PREFIX"_"$INPUT_FILE"_blast.out -remote -max_target_seqs 1
		wait
		} 1>/dev/null 2>&1
		# Check for success 
		check_success "Blasting sequences"

	else #if (($THREADS=1)); then

		# run seperately
		{
		blastp -task blastp-fast -db nr -query "$PREFIX"_all_hits.fa -out "$PREFIX"_summary_blast.txt -remote -max_target_seqs 1 -outfmt '10 sacc qacc evalue bitscore qcovs pident stitle staxids'
		blastp -task blastp-fast -db nr -query "$PREFIX"_all_hits.fa -out "$PREFIX"_"$INPUT_FILE"_blast.out -remote -max_target_seqs 1
		wait
		} 1>/dev/null 2>&1
		# Check for success 
		check_success "Blasting sequences"
 
	fi
	end=`date +%s`
	printf "${ICE}... completed in $(format_time `expr $end - $start`)\n${RESET}"

	find $pwd -maxdepth 1 -name "$PREFIX"_"$INPUT_FILE"_blast.out -exec mv {} $OUTPUT_DIR \;

	printf "Gathering taxonomy information...\n"
	start=`date +%s`
	# Delete everything after ","
	sed 's/,.*/ /' "$PREFIX"_summary_blast.txt > "$F4"

	# run efetch
	{
	efetch_file "$F4" "$F3"
	wait
	} 1>/dev/null 2>&1
	check_success "Gathering taxonomy" 
	# Check for success of gethering taxonomy information
	check_for_taxonomy "$F4" "$F3"
	end=`date +%s`
	printf "${ICE}... completed in $(format_time `expr $end - $start`)\n${RESET}"

	sed -i '1s/^/Acc. No.\tSpecies\tTaxonomy\tPubMedID\n/' "$F3"

	# Replace all ","
	sed 's/,/ \t/g' "$PREFIX"_summary_blast.txt > "$F4"
	sed -i '1s/^/Subject accession\tQuery accesion\tevalue\tbitscore\tQuery coverage\tPercentage of identical matches\tTitle\tstaxids\n/' "$F4"

	paste  "$F1" "$F4" > "$F2" 
	paste  "$F2" "$F3" > "$PREFIX"_Summary.tsv
	check_success "Creation of summary table"

printf "Creating taxonomy plot...\n"
#sed -i -e 's/^([^:]*.[^:]*):.*$/\1/' "$F4"

check_success "Creation of taxonomy plot"

	rm -f "$F3"
	rm -f "$F2"

#auskommentiert für Debugging:
#	rm -f "$PREFIX"_summary_blast.txt
	cp "$PREFIX"_Summary.tsv "$PREFIX"_Summary.xls
	find $pwd -maxdepth 1 -name "$PREFIX""_Summary*" -exec mv {} $OUTPUT_DIR \;

#CONSTRUCTIONSITE :D 13.07.2021
# both don't work as I need them to but kind of half way
#paste -d ";" - - <"$PREFIX"_intermediate_hits.fa | awk 'BEGIN{OFS=FS=";"}{print $1,$2,$3,$4"\n"$5>$2".fa"}'
#paste -d ";" - - <"$PREFIX"_intermediate_hits.fa | awk 'BEGIN{OFS=FS=";"}{if($2=="partial=00"){print $1,$2,$3,$4"\n"$5>"non_partial.fa"}else{print $1,$2,$3,$4"\n"$5>"partial.fna"}}'

	rm "$PREFIX"_intermediate_hits.fa
fi

printf "${ICE}The summary table contains the signal peptide predictions of a
 SignalP 'gram-' analysis.\n${RESET}"

#==============================================================================
# Cladogram
#==============================================================================
if [ "${FLAGS_cladogram}" = "yes" ]; then 
	if  [[ 100 -gt "$HITS_COUNT" ]]; then 
		#HITS_COUNT<100
		# Integration of RAxML
		printf "=======================================================================
Preparing files for cladogram...\n"
		start=`date +%s`
		{
		t_coffee "$PREFIX"_all_hits.fa -run_name "$PREFIX"_hits.aln
		} 1>/dev/null 2>&1
		# -mode accurate
		
		printf "Reformating Clustal to FASTA...\n"
		# Clustal to Fasta reformatting # -in "$F0"
		t_coffee -other_pg seq_reformat -in "$PREFIX"_hits.aln -output fasta_aln > "$PREFIX"_hits.fasta
		#t_coffee -other_pg seq_reformat -in sequences.aln -output fasta_aln > "$PREFIX"_"$INPUT_FILE"_"aln.fasta" 
		check_success "MSA (for cladogram)"
		{
		raxmlHPC -T 2 -f a -x 445 -p $PREFIX -s "$PREFIX"_hits.fasta -m PROTGAMMAWAG -n "$PREFIX"_hits -k -N 100
		} 1>/dev/null 2>&1
		check_success "RAxML"
		end=`date +%s`

		rm RAxML_info."$PREFIX"_hits
		rm RAxML_bipartitionsBranchLabels."$PREFIX"_hits
		rm RAxML_bestTree."$PREFIX"_hits
		rm RAxML_bootstrap."$PREFIX"_hits

		find $pwd -maxdepth 1 -name "*$PREFIX" mv {} $OUTPUT_DIR \;

		printf "${ICE}... completed in $(format_time `expr $end - $start`)\n${RESET}"
		printf "${ORANGE} The created tree files in netwick format will have to
 be visualized with MEGAX or a comparable program.\n${RESET}"
	else
		printf "${RED} As the number of hits for this project exceeds 100 
sequences ($HITS_COUNT) no files for a cladogram are being created.\n${RESET}"
	fi
fi
#==============================================================================
# Finishing the job stuff
#==============================================================================
printf "=======================================================================
Finishing...\n"

printf "${ICE}Script completed in $(format_time $SECONDS)\n${RESET}"
#}
# Changes in "$PREFIX"_log_file.txt: Delete all occurences of the color codes
sed -i -e 's/\x1B[^m]*m//g' "$LOGFILE" #"$PREFIX"_log_file.txt > "$F1" #$"$USER"_log_file.txt
#find $pwd -maxdepth 1 -name "$PREFIX"_log_file.txt -exec mv {} $OUTPUT_DIR \;

# Remove all file that are not needed
#rm 0
rm -f "$PREFIX"_temporaryfile
rm -f "$PREFIX"_temporaryfile1
rm -f "$F1"
rm -f "$F2"
rm -f "$F4"
rm -f "$F5" 
#rm -f "$F0" 

# Moving all produced files to the created and assigned folder
# All files to the sub directory
find $pwd -maxdepth 1 -name "$PREFIX*" -exec mv {} $OUTPUT_DIR \; #"$WORKING_DIR/Project_Results/$RIGHT_NOW" \;
#} 1>/dev/null 2>&1

printf "End\n"
exit
