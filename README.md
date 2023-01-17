# AHATool-container
Container for AHATool an Automatic HMM and Analysis Tool.

V.0.12 [26.Jan 2021]
Microbiology and Biotechnology - Streit's lab
University of Hamburg (D)
Developed by Nele Schulte (and P.Pérez-García)

Maintainer Albert Cañellas-Sole albert.canellas@bsc.es


## Introduction
AHATool.sh - Automatioc Hmmsearch and Analysis Tool
Version: 0.12 - 26.Jan 2021
Author: Nele Schulte

The files and the tool in this directory provide an automatisation of the HMM
creation an search. As well as the summary of the results. The MSA and HMM are
built using the NCBI non-redundant database.
(ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/)

Please do not attempt to view the non-redundant database file with a text 
editor since it contain millions of lines.



## Execution instructions
The tool can be executed using no parameters and thereby relying on the default
settings or inserting customized information on how where and what to process.

```
flags:
	-p,--prefix:  The prefix the tool will use for produced files.
        (default: 'yyddmmhhmm') [1]
  	-i,--inputFASTA:  the input FASTA file. (default: 'pet.fasta')
  	-s,--start:  start of execution (search or build). (default: 'build')
  	-d,--database:  database options: 1. nr_db; 2. custom_db (default: 'nr.fa')
  	-e,--evalue:  e-value (recommended: 1e-10). (default: 0.0000000001)
  	-t,--threads:  processor options: 1, 2, 4 (default: 2)
  	-h,--help:  show this help (default: false)

	[1] as in '2101271124' refering to 27.01.20201 - 11:24 a.m.
```
### Execution with container:
In order to use this software the user will need to create a volumme into container with the folder that have the database.
And to run the container with docker use the following command:
´´´
docker run -v /folder/with/database:/home/projects AHATool.sh -i input.fasta -d nr.fa
´´´

### Examples:
- Example 0: No parameters (default settings)
bash AHATool.sh

- Example 1: Setting the prefix for poduced files. (default: name of user)
bash AHATool.sh -p nele_12.12.2020

- Example 2: Setting the directory to work in (default: current/ this directory)
bash AHATool.sh /home/[user]/Documents/Analysis_with_AHA 
-or-
bash AHATool.sh i PETases.fasta -e 0.000001 --noSignalP  -d nr.fa /home/[user]/Documents/Analysis_with_AHA

> NOTE: The order of the input flags and arguments is not important.
> 
> It is recomended nevertheless to input the directory to work at the end of the iput parameters as this is the only parameter that comes as an argument.
> 
> The directory to work in (with default: current/ this directory) is not indicated with a flag but stands alone.

- Example 3: Setting the input FASTA file (default: only fasta file in input directory)
bash AHATool.sh -i PETases.fasta
-or-
bash AHATool.sh -i PETases.fa
-or-
bash AHATool.sh -i PETases.faa
-or-
bash AHATool.sh /home/[user]/Documents/Analysis_with_AHA -i PETases.fasta
-or-
bash AHATool.sh -i PETases.fasta /home/[user]/Documents/Analysis_with_AHA

- Example 4: Selecting start of execution (hmmbuild or hmmsearch)
bash AHATool.sh -s build
-or-
bash AHATool.sh -s search
-or-
bash AHATool.sh /home/[user]/Documents/Analysis_with_AHA -i PETases.fasta -s build
-or-
bash AHATool.sh -i PETases.fasta /home/[user]/Documents/Analysis_with_AHA -s build
-or-
bash AHATool.sh i PETases.fasta -s search /home/[user]/Documents/Analysis_with_AHA 

- Example 5: Setting the e-value (default: 1e-10)
bash AHATool.sh -e 0.000000001
-or-
bash AHATool.sh -i PETases.fasta -e 0.000001 -s build

- Example 6: Choosing database (default: nr.fa)
(two options: 1. nr_db; 2. custom_db)
bash AHATool.sh -d custom
-or-
bash AHATool.sh -d nr.fa
-or-
bash AHATool.sh i PETases.fasta -e 0.000001 --noSignalP  -d nr.fa

- Example 7: Usage of threads: Multithreadding (default: 2)
bash AHATool.sh -t 4
-or-
bash AHATool.sh -t 2
-or-
bash AHATool.sh i PETases.fasta -e 0.000001 -t 1

>NOTE: SignalP does not predict TAT (Twin-arginine translocation) signal peptides very well and signal peptides of bacterial lipoproteins not at all.
>
>Not all secretory proteins carry signal peptides. There can be non-classical secretory pathways or no sequence motif at all.
>Furthermore Gram-negative bacteria are known to have several secretion systems (type I, III, IV and VI) that function without signal peptides. 
>
>For prediction of such proteins use the SecretomeP server (http://www.cbs.dtu.dk/services/SecretomeP/) or the 
>MacSyFinder (https://galaxy.pasteur.fr/root?tool_id=toolshed.pasteur.fr/repos/odoppelt/txsscan/TXSScan/1.0.2).



## Files in this directory
Two kinds of files are available in this directory system. 

The first kind is explanatory and example files:

'pet.fasta'
Examplary input file.

'pet.fasta.aln'
Examplary MSA file.

'pet.fasta.html', 'pet.fasta.hmm', 'pet.fasta.dnd' and 'pet.fasta.fasta'
Exemplary HMM files.

The second kind of files is essential for running the tool (related to
the tool or the tool itself):

'nr.fa'
The non-redundant database. An update of this db is executed by the tool.

'nr.fa.ssi'
Index file for the nr db. An update of this index is executed by the tool.

'custom_db.fa'
The custom database. An update of this db has to be performed by the user.

'custom_db.ssi'
Index file for the custom db. Update has to be performed by the user.
example how to make ssi file

'AHATool.sh'
The tool itself. (Main code)

'SOFTWAREneeded.txt'
List of software needed for the correct execution of the tool.

'update_blastdb.pl'
Script used by 'Commands.sh' to download the pre-formatted nr db.

'shflags'
Command-line flag library.

Only the non-redundant database file has to be updated regularly. The tool will
exclaim a recommendation for the update when necessary and possible.

Each version of the database files is paired with an index file. For the nr db 
it is created when the download and entire update of the non-redundant 
database file has been a success. Namely: 'nr.fa.ssi'
For the custom db a creation of the index file needs to be done by the user
themslefs via: 'esl-sfetch --index custom_db.fasta'



## Input files
'[input].fasta'
Input file.

Optional:
'[input].fasta.aln'
MSA file.

'[input].fasta.html', '[input].fasta.hmm', '[input].fasta.dnd' and '[input].fasta.fasta'
HMM files.



## Output files
'[input file]_log_file.txt'
Logfile; Documenting the console output of the run.

'[input file].dnd'
Output result from t_coffee MSA; File Type= GUIDE_TREE; Format= newick   

'[input file].html'
Output result from t_coffee MSA; File Type= MSA; Format= html

'[input file].aln'
Output result from t_coffee MSA; File Type= MSA; Format= aln

'[input file].aln.fasta'
Reformated file: From Clusal to FASTA.

'temporary_[no.].txt'
Temporary files that will be deleted after execution.

'SummaryTable.tsv' and 'SummaryTable.xls'
Summary Table with tab separated columns. 
As followed when using the nr.fa database:

	1.  Column: target name
	2.  Column: accession number
	3.  Column: query name
	4.  Column: accession version
	5.  Column: e-value
	6.  Column: score
	7.  Column: bias
	8.  Column: e-value
	9.  Column: score  
	10. Column: bias    
	11. Column: exp - expected no. of domains in the seq.
	12. Column: reg - no. of discrete regions identified by this posterior decoding step
	13. Column: clu - no. of regions that had to be subjected to stochastic trace back clustering
	14. Column: ov - no. of envelopes identified by stochastic trace back clustering that overlap with other envelopes
	15. Column: env - total no. of domain envelopes identified (either by the simple method or by stochastic trace back clustering) 
	16. Column: dom - no. of domains in envelopes (before any significance threshold)
	17. Column: rep - no. of domains satisfying reporting thresholds.
	18. Column: inc - no. of domains satisfying inclusion thresholds.
	19. Column: description of target
	20. Column: size (aa)
	21. Column: SignalP Prediction 
	22. Column: CS pos.
	23. Column: accession number
	24. Column: species
	25. Column: taxonomy
	26. Column: sequence

Differences when using a custom database: 

	19. Column: SignalP Prediction
	20. Column: CS pos.
	21. Column: size (aa)
	22. Column: Start
	23. Column: End
	24. Column: Direction
	25. Column and following are database specific
 
>NOTE: When using a custom database only column 01 to 18 are of identical order.

'[user name]_all_hits.fa'
FASTA file with extracts of arbitrary subsets of sequences or subsequences from hmmsearch. 

'[user name]_partial_hits.fa'
FASTA file with hits in '[user name]_all_hits.fa' that have partial=[10,01,11].
-> Only when working with a custom database.

'[user name]_entire_hits.fa'
FASTA file with hits in '[user name]_all_hits.fa' that have partial=00.
-> Only when working with a custom database.

'[user name]_[input file].hmm'
A profile HMM created from [input file] by hmmbuild.

'[user name]_[input file].hmm.aln'
Stockholm markup to pick up a variety of information. It saves all the domains (subsequences) that were significant (passed the "inclusion thresholds") as one alignment, in Stockholm format.

'[user name]_[input file].hmm.out'
The output alignment file.

'[user name]_[input file].hmm.tbl'
This is the target hits table produced by the hmmer output option --tblout.
It consists of one line for each different query/target comparison that met the reporting thresholds, ranked by increasing E-value (decreasing statistical significance).

'[user name]_gram+_summary.signalp5' and '[user name]_gram-_summary.signalp5' and '[user name]_arch_summary.signalp5'
This is a tabular file containing the protein prediction, the associated likelihood probability (LP) for the four signal peptides (SP(Sec/SPI), LIPO(Sec/SPII), TAT(Tat/SPI), OTHER) and the cleavage site position with associated LP. 
NOTE: if the cleavage site position is "?", it means that the cleavage site is out of range due to a probable protein fragment as input.)

'[user name]_arch_mature.signalp5' and '[user name]_gram+_mature.signalp5' and '[user name]_gram-_mature.signalp5'
FASTA file with hit sequences edited such that a FASTA file with only mature sequences emerges. (signal peptide removed from sequences)
Mature FASTA file only contains the sequences for which signal peptides could be predicted.
________________________________________________________________________________
Universität Hamburg
Faculty of Mathematics, Informatics and Natural Sciences
Department Biology
Institute of Plant Science and Microbiology
Microbiology and Biotechnology
Ohnhorststr. 18
22609 Hamburg
Germany
________________________________________________________________________________
