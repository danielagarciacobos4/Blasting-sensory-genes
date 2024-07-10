## AIM OF BLASTING

The goal of these notes is to document the steps leading up to BLASTing opsin genes in different snake species of the Thamnophini tribe that differ in their habitats (eg. aquatic, terrestrial, and semi-fossorial). This script was useful for locating the coordinates of the opsin genes in the genome, however, I was unable to extract the CDS regions for selection analysis. To extract CDS regions I had to use haLiftover, using the HAL file that resulted from the CACTUS alignment (genome assemblies and alignments where done by Leroy Nunez using the annotated genome of _Thamnophis elegans_ as reference uploaded in [NCBI](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/35005/). 

## General Information and File Setup

Here, we will outline where the assembled genome is stored on the American Museum of Natural History's Huxley cluster. 
Relevant file locations for these analyses:

_Genomes- from HAL file_
```
/home/dgarcia/nas4/phd/HAL_Thamnophini/Fastas
```
_Snakes opsins database (step #1)_
```
/home/dgarcia/nas4/phd/opsins/BLAST/initial_Blast/proteins_database
```
_Scripts for initial BLAST (step ##)_
```
/home/dgarcia/nas4/phd/opsins/BLAST/initial_Blast/Scripts
```
_Results of intitial BLAST (step ##)_
```
/home/dgarcia/nas4/phd/opsins/BLAST/initial_Blast/Results
```
_List with best BLAST hits from initial BLAST_
```
/home/dgarcia/nas4/phd/opsins/BLAST/second_Blast/lists
```
_Scripts for second BLAST_
```
/home/dgarcia/nas4/phd/opsins/BLAST/second_Blast/scripts
```
_Results for second BLAST_
```
/home/dgarcia/nas4/phd/opsins/BLAST/second_Blast/results
```
_BITS extracted from genome_
```
/home/dgarcia/nas4/phd/opsins/BLAST/second_Blast/BITS/initial_bits
/home/dgarcia/nas4/phd/opsins/BLAST/second_Blast/BITS/final_bits
```

## Step 1) Create a curated database of snake opsins genes

Script for pulling these genes from NCBI using [the edirect package](https://www.nlm.nih.gov/dataguide/edirect/documentation.html) in commandline. This manual is aslo helpfull [Entrez Direct Examples](https://www.ncbi.nlm.nih.gov/books/NBK565821/) This script uses the esearch function to query for particular taxa or gene name, and then pipes the search results into the efetch function, which saves the results into a fasta locally. Our cluster at AMNH does not have the edirect package, so we has to install it in our miniconda environment. Just a reminder of how to use [PBS scripts](https://latisresearch.umn.edu/creating-a-PBS-script).

```
# script for pulling opsin gene protein seqs from NCBI for tblastn step within the Huxley cluster- AMNH (using a PBS submission language)

#!/bin/bash
#PBS -q batch
#PBS -l mem=30gb
#PBS -m abe
#PBS -M dgarcia@amnh.org
#PBS -q batch
#PBS -N opsin
#PBS -l walltime=6:00:00

source ~/.bash_profile
conda activate edirect

# script for pulling opsin gene protein seqs from NCBI for tblastn step
## Retrieve opsin protein and mRNA sequences from NCBI

cd /home/dgarcia/nas4/thamnophini_genomes/protein_fasta_files

esearch -db protein -query "opsin AND snakes [ORGN]" | efetch -format fasta > opsin_protein_v2.fasta
```
Important notes to consider about the previous script:
- We are filtering the opsin gene matches for snakes by the command "opsin AND snakes [ORGN]"
- We are indicating the path where we are going to save our results by inserting "cd /home/dgarcia/nas4/thamnophini_genomes/protein_fasta_files"
- We can count the number of matches by running "grep ">" -c opsin_protein_v2.fasta". In this case we found a match of 469
- Our outcome look something like this:
  
<img width="594" alt="Screenshot 2024-02-06 at 5 07 51 PM" src="https://github.com/danielagarciacobos4/Blasting-sensory-genes/assets/67153479/c9a64bdf-3af9-47fc-83c3-8f372265ed3e">

## Step 2) Change our protein fasta file into a database file

Blast will not run against a fasta file. This is why we have to transform our fasta file into a database using the following command: 

```
makeblastdb -in <database>.fa -dbtype prot
makeblastdb -in /home/dgarcia/nas4/thamnophini_genomes/protein_fasta_files/opsin_protein_snakes.fasta -dbtype prot

```
The output of this command will be different types of files. Something like this: 

<img width="726" alt="Screenshot 2024-02-07 at 12 05 25 PM" src="https://github.com/danielagarciacobos4/Blasting-sensory-genes/assets/67153479/25c0286c-905d-4693-bb04-238042438cf4">


## Step 3) Run blast Genome (query) vs Opsins (database)

This is an example using the _Liodytes alleni_ genome and opsin protein database. Scripts name: L.alleni_ib.sh
_Credits to Joe Arguelles for elaborating this script_
```
#!/bin/bash
#PBS -q batch
#PBS -l select=1:ncpus=25
#PBS -m abe
#PBS -o /home/dgarcia/nas4/phd/opsins/BLAST/initial_Blast/TMP
#PBS -e /home/dgarcia/nas4/phd/opsins/BLAST/initial_Blast/TMP
#PBS -M dgarcia@amnh.org
#PBS -N L.alleni_IB

#PBS -l walltime=9999:00:00


#SET VARIABLES
WORKDIR=/home/dgarcia/nas4/phd/opsins/BLAST/initial_Blast/Results #working directory
SPECIES="L.alleni_June16" #species name/abbreviation for output
DB=/home/dgarcia/nas4/phd/opsins/BLAST/initial_Blast/proteins_database/opsin_protein_snakes.fasta #path to nr database (if local) or just 'nr' if running online
NAME="initialBLAST" #keyword for BLAST output
IN=/home/dgarcia/nas4/phd/HAL_Thamnophini/Fastas/Liodytes_alleni.fasta #infile, auto-generated
OUT=${WORKDIR}/${SPECIES}_${NAME}.txt #output name, auto-generated


cd $WORKDIR
module load ncbi-blast-2.12.0+

blastx -query $IN -db $DB -out $OUT \
-evalue 1e-10 -max_target_seqs 100 \
-outfmt "6 qseqid sseqid pident evalue bitscore qstart qend sstart send sacc stitle"

```

Things to consider for the previous script: 
- TMP folder will hold all my output errors and notifications. This way, once the job is finished I can just eliminate this folder (#PBS -o and #PBS -e),
- outfmt 6 is a type of format in which my results will be executed. For more information search [the manual](https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.options_common_to_all_blast/)
- more on the output format 6 [outfmt 6](https://www.metagenomics.wiki/tools/blast/blastn-output-format-6). Note for next time I do blast: include the "sacc" command to get more straight forward genes names in the output.

An example of number of hits for some of the species I ran the code for:
- _Liodytes alleni_: 1554 hits
- _Nerodia clarkii_: 1533 hits
- _Clonophis kirtlandii_: 1529 hits
- _Natrix natrix_: 1559 hits 
- _Regina grahamii_: 1522 hits
- _Tropidoclonion lineatum_: 1539

## Step 4) Filter the best hits from step 3 (sort and eliminate duplicated coordinates)
The ~ 1500 hits in the genomes will contain duplicated coordinated or very nearby coordinates. We want to sort and reduce the best unique coordinates. For this, we need to sort and eliminate duplicates. The following steps are done in EXCEL (ideally we can find another way to optimize this step). 
  1) Rounded to 1000 the start coordinate =MROUND(cell,1000)
  2) Sort the file by: a) Chromosome name, b) rounded coordinate, c) best bit score (largest to smallest), d) e-value (smallest to largest)
  3) Eliminate all duplicates considering the rounded coordinate
  4) Move 




## Second blasting tblastn
As we can see, our results show a great amount of hits (around 7,000), which makes it difficult to analyze. For this reason, we will try to filter the amount of hits by doing a different search in BLAST (recommended by Jeff W.). In this second approach we will change two main things: 
- 1) We will use the snake's opsin protein fasta file as the query and the whole genome assembly as the database (in the steps before, we were using the genome as query and protein data files as the database)
- 2) We will perform a tblastn command (in the steps before, we were using an blastx command). Just to have a reminder of the possible blast searches available dependings of the datasets (this image is not mine, [credits] https://slideplayer.com/slide/13408600/): 
![Screenshot 2024-02-23 at 6 51 31 PM](https://github.com/danielagarciacobos4/Blasting-sensory-genes/assets/67153479/213dddb9-9ade-4b8f-af6f-395576110481)


```
makeblastdb -in /home/dgarcia/nas4/thamnophini_genomes/Nerodia_clarkii_AMNH_R500948_sma.fa -dbtype nucl
```
##include a reminder to search for cds regions with the reverse compliment. Could this be the reason it did not work before?

. 
. 
. 
. 
. 
. 
. 



_I will leave this here in case I have to go back to this script that will grab all the genes identified for Natricidae_

```
# script for pulling Nerodia gene protein seqs from NCBI for tblastn step

#!/bin/bash
#SBATCH --job-name=ncbi_download
#SBATCH -o %A_%a.2024-02-02_NCBI_Download.out
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH -c 1
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 02:00:00
#SBATCH --account=kawahara
#SBATCH --qos=kawahara

module load edirect/12.2
module load seqkit

# Enable auto abbreviation resolution
export EDIRECT_DO_AUTO_ABBREV="true"

# Set the NCBI API key
export NCBI_API_KEY="25de463246b0746c402efd4ffdc45f70fe09"

## Retrieve opsin protein and mRNA sequences from NCBI
esearch -db protein -query "Natricidae [ORGN]" -api_key $NCBI_API_KEY | efetch -format fasta > Natricidae_protein.fasta
```


