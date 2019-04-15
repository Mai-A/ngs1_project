
source activate ngs1
mkdir ~/workdir/ngs_ass
cd ~/workdir/ngs_ass

#Data Download
wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR879/SRR8797509/SRR8797509.sra

#Prepare the data
#Part1
seqkit split2 -1 ./SRR8797509_X.part_001/SRR8797509_1.part_001.fastq -2 ./SRR8797509_X.part_001/SRR8797509_2.part_001.fastq -p 5 -O ./out_unshuffled5 -f
#Part2
seqkit split2 -1 ./shuffled_srr8797509/shuffled_SRR8797509_1.part_001.fastq -2 ./shuffled_srr8797509/shuffled_SRR8797509_2.part_001.fastq -p 5 -O ./out_shuffled5 -f

#FASTQ Quality Control

fastqc -t 1 -f fastq -noextract ./out_shuffled5/shuffled_SRR8797509_2.part_001.part_002.fastq -O ./fastqfile
fastqc -t 1 -f fastq -noextract ./out_shuffled5/shuffled_SRR8797509_1.part_001.part_001.fastq -O ./fastqfile
fastqc -t 1 -f fastq -noextract ./out_unshuffled5/SRR8797509_1.part_001.part_001.fastq -O ./fastqfile
fastqc -t 1 -f fastq -noextract ./out_unshuffled5/SRR8797509_2.part_001.part_001.fastq -O ./fastqfile
fastqc -t 1 -f fastq -noextract ./out_shuffled5/shuffled_SRR8797509_2.part_001.part_002.fastq -O ./fastqfile
fastqc -t 1 -f fastq -noextract ./out_shuffled5/shuffled_SRR8797509_1.part_001.part_001.fastq -O ./fastqfile
fastqc -t 1 -f fastq -noextract ./out_unshuffled5/SRR8797509_1.part_001.part_001.fastq -O ./fastqfile
fastqc -t 1 -f fastq -noextract ./out_unshuffled5/SRR8797509_2.part_001.part_001.fastq -O ./fastqfile
multiqc -z -o ./fastqfile ./fastqfile

#Trimming
#Mild Trimming for SX_1. {unshuffled}
f1="$home/workdir/ngs_ass/out_unshuffled5/SRR8797509_1.part_001.part_001.fastq"
f2="$home/workdir/ngs_ass/out_unshuffled5/SRR8797509_2.part_001.part_001.fastq"
newf1="$home/workdir/ngs_ass/unshuffled_trim/SRR8797509_1.part_001.part_001.fastq"
newf2="$home/workdir/ngs_ass/unshuffled_trim/SRR8797509_2.part_001.part_001.fastq"
newf1U="$home/workdir/ngs_ass/unshuffled_trim/SRR8797509_1.part_001.part_001.fastq"
newf2U="$home/workdir/ngs_ass/unshuffled_trim/SRR8797509_2.part_001.part_001.fastq"

adap="/home/ngs-01/miniconda3/envs/ngs1/share/trimmomatic-0.38-1/adapters"
trimmomatic PE -threads 1 -phred33 -trimlog trimLogFile -summary statsSummaryFile  $f1 $f2 $newf1 $newf1U $newf2 $newf2U \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:5:5 MINLEN:36

# aggresive trimming 

f3="$home/workdir/ngs_ass/out_shuffled5/shuffled_SRR8797509_1.part_001.part_001.fastq"
f4="$home/workdir/ngs_ass/out_shuffled5/shuffled_SRR8797509_2.part_001.part_001.fastq"
newf3="$home/workdir/ngs_ass/shuffled_trim/shuffled_SRR8797509_1.part_001.part_001.fastq"
newf4="$home/workdir/ngs_ass/shuffled_trim/shuffled_SRR8797509_2.part_001.part_001.fastq"
newf3U="$home/workdir/ngs_ass/shuffled_trim/shuffled_SRR8797509_1.part_001.part_001.fastq"
newf4U="$home/workdir/ngs_ass/shuffled_trim/shuffled_SRR8797509_2.part_001.part_001.fastq"

adap="/home/ngs-01/miniconda3/envs/ngs1/share/trimmomatic-0.38-1/adapters"

trimmomatic PE -threads 1 -phred33 -trimlog trimLogFile2 -summary statsSummaryFile2  $f3 $f4 $newf3 $newf3U $newf4 $newf4U \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:5:33 MINLEN:36



#BWA Alignment
cd ~/workdir/sample_data
wget https://de.cyverse.org/dl/d/A9330898-FC54-42A5-B205-B1B2DC0E91AE/dog_chr5.fa.gz
gunzip dog_chr5.fa.gz

mkdir -p ~/workdir/ngs_ass/bwa_align/bwaIndex && cd ~/workdir/ngs_ass/bwa_align/bwaIndex
ln -s ~/workdir/sample_data/dog_chr5.fa .
bwa index -a bwtsw dog_chr5.fa

cd ~/workdir/ngs_ass/bwa_align
for x in 1 2 3 4 5;do
R1="$home/workdir/ngs_ass/out_unshuffled5/SRR8797509_1.part_001.part_00$x.fastq"
R2="$home/workdir/ngs_ass/out_unshuffled5/SRR8797509_2.part_001.part_00$x.fastq"
bwa mem bwaIndex/dog_chr5.fa $R1 $R2 > SRRFastqunshuffled.sam;done

#HISAT
cd ~/workdir/sample_data
wget http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa

wget http://genomedata.org/rnaseq-tutorial/annotations/GRCh38/chr22_with_ERCC92.gt 

mkdir -p ~/workdir/ngs_ass/hisat_align/hisatIndex && cd ~/workdir/ngs_ass/hisat_align/hisatIndex
ln -s ~/workdir/sample_data/chr22_with_ERCC92.fa .
hisat2_extract_splice_sites.py ~/workdir/sample_data/chr22_with_ERCC92.gtf > splicesites.tsv
hisat2_extract_exons.py ~/workdir/sample_data/chr22_with_ERCC92.gtf > exons.tsv
hisat2-build -p 1 --ss splicesites.tsv --exon exons.tsv chr22_with_ERCC92.fa chr22_with_ERCC92

cd ~/workdir/ngs_ass/hisat_align
for c in 1 2 3 4 5;do
R1="$home/workdir/ngs_ass/out_shuffled5/shuffled_SRR8797509_1.part_001.part_00$c.fastq"
R2="$home/workdir/ngs_ass/out_shuffled5/shuffled_SRR8797509_2.part_001.part_00$c.fastq"

hisat2 -p 1 -x hisatIndex/chr22_with_ERCC92 --dta --rna-strandness RF -1 $R1 -2 $R2 -S UHR_Rep$C.sam;done


