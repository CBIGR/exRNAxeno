#!/bin/bash
# author Vanessa Vermeirssen
#This pipeline preprocesses total RNA-seq data from a xenograft experiment, where there is a combination of human and murine RNA.
#Based on PicoV2 paired-ended total-RNA-seq (stranded)
#SMARTer® Stranded Total RNA-Seq Kit v2 - Pico Input Mammalian - library prep
#Read 1 generates sequences antisense to the original RNA
#Read 2 contains 3 nucleotides from the adapter (template switching) that should be trimmed
#Run this script using "bash pipeline_exRNAxeno_part1.sh <datadir> <outputdir>"
##folders and files
DATADIR="$1" #read from the command line, 1st argument = input directory of the raw fastq files and folders
DIR="$2" #read from the command line, 2nd argument = output directory of the processed data
# mkdir $DIR ## if directory does not exist yet
cd $DIR

#collect all sample folders (starting with RNA0) in the input data directory
find $DATADIR -maxdepth 1 -type d -printf '%f\n' | grep "RNA0" | sort -n > listsamplefolders.txt

#make job files to be submitted on the computing cluster
for samplefolder in $(cat listsamplefolders.txt)
do
	sample=`echo ${samplefolder:0:9}` #extract the sample name from the samplefolder e.g. RNA004660 from foldername RNA004660L1-205346173

  filename="prep_${sample}.sh"
  echo "#!/bin/bash" > $filename
  echo "#SBATCH -J $sample" >> $filename ## jobname
  echo "#SBATCH -D $DIR" >> $filename
  echo "#SBATCH --cpus-per-task=18" >> $filename ## to run on computing cluster
  echo "#SBATCH --time=4:00:00" >> $filename #this parameter was increased to 6 h for large samples
  echo "#SBATCH --mem=80G" >> $filename #this parameter was omitted and job was submitted on phanpy for large sample - clumpify grep command!
	echo "#SBATCH -o prep_${sample}.out" >> $filename
	echo "#SBATCH -e prep_${sample}.err" >> $filename
  echo "" >> $filename

  echo "mkdir $DIR/$samplefolder" >> $filename ## if directories do not exist yet
  echo "cd ${DIR}/$samplefolder || exit" >> $filename

  ##unzip and merge different lanes, so you end up with 2 paired-end reads per sample
  echo "cp -r ${DATADIR}/${samplefolder}/* ." >> $filename
	echo "gunzip *.fastq.gz" >> $filename
  echo "" >> $filename
  echo "cat *R1*.fastq > ${sample}_1.fastq" >> $filename
  echo "cat *R2*.fastq > ${sample}_2.fastq" >> $filename

	##trim template switching adapter with cutadapt
   #For paired-end reads: cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq
     # -U 3 -> remove a fixed number of bases from the second read in a pair
     # --pair-filter=any -> this dedault option means that a read pair is discarded (or redirected) if one of the reads fulfills the filtering criterion
     # -q 30 ->  trim low-quality bases from 3' end of each read before adapter removal. Applied to both reads if data is paired.
     # -m 35 -> discard reads shorter than LENGTH 35
     # -a and -A -> remove 3' adapter sequences from first and second read respectively, HT-TruSeq adapters - see manual SMARTer Stranded Total RNA-Seq
  echo "module purge" >> $filename
  echo "module load cutadapt/1.18-intel-2018b-Python-3.6.6" >> $filename
  echo "cutadapt -q 30 --pair-filter=any -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  -m 35 -U 3 -o ${sample}_1trim.fastq -p ${sample}_2trim.fastq ${sample}_1.fastq ${sample}_2.fastq" >> $filename
  ###

	##clumpify for removal of duplicated reads using BBMap
	# cutadapt: trim every read to length 60 (-l 60) and gzip file befor processing with clumpify
	# clumpify:
	# dedupe: Remove duplicate reads. For pairs, both must match. By default, deduplication does not occur (dedupe=f)
	# (default allowns=t: No-called bases will not be considered substitutions)
	# Duplicates are removed within clumps that are based on reads sharing a specific kmer.
	# subs: number of substitutions allowed (if subs=0, only exact matches are considered as duplicates, default 2, note that this allows 2 substitutions in both ends of the pair)
	# k: kmer size, reads that share a kmer of given length end in same clump (default 31), lower k will result in more reads per clump and thus more chance to eliminate duplicates
	# passes: number of times a different kmer is selected for seeding clumps (default 1), eventually, any pair of duplicates will land in same clump given enough passes if they share a single kmer
	# For plasma data (with a lot of duplicates), we have observed that the plateau is reached after about 10 passes but additional passes do not take much longer. Therefore we run with 20 passes.
	# See https://www.biorxiv.org/content/10.1101/2021.05.11.442610v1
	echo "module purge" >> $filename
	echo "module load BBMap/38.26-foss-2018b" >> $filename
	echo "module load cutadapt/1.18-foss-2018b-Python-3.6.6" >> $filename
	echo "mkdir $DIR/$samplefolder/clumpout" >> $filename
	echo "cutadapt --pair-filter=both -l 60 -o clumpout/${sample}_clumtrim1.fastq -p clumpout/${sample}_clumtrim2.fastq ${sample}_1trim.fastq ${sample}_2trim.fastq" >> $filename
	echo "gzip clumpout/${sample}_clumtrim*.fastq" >> $filename
	echo "clumpify.sh usejni=t in=clumpout/${sample}_clumtrim1.fastq.gz in2=clumpout/${sample}_clumtrim2.fastq.gz out=clumpout/${sample}_tempclumped_1.fastq.gz out2=clumpout/${sample}_tempclumped_2.fastq.gz dedupe subs=2 k=31 passes=20" >> $filename
	echo "gunzip clumpout/${sample}_tempclumped*.fastq.gz" >> $filename
	echo "grep '^@' clumpout/${sample}_tempclumped_1.fastq > clumpout/${sample}_tempnames_1.txt; grep '^@' clumpout/${sample}_tempclumped_2.fastq > clumpout/${sample}_tempnames_2.txt" >> $filename
	echo "grep -Ff clumpout/${sample}_tempnames_1.txt -A3 ${sample}_1trim.fastq | grep -v -- '^--$' > ${sample}_1clumped.fastq" >> $filename
	echo "grep -Ff clumpout/${sample}_tempnames_2.txt -A3 ${sample}_2trim.fastq | grep -v -- '^--$' > ${sample}_2clumped.fastq" >> $filename
	echo "rm -r clumpout/*fastq*" >> $filename

	##sequencing depth
	echo "#Count number of lines per fastq (to get the nr of reads, nr of lines should be divided by 4)" >> $filename
	echo "cat ${sample}_1clumped.fastq | wc -l  >> total_seq_reads_clumped.txt" >> $filename
	echo "cat ${sample}_2clumped.fastq | wc -l  >> total_seq_reads_clumped.txt" >> $filename
	echo "cat ${sample}_1.fastq | wc -l  >> total_seq_reads.txt" >> $filename
	echo "cat ${sample}_2.fastq | wc -l  >> total_seq_reads.txt" >> $filename
###########################
done

##submit the jobs on the computing cluster
module swap cluster/skitty #select the computing cluster
for samplefolder in $(cat listsamplefolders.txt)
do
   cd $DIR #to make sure scripts and error files are put in output directory
   sample=`echo ${samplefolder:0:9}`
   echo $sample
   sbatch prep_${sample}.sh
done
