#!/bin/bash
# author Vanessa Vermeirssen
#This pipeline preprocesses RNA-seq data from a xenograft experiment, where there is a combination of human and murine RNA.
# exRNAxeno_combined makes use of a combined human and mouse reference genome for the read mapping
#Based on PicoV2 paired-ended total-RNA-seq (stranded)
#SMARTer® Stranded Total RNA-Seq Kit v2 - Pico Input Mammalian - library prep
#Run this script using "bash pipeline_exRNAxeno_part2_combined.sh <datadir> <outputdir>"
##folders and files
DATADIR="$1" #read from the command line, 1st argument = input directory of the raw fastq files and folders
DIR="$2" #read from the command line, 2nd argument = output directory of the processed data
## mkdir $DIR ## if directory does not exist yet
cd $DIR
INDEX="../STAR_index94/STAR_combi_index_humanmouse" ##including Sequin and ERCC spike-ins and rDNA
GTF="../STAR_index94/HumanMouse_GRC38.94_spikes.gtf" ##including Sequin and ERCC spike-ins
#masking based on control plasma samples of mouse and human
FILTER_mouseCh_humanCm="../masking/humanmousecontrol_C.bed"

#subsampling size
subsampling_size="491036" #downsample to the smallest sequencing depth

#collect all sample folders (starting with RNA0) in the input data directory
find $DATADIR -maxdepth 1 -type d -printf '%f\n' | grep "RNA0" | sort -n > listsamplefolders.txt

#make job files to be submitted on the computing cluster
for samplefolder in $(cat listsamplefolders.txt)
do
	sample=`echo ${samplefolder:0:9}` #extract the sample name from the samplefolder e.g. RNA004660 from foldername RNA004660L1-205346173

	filename="combined_${sample}.sh"
  echo "#!/bin/bash" > $filename
  echo "#SBATCH -J $sample" >> $filename ## jobname
  echo "#SBATCH -D $DIR" >> $filename
  echo "#SBATCH --cpus-per-task=18" >> $filename ## to run on computing cluster
	echo "#SBATCH --time=72:00:00" >> $filename
  # echo "#SBATCH --mem=100G" >> $filename
	echo "#SBATCH -o combined_${sample}.out" >> $filename
	echo "#SBATCH -e combined_${sample}.err" >> $filename
  echo "" >> $filename
  # echo "mkdir $DIR/$samplefolder" >> $filename ## if directories do not exist yet
  echo "cd ${DIR}/$samplefolder" >> $filename

	##subsampling using seqtk
	#-s100: sets the seed to a fixed number (here 100, but may be any number) in order to have the same result if you rerun this part of the code
  echo "module purge" >> $filename
  echo "module load seqtk/1.3-foss-2018a" >> $filename
  echo "seqtk sample -s100 ${sample}_1clumped.fastq $subsampling_size > ${sample}_1_subsampled.fastq" >> $filename
  echo "seqtk sample -s100 ${sample}_2clumped.fastq $subsampling_size > ${sample}_2_subsampled.fastq" >> $filename
  echo "" >> $filename

	##FastQC for each sample (before and after preprocessing), can be combined with multiQC later on locally
 	echo "module purge" >> $filename
 	echo "module load FastQC/0.11.8-Java-1.8" >> $filename
 	echo "fastqc ${sample}_1_subsampled.fastq" >> $filename
 	echo "fastqc ${sample}_2_subsampled.fastq" >> $filename
 	echo "fastqc ${sample}_1.fastq" >> $filename
 	echo "fastqc ${sample}_2.fastq" >> $filename
 	echo "mkdir FASTQC" >> $filename
 	echo "mv *_fastqc* FASTQC/." >> $filename
 	echo "" >> $filename

	##STAR mapping
	echo "module purge" >> $filename
  echo "ml STAR/2.6.0c-intel-2018a" >> $filename
	echo "mkdir ${sample}_sroutC/" >> $filename
  echo "STAR --runThreadN 10 --outFileNamePrefix ${sample}_sroutC/ --readFilesIn ${sample}_1_subsampled.fastq ${sample}_2_subsampled.fastq --genomeDir ${INDEX} --sjdbGTFfile ${GTF} --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --quantMode GeneCounts TranscriptomeSAM --outSAMprimaryFlag AllBestScore" >> $filename
	# --outSAMmultNmax -1 --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 == these are default options
	# multimapping allowed: --outSAMprimaryFlag AllBestScore and --outSAMmultNmax -1: all multiple alignments with the best score will be output as 'primary alignments'

	##filtering for uniquely mapped reads
	#based on the NH tag of the SAM/BAM file, NH:i:1 are uniquely mapped reads
	#By default, STAR outputs NH HI AS nM SAM attributes
	echo "module purge" >> $filename
	echo "module load SAMtools/1.8-intel-2018a" >> $filename
	echo "samtools view -h ${sample}_sroutC/Aligned.sortedByCoord.out.bam | grep -P '^@|NH:i:1\b' | samtools view -Sb - > ${sample}_sroutC/uniquelymappedC.bam" >> $filename

	##masking based on control mouse only or human only plasma samples
	echo "module purge" >> $filename
	echo "module load SAMtools/1.8-intel-2018a" >> $filename
	echo "module load BEDTools/2.27.1-intel-2018a" >> $filename
	echo "bedtools intersect -v -abam ${sample}_sroutC/uniquelymappedC.bam -b ${FILTER_mouseCh_humanCm} > ${sample}_sroutC/uniquelymapped_filteredC.bam" >> $filename

	##sort bam by read name for HTSeq-count
	echo "module purge" >> $filename
	echo "module load picard/2.21.1-Java-11" >> $filename
	echo "java -jar \${EBROOTPICARD}/picard.jar SortSam I=${sample}_sroutC/uniquelymapped_filteredC.bam O=${sample}_sroutC/uniquelymapped_filteredC_sortedByName.bam SORT_ORDER=queryname" >> $filename

	##total gene read counts by HTSeq-count
	#(default) needs name sorted bam files
	echo "module purge" >> $filename
	echo "module load HTSeq/0.11.0-foss-2018b-Python-2.7.15" >> $filename
	#htseq-count default non-unique ignores all reads mapping to multiple features - so reads mapping equally well to human and mouse are not counted
	#htseq-count parameters: For stranded=yes, the read has to be mapped to the same strand as the feature. For paired-end reads, the first read has to be on the same strand and the second read on the opposite strand. For stranded=reverse, these rules are reversed.
	#SMARTer Stranded Total RNA-seq: read 1 = antisense, read 2 = sense to the original RNA
	echo "htseq-count -f bam -r name -s reverse --additional-attr=gene_name --secondary-alignments=ignore --supplementary-alignments=ignore ${sample}_sroutC/uniquelymapped_filteredC_sortedByName.bam ${GTF} > countspergeneC.tab" >> $filename

############
##quality control of alignment
##strandedness correctness RSeQC
  echo "module purge" >> $filename
  echo "module load Python/2.7.14-intel-2018a" >> $filename
  echo "module load bx-python/0.8.1-intel-2018a-Python-2.7.14" >> $filename
  echo "module load RSeQC/2.6.4-intel-2018a-Python-2.7.14" >> $filename #python script infer_experiment.py is present in this module
	echo "infer_experiment.py -r ${exon} -i ${sample}_sroutC/uniquelymapped_filteredC.bam > RSeQC_outputC.txt" >> $filename
	echo "outC=\`cat RSeQC_outputC.txt | grep '1+-' | cut -d':' -f2\`" >> $filename
	echo "echo ${sample} \$outC > RSeQC_C.txt" >> $filename
  echo "" >> $filename

#concatenate the RSeQC results for each sample - you need to do this at the end!
# echo "cat ${DIR}/*RNA*/RSeQC_C.txt > ${DIR}/StrandedC.txt" >> $filename

##alignment summary statistics using samtools idxstats
	echo "module purge" >> $filename
  echo "module load SAMtools/1.8-intel-2018a" >> $filename
	echo "samtools index ${sample}_sroutC/uniquelymapped_filteredC.bam" >> $filename
	echo "samtools idxstats ${sample}_sroutC/uniquelymapped_filteredC.bam > idx_statC.txt" >> $filename

#remove original fastq files in the output directory
	echo "rm *_1.fastq" >> $filename
  echo "rm *_2.fastq" >> $filename

done

##submit the jobs on the computing cluster
module swap cluster/skitty #select the computing cluster
for samplefolder in $(cat listsamplefolders.txt)
do
   cd $DIR #to make sure scripts and error files are put in output directory
   sample=`echo ${samplefolder:0:9}`
   echo $sample
   sbatch combined_${sample}.sh
done
