---
output: html_document
---
A host-xenograft deconvolution algorithm to distinguish (human) tumoural from (murine) host RNA from liquid biopsies of human xenograft mouse models and as such to specifically analyse tumour-derived exRNA.    

The exRNAxeno framework consists of 1 general preprocessing pipeline and 2 deconvolution pipelines that respectively map to a combined reference genome of human and mouse, or map in parallel to the reference genomes of both species followed by a filtering step. Both deconvolution pipelines are followed by a masking step, where misaligned reads that were detected in human and mouse only samples, are removed.    

The pipelines are specifically designed for plasma samples analyzed by SMARTer® Stranded Total RNA-Seq Kit v2 - Pico Input Mammalian RNA library preparation and subsequent Illumina paired-end sequencing.  
They are tailored to be run on the UGent HPC computing cluster (Slurm Workload Manager).

# pipeline_exRNAxeno_part1.sh
This pipeline performs general preprocessing steps on raw fastq files. It outputs clumped fastq files and the sequencing depth per sample.   
Run this script using ```bash pipeline_exRNAxeno_part1.sh <datadir> <outputdir>```.  
If needed, adapt the computing cluster on line 90.  

1. Unzipping and merging fastq files of different lanes to end up with 2 paired-end reads per sample.
2. Trimming the template switching adapter and low quality bases at 3' end, as well as discarding too short reads by *cutadapt*.
3. Removing duplicate reads by *clumpify*. 
4. Analyzing the sequencing depth per sample.

# pipeline_exRNAxeno_part2_combined.sh
This is the deconvolution pipeline based on mapping to a combined reference genome of human and mouse. It inputs clumped fastq files and outputs a gene count table for human and mouse.  
Run this script using ```bash pipeline_exRNAxeno_part2_combined.sh <datadir> <outputdir>```.  
For STAR mapping, it requires a STAR index of the combined human and mouse reference genomes (line 13 ```INDEX="../STAR_index94/STAR_combi_index_humanmouse"```) and a similar GTF file (line 14 ```GTF="../STAR_index94/HumanMouse_GRC38.94_spikes.gtf"```).     
For masking, a control bed file is needed (line 16 ```FILTER_mouseCh_humanCm="../masking/humanmousecontrol_C.bed"```).  
If needed, adapt the subsampling size on line 19 and the computing cluster on line 124.  

1. Downsampling to the smallest sequencing depth using *seqtk*.
2. Quality control by *FASTQC* on the original, merged fastq files and subsampled fastq files, which can be combined with multiQC later on.
3. Mapping the reads to the combined reference genome by *STAR*.
4. Filtering for uniquely mapped reads based on the *NH:i:1* tag in the SAM/BAM files.
5. Masking for misaligned reads in human and mouse only samples using *intersectBed* and the BED file ```FILTER_mouseCh_humanCm="../masking/humanmousecontrol_C.bed"``` (line 16).
6. Generating gene counts for uniquely mapped reads by *HTSeq-count*.
7. Quality control of the alignment by *RSeQC* and *SAMtools idxstats*.

# pipeline_exRNAxeno_part2_parallel.sh
This is the deconvolution pipeline based on the parallel mapping to the reference genomes of human and mouse. It inputs clumped fastq files and outputs gene count tables for human and mouse.   
Run this script using ```bash pipeline_exRNAxeno_part2_parallel.sh <datadir> <outputdir>```.   
It requires an additional Python script ```selecthumanmouseonly.py``` that is here present in the folder ```scripts``` and should be put upstream of the output directory.
For STAR mapping, it requires STAR indexes of the human and mouse reference genomes (line 13 ```INDEX_human="../STAR_index94/STAR_index_human"``` and line 15 ```INDEX_mouse="../STAR_index94/STAR_index_mouse"```) and similar GTF files (line 14 ```GTF_human="../STAR_index94/Homo_sapiens.GRCh38.94_spikes.gtf"``` and line 16 ```GTF_mouse="../STAR_index94/Mus_musculus.GRCm38.94_spikes.gtf"```).  
For masking, control bed files for human and mouse are needed (line 18 ```FILTER_mouseHnotM="../masking/mousecontrol_HnotM_onlyhuman.bed"``` and line 19 ```FILTER_humanMnotH="../masking/humancontrol_MnotH_onlymouse.bed"```).
If needed, adapt the subsampling size on line 22 and the computing cluster on line 150.  

1. Downsampling to the smallest sequencing depth using *seqtk*.
2. Quality control by *FASTQC* on the original, merged fastq files and subsampled fastq files, which can be combined with multiQC later on.
3. Mapping the reads to both reference genomes in parallel by *STAR*.
4. Filtering for uniquely mapped and `human preferred over mouse` as well as `mouse preferred over human` reads based on the edit distance i.e. the sum of the *NM tag* and the *soft-clipping in the CIGAR string* from both read pairs in the SAM/BAM file.
5. Masking for misaligned reads in human and mouse only samples using *intersectBed* and the BED files ```FILTER_mouseHnotM="../masking/mousecontrol_HnotM_onlyhuman.bed"FILTER_mouseHnotM="../masking/mousecontrol_HnotM_onlyhuman.bed"``` (line 18) and ```FILTER_humanMnotH="../masking/humancontrol_MnotH_onlymouse.bed"``` (line 19).
6. Generating gene counts for uniquely mapped reads by *HTSeq-count*.
7. Quality control of the alignment by *RSeQC* and *SAMtools idxstats*.
