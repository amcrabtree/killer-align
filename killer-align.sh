#!/bin/bash

## alignment and consensus sequence pipeline
## Author: Angela Crabtree

########## OPTIONS ############
while getopts "s:1:2:r:o:pth" opt; do
	case ${opt} in
		s) sample=$OPTARG ;;
		1) fread=$OPTARG ;;
		2) rread=$OPTARG ;;
		r) genomes=$OPTARG ;;
		o) usroutdir=$OPTARG ;; 
		p) fastp=1 ;;
		t) mode="test" ;;
		h)
			printf "\n\n-------------------------------------------------------\n"
			printf "\nOptions:\n"
			printf "\n"
			printf "   -s [arg]	sample name (required)\n"
			printf "   -1 [arg]	forward read file (required)\n"
			printf "   -2 [arg]	reverse read file (required)\n"
			printf "   -r [arg]	reference genome file (required)\n"
			printf "   -o [arg]	destination of output folder (required)\n"
			printf "   -p		clean/trim reads with fastp before aligning\n"
			printf "   -t		test (ensures required CL apps are working)\n"
			printf "   -h		help\n"
			printf "\n-------------------------------------------------------\n\n\n"
			exit 0
			;;
	esac
done

########## ESTABLISH BACKGROUND LOGISTICS ############

## LOAD MODULES (SERVER USE ONLY)
local=$(module load 2>&1 >/dev/null | grep 'command not found' | wc -l)
if [[ $local != 1 ]]; then
	module load java
	module load python
	module load bwa
	module load samtools
	module load R
fi

## SET UP DIRECTORIES
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"   # location of this script file
if [ -n "$mode" ]; then
	printf "\n\tOutput files will be saved to your current directory.\n"
	usroutdir=$(pwd)
	genome="${scriptdir}"/test/M1.fasta
	sample="TEST_STRAIN"
	fread="${scriptdir}"/test/NCYC190_S3_L001_R1_001.fastq.gz
	rread="${scriptdir}"/test/NCYC190_S3_L001_R2_001.fastq.gz
fi	
appdir=~/bin
outdir="${usroutdir}"/"${sample}"/mapping
mkdir -p ${outdir}/stats ${outdir}/img ${outdir}/logs ${outdir}/vcf ${outdir}/bam
report_file=${outdir}/${sample}_mapping-report.md
report_html=${outdir}/${sample}_mapping-report.html
bcflog=${outdir}/logs/bcftools.log
fastp_html=${outdir}/stats/${sample}_fastp.html

#################### OPTIONAL READ TRIMMING WITH FASTP ###################

if [[ $fastp = 1 ]]; then
	printf "\n\n*** Cleaning reads with Fastp ***\n\n"
	# run fastp
	fastp \
		--trim_poly_x --trim_front1 6 --trim_tail1 6 \
		-i "${fread}" \
		-I "${rread}" \
		-o ${outdir}/${sample}_R1_filt.fastq.gz \
		-O ${outdir}/${sample}_R2_filt.fastq.gz \
		-h "$fastp_html"
	rm fastp.json
	fread=${outdir}/${sample}_R1_filt.fastq.gz
	rread=${outdir}/${sample}_R2_filt.fastq.gz
fi

## start report output
touch $report_file
printf "Reference Mapping Report\n" > $report_file
printf "==================================\n\n" >> $report_file
printf "## ${sample}\n\n" >> $report_file
printf "Forward reads file: ${fread}\n\n" >> $report_file
printf "Reverse reads file: ${rread}\n\n" >> $report_file
printf "* * * *\n\n" >> $report_file

############ Start loop to process each genome #################
IFS=',' 
for g in ${genomes}; do 
	ext=$(basename "${g}" | cut -d"." -f2)
	if [ "$ext" != "fasta" ]; then 			# if no file extension given
		refname="${g}"
		genome=~/refseq/"${g}".fasta
	else									# if full file/path is given
		refname=$(basename "${g}" | cut -d"." -f1)   # store name of reference file, w/o extension
		genome="${g}"
	fi
	## ensure the genome file is not in DOS format
	dos2unix ${genome}
	
	## establish genome-specific files
	fastp_html=${outdir}/stats/${sample}_fastp.html
	samfile=${outdir}/bam/${sample}_${refname}.sam
	bamfile=${outdir}/bam/${sample}_${refname}.bam
	sorted_bam=${outdir}/bam/${sample}_${refname}_sorted.bam
	flagstat=${outdir}/stats/${sample}_${refname}_mapping-stats.txt
	readdepth=${outdir}/stats/${sample}_${refname}_read_depth.txt
	plot=${outdir}/img/${sample}_${refname}_read_depth.jpeg
	raw_bcf=${outdir}/vcf/${sample}_${refname}.bcf
	variants=${outdir}/vcf/${sample}_${refname}_unfiltered.vcf
	final_variants=${outdir}/vcf/${sample}_${refname}.vcf
	norm_bcf=${outdir}/vcf/${sample}_${refname}_norm.bcf
	filt_bcf=${outdir}/vcf/${sample}_${refname}_norm_flt.bcf
	consensus=${outdir}/${sample}_${refname}_consensus_untrimmed.fasta
	final_consensus=${outdir}/${sample}_${refname}_consensus.fasta
	aln=${outdir}/stats/${sample}_${refname}_alignment.fasta
	misc_stats=${outdir}/stats/${sample}_${refname}_misc-stats.txt
	
	############################## ALIGNMENT ###################################

	### Align reads with bwa
	printf "\n*** Aligning reads with BWA MEM ***\n"
	bwa mem \
		-R '@RG\tID:'${sample}'\tSM:'${sample}'\tLB:1' \
		${genome} \
		"${fread}" \
		"${rread}" \
		> $samfile \
		2> ${outdir}/logs/bwa_stderr.log
	# if index files are missing, automatically create them
	if grep -q "fail to locate the index files" ${outdir}/logs/bwa_stderr.log; 
	then
		## create index files
		printf "\n*** BWA - creating index files ***\n"
		bwa index "${genome}"
		bwa faidx "${genome}"
		# this index file is useful for running the GATK variant caller
		java -Dpicard.useLegacyParser=false -jar picard.jar CreateSequenceDictionary \
			-REFERENCE "${genome}"
		### Align reads with bwa
		bwa mem \
			-R '@RG\tID:'${sample}'\tSM:'${sample}'\tLB:1' \
			${genome} \
			"${fread}" \
			"${rread}" \
			> $samfile \
			2> ${outdir}/logs/bwa_stderr.log
	fi

	############################## REFORMATTING ###################################

	## converts .sam to .bam
	printf "\n*** Samtools - View ***\n"
	samtools view -bS $samfile -o $bamfile

	## produces a sorted .bam file that can then be used by other programs (like bcftools)
	printf "\n*** Samtools - Sort ***\n"
	samtools sort $bamfile -o $sorted_bam 

	## produces a BAM index file (.bai) that could be useful for variant calling, etc.
	printf "\n*** Samtools - Index ***\n"
	samtools index $sorted_bam 

	## Print stats on how well the alignment worked
	printf "\n*** Samtools - Flagstat ***\n"
	samtools flagstat $bamfile > $flagstat

	## Print read depth information
	printf "\n*** Samtools - Read Depth ***\n"
	samtools depth $sorted_bam > $readdepth
	
	if [ -s $readdepth ] # procede only if reads mapped to the reference sequence
	then
		## Graph read depth 
		printf "\n*** RScript - Read Depth Graph ***\n"
		cd "${scriptdir}"/bin
		reflen=$(expr $(sed '1d' ${genome} | wc -m) - $(sed '1d' ${genome} | wc -l))
		printf "\n Reflen = $reflen\n\n"
		chmod a+x read_depth_graphs.R
		./read_depth_graphs.R $readdepth $plot $reflen
		
		############################## VARIANT CALLING ################################

		## produces read coverage info file for downstream variant calling
		bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam

		printf "\n*** Call Variants with BCFtools ***\n"

		## call variants
		bcftools call --ploidy 1 -mv -o $variants $raw_bcf

		## filter VCF (variant call format) file based on quality scores 
		vcfutils.pl varFilter $variants > $final_variants

		############################## MAKE CONSENSUS SEQUENCE ################################

		## normalize indels
		printf "\n*** BCFtools - norm ***\n"
		bcftools norm -f ${genome} $final_variants -Ob -o $norm_bcf

		## filter adjacent indels within 5bp
		printf "\n*** BCFtools - filter ***\n"
		bcftools filter \
			--IndelGap 5 $norm_bcf \
			-Ob -o $filt_bcf

		## index normalized, filtered variants for use in making consensus
		printf "\n*** BCFtools - index ***\n"
		bcftools index $filt_bcf

		## generate consensus file using reference file and variant call file (vcf/bcf)
		printf "\n*** BCFtools - consensus ***\n"
		cat ${genome} | bcftools consensus $filt_bcf > $consensus 2> $bcflog

		#"\n*** Calculating Nucleotide Alignment Identity ***\n"	
		## generate the truncated version of this file, based on the first and last 
		## positions listed in the read depth file (I made a python program to do this)
		printf "\n*** Python - trimming consensus ***\n"
		cd ${scriptdir}/bin
		chmod a+x assemblyRange.py
		./assemblyRange.py $readdepth $consensus $final_consensus

		## generate a reference sequence matching the lengths of the truncated file (for id calc)
		./assemblyRange.py $readdepth ${genome} ${outdir}/ref_trunc.fasta

		## produce a pairwise alignment file using muscle
		printf "\n*** Creating Pairwise Alignment File ***\n"
		awk '{print}' ${genome} $final_consensus > ${outdir}/muscle_input.fasta
		muscle -in ${outdir}/muscle_input.fasta -out $aln &> ${outdir}/logs/muscle_output.log
			
		## print percent identity of alignment
		printf "\n*** Alignment Info ***\n"
		chmod a+x percent_id.py
		ntnum=$(wc -l ${readdepth} | awk '{ print $1 }')
		./percent_id.py $final_consensus ${outdir}/ref_trunc.fasta ${genome} ${ntnum} > $misc_stats

		#################################### STORE STATS ###################################

		# calculate and store values
		nummap=$(head -n 5 ${flagstat} | tail -n 1 | sed -n "s/^\(\w\+\).\+/\1/p")
		pmap=$(head -n 5 ${flagstat} | tail -n 1 | sed -n "s/^.\+(\(.\+\)%.\+/\1/p")
		pid=$(head -n 2 $misc_stats | tail -n 1)
		refmap=$(head -n 4 $misc_stats | tail -n 1)
		conslen=$(head -n 6 $misc_stats | tail -n 1)
		nvar=$(tail -n 1 ${bcflog} | sed -n "s/^Applied\s\(.\+\)\svariants/\1/p")

		java_v=$(java --version | head -n 3 | tail -n 1 )
		python_v=$(python3 --version)
		fastp_v=$(fastp --version 2>&1 >/dev/null | head)
		rscript_v=$(Rscript --version 2>&1 >/dev/null | head)
		bwa_v=$(bwa 2>&1 >/dev/null | head -n 3 | tail -n 1)
		samtools_v=$(samtools --version | head -n 1)
		muscle_v=$(muscle -version | head -n 1)

		################################ APPEND MARKDOWN FILE #####################################

		printf "### Reads mapped to ${refname}\n\n" >> $report_file
		printf "Reference genome file: ${genome}\n\n" >> $report_file

		printf "### Summary Table for ${refname} Mapping\n\n" >> $report_file
		printf "| Metric                                           | Value               |\n" >> $report_file
		printf "| :----------------------------------------------- | :------------------ |\n" >> $report_file
		printf "| Percent of reads that mapped to reference:       | $pmap %%           |\n" >> $report_file
		printf "| Number of reads that mapped to reference:        | $nummap            |\n" >> $report_file
		printf "| Percent identity (identical nt matches):         | $pid %%            |\n" >> $report_file
		printf "| Percent of reference genome with >0 read depth:  | $refmap %%         |\n" >> $report_file
		printf "| Consensus sequence length (bp):                  | $conslen           |\n" >> $report_file
		printf "| Number of variants in consensus:                 | $nvar             |\n" >> $report_file

		# read coverage image
		printf "![](img/${sample}_${refname}_read_depth.jpeg){ width=40%% }\n\n" >> $report_file
		printf "* * * *\n\n" >> $report_file

		############################## FILE CLEANUP ###################################

		rm ${outdir}/ref_trunc.fasta 
		rm ${outdir}/muscle_input.fasta 
		rm $consensus
		rm $samfile
		rm $bamfile
		rm $raw_bcf
		rm $norm_bcf
		rm $filt_bcf
		rm $misc_stats
		rm ${outdir}/logs/*.log
		rm ${filt_bcf}.csi
		
	else
		printf "### Reads mapped to ${refname}\n\n" >> $report_file
		printf "Reference genome file: ${genome}\n\n" >> $report_file
		printf "*No reads mapped to reference genome.*\n\n" >> $report_file
		printf "* * * *\n\n" >> $report_file
	fi

done

# print all the versions of the programs used in this script
printf '%s\n\n' "### Dependencies:" >> $report_file
printf '%s\n\n' "- ${java_v}" >> $report_file
printf '%s\n\n' "- ${python_v}" >> $report_file
printf '%s\n\n' "- ${rscript_v}" >> $report_file
printf '%s\n\n' "- BWA ${bwa_v}" >> $report_file
printf '%s\n\n' "- ${samtools_v}" >> $report_file
printf '%s\n\n' "- ${muscle_v}" >> $report_file
if [[ "$fastp" = "1" ]]; then 
	printf '%s\n\n' "- ${fastp_v}" >> $report_file
	printf '\t\t%s\n\n' "$(tail -n 1 ${fastp_html} | sed -n "s/^<div id='footer'> <p>\(.\+\)-i.\+/\1/p")" >> $report_file
fi

printf '\n%s\n\n' "*Script by Angela Crabtree*" >> $report_file

# convert markdown to html using pandoc
pandoc --from markdown --to html $report_file > $report_html

# remove any additional files
rm -r ${outdir}/logs
rm -f ${outdir}/*_filt.fastq.gz

printf "\n\tAll done!\n\n"