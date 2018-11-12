#!/bin/bash -eu
# -e: Exit immediately if a command exits with a non-zero status.
# -u: Treat unset variables as an error when substituting.
#Variables

	#File Types
		fz='.fastq.gz'

	#Directories
		maindir='/home/tom/Documents/University/Bioinformatics/CAIN/RNA-Seq/FromFASTQ'
		sal='/home/tom/Documents/University/Bioinformatics/CAIN/RNA-Seq/FromFASTQ/Salmon'
		ftpdir='ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR576'
		D0='/home/tom/Documents/University/Bioinformatics/CAIN/RNA-Seq/FromFASTQ/D0'
		D3='/home/tom/Documents/University/Bioinformatics/CAIN/RNA-Seq/FromFASTQ/D3'
		D5='/home/tom/Documents/University/Bioinformatics/CAIN/RNA-Seq/FromFASTQ/D5'

	#New file name array
		nfn=('Rep1L1'
			'Rep1L2'
			'Rep2L1'
			'Rep2L2')

	#Accession Number Array
		arr=('ERR576157_1'
			'ERR576157_2'
			'ERR576173_1'
			'ERR576173_2'
			'ERR576125_1'
			'ERR576125_2'
			'ERR576141_1'
			'ERR576141_2'
			'ERR576159_1'
			'ERR576159_2'
			'ERR576175_1'
			'ERR576175_2'
			'ERR576127_1'
			'ERR576127_2'
			'ERR576143_1'
			'ERR576143_2'
			'ERR576161_1'
			'ERR576161_2'
			'ERR576177_1'
			'ERR576177_2'
			'ERR576129_1'
			'ERR576129_2'
			'ERR576145_1'
			'ERR576145_2')

#Salmon

	# Build Salmon Index
		cd $sal
		salmon --no-version-check \
		index \
		-t $maindir/ecl8transcriptome.fa \
		-i $sal/ecl8_index
	#Quantifying the FASTQs
		cd $sal
		mkdir quants
		#D0
		for samp in "${nfn[@]}"
		do
			espeak -s 140 "Processing sample D0 ${samp}"
			salmon --no-version-check \
			quant -i ecl8_index -l A \
         		-1  $D0/FASTQ/D0${samp}R1$fz \
         		-2  $D0/FASTQ/D0${samp}R2$fz \
         		-p 8 -o quants/D0${samp}_quant
		done
		echo "Day 0 Complete"
		#D3
		for samp in "${nfn[@]}"
		do
			espeak -s 140 "Processing sample D3 ${samp}"
			salmon --no-version-check \
			quant -i ecl8_index -l A \
         		-1  $D3/FASTQ/D3${samp}R1$fz \
         		-2  $D3/FASTQ/D3${samp}R2$fz \
         		-p 8 -o quants/D3${samp}_quant
		done
		echo "Day 3 Complete"
		#D5
		for samp in "${nfn[@]}"
		do
			espeak -s 140 "Processing sample D5 ${samp}"
			salmon --no-version-check \
			quant -i ecl8_index -l A \
         		-1  $D5/FASTQ/D5${samp}R1$fz \
         		-2  $D5/FASTQ/D5${samp}R2$fz \
         		-p 8 -o quants/D5${samp}_quant
		done
		echo "All quants complete!!!"