#!/bin/bash
cat << "EOF"
"
O---o               -/ossssssssso+/.               O---o
 O-o            `+syo+:::::::::::/+oyo+`            O-o
  O           -sy+/////:://///////////+ys:           O
 o-O        -hmshhhhy+++sssosyyy+++hhhhhyNy-        o-O
o---O      sy/hhMMMMMd: ./oyyho``oNMMMMMom/hs      o---O
O---o    `ho::hhMMMMMMMd:/yohm-oNMMMMMMMom::oy     O---o
 O-o     hs///hyMMMMMMMMMdsysyNMMMMMMMMmom///yy     O-o
  O     /h:::-y:.sMMMMMh//odho+yNMMMMm/ om--::d-     O
 o-O    d+////h:  .sNMN:...+...:sMMm/   om////oy    o-O
o---O   h/:/:/y:    `sm:...+.-./sm/     om::/:/d   o---O
O---o   d:::::h:    .ym:.-.+`.-/yN+`    om::::/d   O---o
 O-o    d+////y:  .yMMN+-:-s::/omMMNo`  sm////ss    O-o
  O     +y--::y:-hMMMMMMMMMsmMMMMMMMMN+`sd-::-d.     O
 o-O     ds//+ydMMMMMMMMMd: `oNMMMMMMMMNsm///ys     o-O
o---O    `d+::yNMMMMMMMhmhsdhomhNMMMMMMMsm-:sy     o---O  
O---o     `yy/yNMMMMMh- :+mmmho:`+NMMMMMsm+ho      O---o
 O-o        :ymyMMMh-   omMdMNd:  `+NMNhhmy.        O-o
  O            :hdhyo+//:ymd+ymds++osyhhNy-          O
 o-O            .+sso+/+oos+hyoo++osys/`            o-O
o---O               ./ssssymdsssso/.               o---O
O---o                                              O---o
 O-o                 Schneiders Lab                 O-o
  O                 RNA-Seq Pipeline                 O
 o-O                                                o-O
o---O                  TDS - 2018                  o---O
EOF

#Accession Number Input
read -p "Was each sample was run over two lanes? [Y/n] " multilane
read -p "Desired number of threads: " threads
read -p "Desired kmer length: " kmerlen
read -p "Number of bootstrapping iterations (at least 30 if performing DE analysis in Sleuth): " bootstraps
    printf "

Enter sample accession number followed by 'Enter'

NB. If samples were run over two lanes, enter the accession number for Lane 1, then Lane 2.

When finished press 'Ctrl+D':

"

    while read line
    do
        ana+=( $line )
    done
    printf -- '\n%s ' "The analysis will include following accession numbers and downloads will begin:

${ana[@]}"

    read -r -p "

Are You Sure? [Y/n] " input
    
    case $input in
        [yY][eE][sS]|[yY])
    echo "Yes"
    ;;
    
        [nN][oO]|[nN])
    echo "No"
    exit 1
    ;;
    esac

# Make New Directory
    mkdir kallisto_$(date +%d.%m.%Y_%H:%M:%S) && cd "$_"

# Downloads
    #Accession
        for i in "${ana[@]}"
        do
            curl \
            ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${i::-3}/${i}/${i}_1.fastq.gz \
            | gunzip  >> ${i}_1.fastq
        done
        for i in "${ana[@]}"
        do  
            curl \
            ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${i::-3}/${i}/${i}_2.fastq.gz \
            | gunzip  >> ${i}_2.fastq
        done
        wait
        printf "
All FASTQ Files Downloaded!
    "

# FASTQC
    printf "
Running fastQC
    "
	mkdir fastqc_results
	fastqc *.fastq -o fastqc_results
    wait
    google-chrome fastqc_results/*.html
    printf "
fastQC Complete!
    "

# ECL8 Transcriptome
    printf "
Downloading Transcriptome
    "
	curl \
	ftp://ftp.ensemblgenomes.org/pub/bacteria/release-41/fasta/bacteria_9_collection/klebsiella_pneumoniae_subsp_pneumoniae_ecl8/cdna/Klebsiella_pneumoniae_subsp_pneumoniae_ecl8.ASM31538v1.cdna.all.fa.gz \
    ftp://ftp.ensemblgenomes.org/pub/bacteria/release-41/fasta/bacteria_9_collection/klebsiella_pneumoniae_subsp_pneumoniae_ecl8/ncrna/Klebsiella_pneumoniae_subsp_pneumoniae_ecl8.ASM31538v1.ncrna.fa.gz \
    | gunzip | cat >> ecl8transcriptome.fa
    printf "
Transcriptome downloaded!
    "

#Gene List
    grep -e ">" ecl8transcriptome.fa > headers.fa
    sed 's/gene://g;s/description://g;s/gene_symbol://g;s/,//g' headers.fa > out.fa
    awk 'BEGIN{RS=">"}{print $4"\t"$7" "$8" "$9" "$10" "$11;}' out.fa | tail -n+2 > genelist.tsv
    sed -i 1i"Gene","Description" genelist.tsv
    mkdir quants

#Bulding kallisto Index
		printf "
Building kallisto Index
"
	kallisto index \
	--make-unique \
    --kmer-size=${kmerlen} \
	-i kallisto_index \
	ecl8transcriptome.fa
		wait
        printf "
kallisto Index Built
"
if [[ "$multilane" == "Y" ]]
    then
    #kallisto One Sample Two-FASTQ
        printf "
Starting quasi-mapping with kallisto
        "
        for (( i=0; i<${#ana[@]} ; i+=2 ))
        do
            kallisto quant \
            -t ${threads} \
			-i kallisto_index \
            -o quants/${ana[i]}_${ana[i+1]} \
            -b ${bootstraps} \
            ${ana[i]}_1.fastq ${ana[i]}_2.fastq \
            ${ana[i+1]}_1.fastq ${ana[i+1]}_2.fastq
        wait
        done
    else
    #kallisto One Sample One FASTQ
	    for i in "${ana[@]}"
	    do 
            kallisto quant \
            -t ${threads} \
			-i kallisto_index \
            -o quants/${i} \
            -b ${bootstraps} \
            ${i}_1.fastq ${i}_2.fastq
        wait
        done
fi
printf "
kallisto complete!
        "
sleep 1
#Merge the abundance files
    cd quants
    for D in *; do [ -d "${D}" ] && \
	cut -f4 -d$'\t' ${D}/abundance.tsv | tail -n +2 > ${D}/${D}_headless.tsv
	echo -e "${D}" | cat - ${D}/${D}_headless.tsv > ${D}.tsv&
	done
    cd ../
    mv genelist.tsv $PWD/quants/genelist.tsv
    cd quants
	paste \
    *.tsv \
	> kallisto_counts.tsv
    cd ../
    mv $PWD/quants/kallisto_counts.tsv $PWD/kallisto_counts.tsv

# File Cleanup
    cd quants
    for D in *; do [ -d "${D}" ] && \
	do
    rm ${D}/${D}_headless.tsv&
    rm ${D}.tsv
    done
    cd ../
    tar -czvf FASTQfiles.tar.gz *.fastq
    wait
    rm genelist.tsv&
    rm *.fastq&
    rm kallisto_index&
    rm fastqc_results/*.zip&
    wait
    printf "
    Directory Cleaned!
Analysis Complete! Check, then upload kallisto_counts.tsv here: http://degust.erc.monash.edu/upload"