##script for detecting contamination in UCE samples using heterozygosity measures##

#step 1: create reference dictionary by taking the longest contigs in my UCE loci
#note that ALL outgroups have been removed--this is just Pogonomyrmex

for file in /home/lcg65/Desktop/PogoUCE/2024_for_decontam/*.fasta
do
	name=`basename $file .fasta`
	cat $file | bioawk -c fastx '{ print length($seq), $name }' | sort -k1,1rn | head -1 > $name.longest.txt
	contig_name=`awk '{print $2}' $name.longest.txt`
	echo "$contig_name"
	samtools faidx $file "$contig_name" > longest_contig_per_locus/$name.longest_contig.fasta
done

#it's mostly pogonomyrmex abdominalis...but maybe that's okay

#add file name into header
cd longest_contig_per_locus
for f in *.fasta; do sed -i "s/^>/>${f%.longest_contig.fasta} /g" "${f}"; done

#concatenate the longest contigs into a reference genome
cat *.fasta > pogo_decontam_reference.fasta

#step 2: bwa-samtools loop with cleaned fastq & reference genome

#index constructed reference genome
/programs/bwa-mem2-2.2.1/bwa-mem2 index /home/lcg65/Desktop/PogoUCE/pogo_decontam_reference/pogo_decontam_reference.fasta

#loop cleaned fastq through bwa & samtools--get sorted bam
for file in /home/lcg65/Desktop/PogoUCE/clean-fastq-all/clean-fastq-final/*READ1.fastq.gz;
do
  name=$(basename "$file" -READ1.fastq.gz)
  echo "Assembling plastome for $name"
  forward=$name"-READ1.fastq.gz";
  reverse=$name"-READ2.fastq.gz";
  bam=/local/workdir/lcg65/"$name".bam;
  /programs/bwa-mem2-2.2.1/bwa-mem2 mem -t 20 /home/lcg65/Desktop/PogoUCE/pogo_decontam_reference/pogo_decontam_reference.fasta /home/lcg65/Desktop/PogoUCE/clean-fastq-all/clean-fastq-final/"$forward" /home/lcg65/Desktop/PogoUCE/clean-fastq-all/clean-fastq-final/"$reverse" | samtools view -b -h  | samtools sort --threads 20 -o "$bam" - ;
done

#step 3: make a dictionary with reference fasta
#make sure update java is installed with conda
conda install bioconda::java-jdk

#create dictionary file for reference genome
/programs/gatk4/gatk CreateSequenceDictionary R=/home/lcg65/Desktop/PogoUCE/pogo_decontam_reference/pogo_decontam_reference.fasta O=/home/lcg65/Desktop/PogoUCE/pogo_decontam_reference/pogo_decontam_reference.dict

#create index file with faidx for reference genome
samtools faidx /home/lcg65/Desktop/PogoUCE/pogo_decontam_reference/pogo_decontam_reference.fasta

#step 4:angst loop
for file in /local/workdir/lcg65/*bam
do
	name=`basename $file .bam`
	echo "Running angsd on $name to calculate genome wide heterozygostiy"
	/programs/angsd-0.940/angsd -i $file -anc /home/lcg65/Desktop/PogoUCE/pogo_decontam_reference/pogo_decontam_reference.fasta -dosaf 1 -gl 1 -minMapQ 30 -minQ 20
	/programs/angsd-0.940/misc/realSFS angsdput.saf.idx > $name.est.ml
	rm angsdput.saf.pos.gz
	rm angsdput.saf.idx
	rm angsdput.saf.gz
	rm angsdput.arg
done

for file in ./*.est.ml
do
	name=`basename $file .est.ml`
	awk '{print FILENAME, $2/($1 + $2 +$3)}' $file >> heteroszygosity_calculation.txt
done
#######################################################################################
#TRYING THIS AGAIN BUT WITH P. CALIFORNICUS REFERENCE GENOME#

#index P. californicus reference genome
/programs/bwa-mem2-2.2.1/bwa-mem2 index /local/workdir/lcg65/p_cali_genome/GCA_024349325.1_Pcal.v2_genomic.fna

#loop cleaned fastq through bwa & samtools--get sorted bam
for file in /home/lcg65/Desktop/PogoUCE/clean-fastq-all/clean-fastq-final/*READ1.fastq.gz;
do
  name=$(basename "$file" -READ1.fastq.gz)
  echo "Assembling plastome for $name"
  forward=$name"-READ1.fastq.gz";
  reverse=$name"-READ2.fastq.gz";
  bam=/local/workdir/lcg65/p_cali_bam_files/"$name".bam;
  /programs/bwa-mem2-2.2.1/bwa-mem2 mem -t 20 /local/workdir/lcg65/p_cali_genome/GCA_024349325.1_Pcal.v2_genomic.fna /home/lcg65/Desktop/PogoUCE/clean-fastq-all/clean-fastq-final/"$forward" /home/lcg65/Desktop/PogoUCE/clean-fastq-all/clean-fastq-final/"$reverse" | samtools view -b -h  | samtools sort --threads 20 -o "$bam" - ;
done

#step 3: make a dictionary with reference fasta
#make sure update java is installed with conda
conda install bioconda::java-jdk

#create dictionary file for reference genome
/programs/gatk4/gatk CreateSequenceDictionary R=/local/workdir/lcg65/p_cali_genome/GCA_024349325.1_Pcal.v2_genomic.fna O=/local/workdir/lcg65/p_cali_genome/GCA_024349325.1_Pcal.v2_genomic.dict

#create index file with faidx for reference genome
samtools faidx /local/workdir/lcg65/p_cali_genome/GCA_024349325.1_Pcal.v2_genomic.fna

#step 4:angst loop
for file in /local/workdir/lcg65/p_cali_bam_files/*bam
do
	name=`basename $file .bam`
	echo "Running angsd on $name to calculate genome wide heterozygostiy"
	/programs/angsd-0.940/angsd -i $file -anc /local/workdir/lcg65/p_cali_genome/GCA_024349325.1_Pcal.v2_genomic.fna -dosaf 1 -gl 1 -minMapQ 30 -minQ 20
	/programs/angsd-0.940/misc/realSFS angsdput.saf.idx > $name.est.ml
	rm angsdput.saf.pos.gz
	rm angsdput.saf.idx
	rm angsdput.saf.gz
	rm angsdput.arg
done

#step 5: quantify heterozygosity
for file in ./*.est.ml
do
	name=`basename $file .est.ml`
	awk '{print FILENAME, $2/($1 + $2 +$3)}' $file >> heteroszygosity_calculation.txt
done
