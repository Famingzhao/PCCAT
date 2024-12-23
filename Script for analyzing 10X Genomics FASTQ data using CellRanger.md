### Script for analyzing 10X Genomics FASTQ data using CellRanger

- Example data set: SRP304942

- Data set was published in https://www.sciencedirect.com/science/article/pii/S030228382200001X?via%3DihubPre-existing 

- Cheng Q, Butler W, Zhou Y, Zhang H, Tang L, Perkinson K, Chen X, Jiang XS, McCall SJ, Inman BA, Huang J. Pre-existing Castration-resistant Prostate Cancer-like Cells in Primary Prostate Cancer Promote Resistance to Hormonal Therapy. Eur Urol. 2022 May;81(5):446-455. doi: 10.1016/j.eururo.2021.12.039. Epub 2022 Jan 17. Erratum in: Eur Urol. 2023 Jun;83(6):e170-e171. PMID: 35058087; PMCID: PMC9018600.

  

#### Step1. FASTQ Download：

```shell
mkdir 1.rawdata
conda activate 10x
########################### bam files
cat >download_file
/vol1/run/SRR136/SRR13644607/possorted_genome_bam12.bam.1
/vol1/run/SRR136/SRR13644608/possorted_genome_bam11.bam.1
/vol1/run/SRR136/SRR13644609/possorted_genome_bam10.bam.1
/vol1/run/SRR136/SRR13644610/possorted_genome_bam9.bam.1
/vol1/run/SRR136/SRR13644611/possorted_genome_bam8.bam.1
/vol1/run/SRR136/SRR13644612/possorted_genome_bam6.bam.1
/vol1/run/SRR136/SRR13644613/possorted_genome_bam4.bam.1
/vol1/run/SRR136/SRR13644614/possorted_genome_bam3.bam.1
/vol1/run/SRR136/SRR13644615/possorted_genome_bam13.bam.1
/vol1/run/SRR136/SRR13644616/possorted_genome_bam2.bam.1
/vol1/run/SRR136/SRR13644617/possorted_genome_bam1.bam.1

nohup ascp -v -QT -l 300m -P33001 -k1 -i ~/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh --mode recv --host fasp.sra.ebi.ac.uk --user era-fasp --file-list download_file ./ >download.log 2>&1 &
```



#### Step2. Cellranger bamtofastq

```shell
##3.3 name list
cat >name.list.txt
SRR13644607 possorted_genome_bam12.bam.1
SRR13644608 possorted_genome_bam11.bam.1
SRR13644609 possorted_genome_bam10.bam.1
SRR13644610 possorted_genome_bam9.bam.1
SRR13644611 possorted_genome_bam8.bam.1
SRR13644612 possorted_genome_bam6.bam.1
SRR13644613 possorted_genome_bam4.bam.1
SRR13644614 possorted_genome_bam3.bam.1
SRR13644615 possorted_genome_bam13.bam.1
SRR13644616 possorted_genome_bam2.bam.1
SRR13644617 possorted_genome_bam1.bam.1


##3.4 create file
mkdir fastq_file

##3.5 shell
cat name.list.txt |while read id
do
srr=$(echo $id | cut -d " " -f 1)
bam=$(echo $id | cut -d " " -f 2)
echo "nohup cellranger bamtofastq --nthreads 10 --traceback ${bam} ./fastq_file/${srr} >${srr}.log 2>&1 &"
done > bamtofastq.sh

#run
bash bamtofastq.sh 
```



#### Step3.Cellranger count

```shell
## shell
cat name.list.txt | while read id
do
srr=$(echo $id | cut -d " " -f 1)
ref=~/software_install/10x_refernce/refdata-gex-GRCh38-2020-A
sample_name=${srr}_out_file
echo "cellranger count --id=$sample_name \
--transcriptome=$ref \
--fastqs=${srr} \
--sample=bamtofastq \
--nosecondary \
--localmem=20 \
--localcores=30"
done>file.sh

#run
nohup bash file.sh>Run_count.log 2>&1 &
```
