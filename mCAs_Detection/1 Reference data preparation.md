# **Reference Data Preparation**

## **1. Download** 
You can find the required GRCh37 resources at https://github.com/freeseek/mocha#download-resources-for-grch37 or you can download them as follows

### Human genome reference

    wget -O- ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz | \
    gzip -d > /public/home/naq/ref/GRCh37/human_g1k_v37.fasta
    samtools faidx /public/home/naq/ref/GRCh37/human_g1k_v37.fasta

### Genetic map

    wget -P /public/home/naq/ref/GRCh37 https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/genetic_map_hg19_withX.txt.gz

### 1000 Genomes project phase 3

    cd /public/home/naq/ref/GRCh37
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{{1..22}.phase3_shapeit2_mvncall_integrated_v5b,X.phase3_shapeit2_mvncall_integrated_v1c,Y.phase3_integrated_v2b}.20130502.genotypes.vcf.gz{,.tbi}
    for chr in {1..22} X Y; do
      bcftools view --no-version -Ou -c 2 ALL.chr${chr}.phase3*integrated_v[125][bc].20130502.genotypes.vcf.gz | \
      bcftools annotate --no-version -Ou -x ID,QUAL,FILTER,INFO,INFO/END,^FMT/GT | \
      bcftools norm --no-version -Ou -m -any | \
      bcftools norm --no-version -Ou -d none -f $HOME/GRCh37/human_g1k_v37.fasta | \
      bcftools sort -Ob -T ./bcftools-sort.XXXXXX | \
      tee ALL.chr${chr}.phase3_integrated.20130502.genotypes.bcf | \
      bcftools index --force --output ALL.chr${chr}.phase3_integrated.20130502.genotypes.bcf.csi
    done

### List of common germline duplications and deletions

    wget -P /public/home/naq/ref/GRCh37 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz{,.tbi}
    bcftools query -i 'AC>1 && END-POS+1>10000 && SVTYPE!="INDEL" && (SVTYPE=="CNV" || SVTYPE=="DEL" || SVTYPE=="DUP")' \
     -f "%CHROM\t%POS0\t%END\t%SVTYPE\n" /public/home/naq/ref/GRCh37/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz > /public/home/naq/ref/GRCh37/cnps.bed

### Minimal divergence intervals from segmental duplications (make sure your bedtools version is 2.27 or newer)

    wget -O- http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz | gzip -d | awk '!($2=="chrX" && $8=="chrY" || $2=="chrY" && $8=="chrX") {print $2"\t"$3"\t"$4"\t"$30}' > genomicSuperDups.bed

    awk '{print $1,$2; print $1,$3}' genomicSuperDups.bed | \
     sort -k1,1 -k2,2n | uniq | \
     awk 'chrom==$1 {print chrom"\t"pos"\t"$2} {chrom=$1; pos=$2}' | \
     bedtools intersect -a genomicSuperDups.bed -b - | \
     bedtools sort | \
     bedtools groupby -c 4 -o min | \
     awk 'BEGIN {i=0; s[0]="+"; s[1]="-"} {if ($4!=x) i=(i+1)%2; x=$4; print $0"\t0\t"s[i]}' | \
     bedtools merge -s -c 4 -o distinct | \
     sed 's/^chr//' | grep -v gl | bgzip > /public/home/naq/ref/GRCh37/segdups.bed.gz && \
     tabix -f -p bed /public/home/naq/ref/GRCh37/segdups.bed.gz

### 1000 Genomes project phase 3 imputation panel for beagle5 and impute5

    cd /public/home/naq/ref/GRCh37
    pfx="ALL.chr"
    sfx=".phase3_integrated.20130502.genotypes"
    for chr in {{1..22},X}; do imp5Converter --h $pfx$chr$sfx.bcf --o $pfx$chr$sfx --r $chr; done
    wget -O bref3.jar http://faculty.washington.edu/browning/beagle/bref3.28Jun21.220.jar
    for chr in {1..22}; do bcftools view --no-version $pfx$chr$sfx.bcf | java -jar bref3.jar > $pfx$chr$sfx.bref3; done
    chr=X; bcftools +fixploidy --no-version $pfx$chr$sfx.bcf | \
    sed 's/0\/0/0|0/g;s/1\/1/1|1/g' | java -jar ~/bin/bref3.jar > $pfx$chr$sfx.bref3

### Download cytoband file

    wget -P /public/home/naq/ref/GRCh37 http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz

## **2. Data Process** 

This step is to create every chromosome files for phasing and imputation from reference genome.

    source ${config_path}/config.sh

    for chr in {1..22} X; do
      zcat ${map} | sed 's/^23/X/' | awk -v chr=${CHR} '$1==chr {print $2,$3,$4}' > ${phaseref_path}/genetic_map.chr${CHR}.txt
    done

    for chr in {1..22} X; do
     zcat ${map} | sed 's/^23/X/' | awk -v chr=${CHR} '$1==chr {print chr,".",$4,$2}' > ${imputeref_path}/genetic_map.chr${CHR}.txt
    done
