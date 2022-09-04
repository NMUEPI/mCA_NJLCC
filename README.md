# Mosaic-chromosomal-alteration

This repository contains custom scripts used for raw IDAT files from the Illumina GSA processing and mosaic chromosomal alterations (mCAs) detecting. 

## **Reference**
### **Download** 
You can find the required GRCh37 resources at https://github.com/freeseek/mocha#download-resources-for-grch37 or you can download them as follows

#### Human genome reference

    wget -O- ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz | \
    gzip -d > /public/home/naq/ref/GRCh37/human_g1k_v37.fasta
    samtools faidx /public/home/naq/ref/GRCh37/human_g1k_v37.fasta

#### Genetic map

    wget -P /public/home/naq/ref/GRCh37 https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/genetic_map_hg19_withX.txt.gz

#### 1000 Genomes project phase 3

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

#### List of common germline duplications and deletions

    wget -P /public/home/naq/ref/GRCh37 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz{,.tbi}
    bcftools query -i 'AC>1 && END-POS+1>10000 && SVTYPE!="INDEL" && (SVTYPE=="CNV" || SVTYPE=="DEL" || SVTYPE=="DUP")' \
     -f "%CHROM\t%POS0\t%END\t%SVTYPE\n" /public/home/naq/ref/GRCh37/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz > /public/home/naq/ref/GRCh37/cnps.bed

#### Minimal divergence intervals from segmental duplications (make sure your bedtools version is 2.27 or newer)

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

#### 1000 Genomes project phase 3 imputation panel for beagle5 and impute5

    cd /public/home/naq/ref/GRCh37
    pfx="ALL.chr"
    sfx=".phase3_integrated.20130502.genotypes"
    for chr in {{1..22},X}; do imp5Converter --h $pfx$chr$sfx.bcf --o $pfx$chr$sfx --r $chr; done
    wget -O bref3.jar http://faculty.washington.edu/browning/beagle/bref3.28Jun21.220.jar
    for chr in {1..22}; do bcftools view --no-version $pfx$chr$sfx.bcf | java -jar bref3.jar > $pfx$chr$sfx.bref3; done
    chr=X; bcftools +fixploidy --no-version $pfx$chr$sfx.bcf | \
    sed 's/0\/0/0|0/g;s/1\/1/1|1/g' | java -jar ~/bin/bref3.jar > $pfx$chr$sfx.bref3

#### Download cytoband file

    wget -P /public/home/naq/ref/GRCh37 http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz

### **Data Processing** 

This step is to create every chromosome files for phasing and imputation from reference genome.

    source ${config_path}/config.sh

    for chr in {1..22} X; do
      zcat ${map} | sed 's/^23/X/' | awk -v chr=${chr} '$1==chr {print $2,$3,$4}' > ${phaseref_path}/genetic_map.chr${chr}.txt
    done

    for chr in {1..22} X; do
     zcat ${map} | sed 's/^23/X/' | awk -v chr=${chr} '$1==chr {print chr,".",$4,$2}' > ${imputeref_path}/genetic_map.chr${chr}.txt
    done

## **Overview**

This repository has two directories:

Part1 is to convert the raw intensity data (IDAT) files from the Illumina Infinium Global Screening Array into VCF files. 

    ${tool_path}/iaap-cli \
    gencall \
    ${bpm_manifest_file} \
    ${egt_cluster_file} \
    ${path_to_output_folder} \
     --idat-folder $path_to_idat_folder \
     --output-gtc \
     --gender-estimate-call-rate-threshold -0.1

    bcftools +gtc2vcf \
     --no-version -Ou \
     --bpm ${bpm_manifest_file} \
     --csv ${csv_manifest_file} \
     --egt ${egt_cluster_file} \
     --gtcs ${path_to_gtc_folder} \
     --fasta-ref ${ref} \
     --extra ${vcf_path}/${out_prefix}.tsv | \
    bcftools sort -Ou -T ${output_path}/bcftools-sort.XXXXXX | \
    bcftools norm --no-version -o ${vcf_path}/${out_prefix}.vcf -c x -f ${ref} && \  
    bgzip ${vcf_path}/${out_prefix}.vcf && \
    bcftools index -f ${vcf_path}/${out_prefix}.vcf.gz 

Part2 is to detect mosaic chromosomal alterations (mCAs). 

### **Preprocessing**

    vcftools --gzvcf ${vcf_path}/${out_prefix}.unphased.vcf.gz --positions ${output_path}/${out_prefix}_updateQC.list --recode --recode-INFO-all --out ${vcf_path}/${pfx}.unphased.filter.vcf
    mv ${vcf_path}/${out_prefix}.unphased.filter.vcf.recode.vcf ${vcf_path}/${out_prefix}.unphased.filter.vcf
    rm ${vcf_path}/${out_prefix}.unphased.filter.vcf.gz
    bgzip ${vcf_path}/${out_prefix}.unphased.filter.vcf 
    bcftools index -f ${vcf_path}/${out_prefix}.unphased.filter.vcf.gz

    echo '##INFO=<ID=JK,Number=1,Type=Float,Description="Jukes Cantor">' | \
    bcftools annotate --no-version -Ou -a ${dup} -c CHROM,FROM,TO,JK -h /dev/stdin ${vcf_path}/${out_prefix}.unphased.filter.vcf.gz | \
    bcftools view --no-version -Ou -G | \
    bcftools annotate --no-version -Ob -o ${vcf_path}/${out_prefix}.exclude.bcf \
     -i 'FILTER!="." && FILTER!="PASS" || INFO/JK<.02' \
     -x ^INFO/JK && \
    bcftools index -f ${vcf_path}/${out_prefix}.exclude.bcf


    bcftools isec --no-version -Ou --complement --exclude "N_ALT>1" --write 1 ${vcf_path}/${out_prefix}.unphased.filter.vcf.gz ${vcf_path}/${out_prefix}.exclude.bcf | \
    bcftools annotate --no-version -Ou --remove ID,QUAL,INFO,^FMT/GT  | \
    bcftools +scatter --no-version -Ob --output ${vcf_path} --scatter $(echo {{1..22},X} | tr ' ' ',') --prefix ${out_prefix}.chr


### **Phasing**

    bcftools index --force ${vcf_path}/${out_prefix}.chr${CHR}.bcf

    shapeit4 \
     --thread ${phase_threads} \
     --input ${vcf_path}/${out_prefix}.chr${CHR}.bcf \
     --reference ${panel_pfx}${CHR}${panel_sfx}.bcf \
     --map ${vcf_path}/genetic_map.chr${CHR}.txt \
     --region CHR \
     --output ${phase_path}/${out_prefix}.chr${CHR}.pgt.bcf

    bcftools concat --no-version -Ob -o ${phase_path}/${out_prefix}.pgt.bcf ${phase_path}/${out_prefix}.chr{{1..22},X}.pgt.bcf && \
    bcftools index -f ${phase_path}/${out_prefix}.pgt.bcf

    bcftools annotate --no-version -Ob -o ${phase_path}/${out_prefix}.bcf --annotations ${phase_path}/${out_prefix}.pgt.bcf --columns -FMT/GT ${vcf_path}/${out_prefix}.unphased.filter.vcf.gz && \
    bcftools index -f ${phase_path}/${out_prefix}.bcf

### **Imputation**

    impute5 \
     --h ${panel_pfx}${CHR}${panel_sfx}.bcf \
     --m ${imputeref_path}/genetic_map.chr${CHR}.txt \
     --g ${phase_path}/${out_prefix}.pgt.bcf \
     --r ${CHR} \
     --o ${impute_path}/${out_prefix}.chr${CHR}.imp.bcf \
     --l ${impute_path}/${out_prefix}.chr${CHR}.log \
     --threads ${impute_threads}

    bcftools concat --no-version -Ob -o ${impute_path}/${out_prefix}.imp.bcf ${impute_path}/${out_prefix}.chr{1..22}.imp.bcf && \
    bcftools index -f ${impute_path}/${out_prefix}.imp.bcf
	
    bcftools concat --no-version -o ${impute_path}/${out_prefix}.imp.vcf.gz ${impute_path}/${out_prefix}.chr{1..22}.imp.bcf && \
    bcftools index -f ${impute_path}/${out_prefix}.imp.vcf.gz

### **mCAs**

    bcftools +mocha \
     --rules ${assembly} ${pase_path}/${out_prefix}.bcf \
     --sex ${sexstat} \
     --call-rate ${callratestat} \
     --no-version \
     --output-type b \
     --output ${mocha_path}/${out_prefix}.bdev.bcf \
     --mosaic-calls ${mocha_path}/${out_prefix}.calls.tsv \
     --genome-stats ${mocha_path}/${out_prefix}.stats.tsv \
     --ucsc-bed ${mocha_path}/${out_prefix}.ucsc.bed \
     && \
    bcftools index -f ${mocha_path}/${out_prefix}.bdev.bcf

If you want to 