# **Detect mosaic chromosomal alterations (mCAs)**

## **1. Preprocessing**

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


## **2. Phasing**

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

## **3. Imputation**

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

## **4. MoChA Detection**

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
