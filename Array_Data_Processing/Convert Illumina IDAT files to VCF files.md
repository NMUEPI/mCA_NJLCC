# **Convert Illumina IDAT files to VCF files**

## **1.  Convert IDAT files to GTC files**

> Convert raw Intensity Data files (.idat) to Genotype Call files (.gtc) with GenomeStudio software.


## **2.  Convert GTC files to VCF files**

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
