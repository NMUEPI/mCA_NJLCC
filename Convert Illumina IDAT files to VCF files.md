# **Convert Illumina IDAT files to VCF files**

**1.  Illumina GenCall**

    ${tool_path}/iaap-cli \
     gencall \
     ${bpm_manifest_file} \
     ${egt_cluster_file} \
     ${path_to_output_folder} \
      --idat-folder $path_to_idat_folder \
      --output-gtc \
      --gender-estimate-call-rate-threshold -0.1

**2.  gtc2vcf**

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
