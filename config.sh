
export PATH="$HOME/bin:$PATH"
export BCFTOOLS_PLUGINS="$HOME/bin"

##########################################
## work folder
##########################################
export work_path=/public/home/naq/project/20210127_LC_mosaic

##########################################
## script
##########################################
export scripts=${work_path}/scripts

##########################################
## config
##########################################
export config=${work_path}/config

##########################################
## thread
##########################################
export phase_threads=2
export impute_threads=8

##########################################
## input folder
##########################################
export output_path=${work_path}/result
export data_path=${work_path}/data
export GSA_manifest_path=${data_path}/GSA_manifest
export vcf_path=${output_path}/vcf
export phaseref_path=${output_path}/phaseref
export imputeref_path=${output_path}/imputeref
export phase_path=${output_path}/phase
export impute_path=${output_path}/impute
export mocha_path=${output_path}/mocha
export path_to_gtc_folder=${output_path}/gtc
export log_path=${work_path}/log
export config_path=${work_path}/config
export scripts_path=${work_path}/scripts


##########################################
## reference
##########################################
export ref="/public/home/naq/ref/GRCh3[78]/GRCh37/human_g1k_v37.fasta"

export map="/public/home/naq/ref/GRCh3[78]/GRCh37/genetic_map_hg19_withX.txt.gz"
export panel_pfx="/public/home/naq/ref/GRCh3[78]/GRCh37/1000Genomes/ALL.chr"
export panel_sfx=".phase3_integrated.20130502.genotypes"

export mhc_reg="6:27486711-33448264"
export kir_reg="19:54574747-55504099"

export rule="GRCh37"

export cnp="/public/home/naq/ref/GRCh3[78]/GRCh37/cnps.bed"
export dup="/public/home/naq/ref/GRCh3[78]/GRCh37/segdups.bed.gz"
export cyto="/public/home/naq/ref/GRCh3[78]/GRCh37/cytoBand.txt.gz"

export bpm_manifest_file=${GSA_manifest_path}/GSA-24v1-0_C1.bpm
export csv_manifest_file=${GSA_manifest_path}/GSA-24v1-0_C1.csv
export egt_cluster_file=${GSA_manifest_path}/GSA-24v1-0_C1_ClusterFile.egt


if [ ! -d ${output_path} ]
then
mkdir ${output_path}
fi

if [ ! -d ${data_path} ]
then
mkdir ${data_path}
fi

if [ ! -d ${GSA_manifest_path} ]
then
mkdir ${GSA_manifest_path}
fi

if [ ! -d ${vcf_path} ]
then
mkdir ${vcf_path}
fi

if [ ! -d ${phaseref_path} ]
then
mkdir ${phaseref_path}
fi

if [ ! -d ${imputeref_path} ]
then
mkdir ${imputeref_path}
fi

if [ ! -d ${phase_path} ]
then
mkdir ${phase_path}
fi

if [ ! -d ${impute_path} ]
then
mkdir ${impute_path}
fi

if [ ! -d ${mocha_path} ]
then
mkdir ${mocha_path}
fi

if [ ! -d ${path_to_gtc_folder} ]
then
mkdir ${path_to_gtc_folder}
fi

if [ ! -d ${log_path} ]
then
mkdir ${log_path}
fi

if [ ! -d ${config_path} ]
then
mkdir ${config_path}
fi

if [ ! -d ${scripts_path} ]
then
mkdir ${scripts_path}
fi
