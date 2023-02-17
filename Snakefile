
configfile: 'config/config.yml'

samples = config['samples']

output = 'output/' #+ config['version'] + '/'

rule all:
    input:
        expand(output + 'infercnv/{sample}/infercnv.png', sample = samples)


rule run_infer_cnv:
    input:
        sce = config['input_rds'],
        genelist = "resources/hg38_gencode_v27.txt"
    params:
        out_dir = output + 'infercnv/{sample}/',
        annotation_column = config['annotation_column'],
        tumour_type = config['tumour_type'],
        normal_type = config['normal_type'],
        cutoff = config['cutoff'],
        denoise = config['denoise'],
    threads: 8
    resources:
        mem_mb = 150000
    container:
        "docker://trinityctat/infercnv"
    output:
        png = output + 'infercnv/{sample}/infercnv.png',
        obj = output + 'infercnv/{sample}/20_HMM_pred.repr_intensitiesHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.infercnv_obj'
    script:
        'scripts/run_infercnv.R'