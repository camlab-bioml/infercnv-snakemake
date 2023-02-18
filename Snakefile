
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
        leiden_res = config['leiden_res']
    threads: 8
    # log: out = "logs/{sample}_stdout.log",
    #      err = "logs/{sample}_stderr.err"
    resources:
        mem_mb = 150000
    container:
        "docker://trinityctat/infercnv"
    output:
        png = output + 'infercnv/{sample}/infercnv.png',
        obj = output + 'infercnv/{sample}/20_HMM_pred.repr_intensitiesHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.infercnv_obj'
    shell:
        'Rscript scripts/run_infercnv.R --sce {input.sce} --out_dir {params.out_dir} --genelist {input.genelist} --annotation_column {params.annotation_column} --tumour_type "{params.tumour_type}" --normal_type "{params.normal_type}" --cutoff {params.cutoff} --denoise {params.denoise} --threads {threads} --leiden_res {params.leiden_res} --sample {wildcards.sample}'