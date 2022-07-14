configfile: "config.yaml"

with open(config['TREAT']) as fp:
    TREAT = fp.read().splitlines()

with open(config['CONTROL']) as fp:
    CONTROL = fp.read().splitlines()

print(TREAT)
print(CONTROL)

rule all: 
     input: 
         expand("{group}_SALMON.{sample}/quant.sf", sample = TREAT, group = config['TREAT_NAME']), 
         expand("{group}_SALMON.{sample}/quant.sf", sample = CONTROL, group = config['CONTROL_NAME']),
         "dexseq_isoformsfeatures.csv"
rule quant:  
      input: 
            r1 = "galore/{sample}.r_1_val_1.fq.gz",
            r2 = "galore/{sample}.r_2_val_2.fq.gz",
      params: 
            index = config['SALMON_INDEX'],
            lib = config['SALMON_LIBRARY'],
            outdir = "{group}_SALMON.{sample}"
      conda: 'env/env-quant.yaml'
      output:
                  "{group}_SALMON.{sample}/quant.sf",
      shell:
                 """
                 salmon quant -i {params.index} --libType {params.lib} -1 {input.r1} -2 {input.r2} -o {params.outdir} -p 3 --writeUnmappedNames --seqBias --gcBias --validateMappings
                 """ 
         
rule isoform: 
      output: 
         "dexseq_isoformsfeatures.csv"
      shell: 
           "Rscript splice.R"
