configfile: "config.yaml"

with open(config['TREAT1']) as fp:
    TREAT1 = fp.read().splitlines()
with open(config['TREAT2'] ) as fp: 
    TREAT2 = fp.read().splitlines () 

with open(config['CONTROL']) as fp:
    CONTROL = fp.read().splitlines()

print(TREAT1)
print(TREAT2)
print(CONTROL)

rule all: 
     input: 
         expand("{group}_SALMON.{sample}/quant.sf", sample = TREAT1, group = config['TREAT1_NAME']),
         expand("{group}_SALMON.{sample}/quant.sf", sample = TREAT2, group = config['TREAT2_NAME']), 
         expand("{group}_SALMON.{sample}/quant.sf", sample = CONTROL, group = config['CONTROL_NAME']),
         "dexseq_isoformsfeatures.csv", 

rule trim:
       input:
           r1 = "{sample}.r_1.fq.gz",
           r2 = "{sample}.r_2.fq.gz"
       output:
           "galore/{sample}.r_1_val_1.fq.gz",
           "galore/{sample}.r_2_val_2.fq.gz"
       conda: 'env/env-trim.yaml'
       shell:
           """ 
           mkdir -p galore
           mkdir -p fastqc
           trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2}
           """

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
         
rule isoform_all:
      input: 
          expand("{group}_SALMON.{sample}/quant.sf", sample = TREAT1, group = config['TREAT1_NAME']), 
          expand("{group}_SALMON.{sample}/quant.sf", sample = TREAT2, group = config['TREAT2_NAME']), 
          expand("{group}_SALMON.{sample}/quant.sf", sample = CONTROL, group = config['CONTROL_NAME']) 
      params: 
           config['CONTROL_NAME'], 
           config['TREAT1_NAME'], 
           config['TREAT2_NAME'],
           config['Nc'], 
           config['Nt1'], 
           config['Nt2']  
      output: 
         "dexseq_isoformsfeatures.csv"
      shell: 
          """ 
          mkdir -p {params[1]}2 
          mkdir -p {params[2]}2
          cp -r {params[1]}_SALMON.* {params[1]}2 
          cp -r {params[0]}_SALMON.* {params[1]}2
         
          cp -r {params[2]}_SALMON.* {params[2]}2
          cp -r {params[0]}_SALMON.* {params[2]}2
          
          Rscript {config[SRC]}/splice_all.R {params[0]} {params[1]} {params[2]} {params[3]} {params[4]} {params[5]}
          cd {params[1]}2
          Rscript {config[SRC]}/splice.R {params[0]} {params[1]} {params[3]} {params[4]} 
          cd ..
          cd {params[2]}2 
          Rscript {config[SRC]}/splice.R {params[0]} {params[2]} {params[3]} {params[5]} 
          """  

