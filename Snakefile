###########################
# module to be loaded
#module load snakemake-5.28.0 
#module load STAR-2.7.6a 
#module load fastqc-0.11.9 
#module load picard-2.23.8
#module load subread-2.0.1
###########################


configfile:
    "config.json"

workdir:
   config['data']

SAMPLES, = glob_wildcards(config['data'] + "/{sample}_1.fq.gz")
RESULTS = config['data'] + '/results'
LOGS = RESULTS + '/logs'


rule all:
    input:
        #idx = directory('index'),
        bam1 = expand("{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        bam2 = expand("{sample}.Aligned.sortedByCoord.markDup.out.bam", sample=SAMPLES),
        fqc1 = expand(LOGS + '/{sample}.fastqc_before_md.log', sample=SAMPLES),
        fqc2 = expand(LOGS + '/{sample}.fastqc_after_md.log', sample=SAMPLES),
        fc1 = RESULTS + '/featureCounts_genes.txt',
        fc2 = RESULTS + '/featureCounts_transcripts.txt'


rule index:
    input:
        fa = config['ref'], # provide your reference FASTA file
        gtf = config['gtf'] # provide your GTF file
    output:
        directory('index') # you can rename the index folder
    threads: 
        24 # set the maximum number of available cores
    shell:
        'mkdir {output} && '
        'STAR --runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {output} '
        '--genomeFastaFiles {input.fa} ' #'--genomeFastaFiles <(zcat {input.fa}) '
        '--sjdbGTFfile {input.gtf} '
        '--sjdbOverhang 100'


rule fastqc_before:
    input:
        R1 = config['data'] + "/{sample}_1.fq.gz",
        R2 = config['data'] + "/{sample}_2.fq.gz",
    output:
        out = directory(RESULTS + "/fastqc/before_md/{sample}")
    log:
        LOGS + '/{sample}.fastqc_before_md.log'
    threads:
        48 # set the maximum number of available cores
    shell:
        """
        mkdir {output.out}
        fastqc {input.R1} {input.R2} -o {output.out} >> {log} 2>&1
        """


rule align_sort:
    input:
        R1 = config['data'] + "/{sample}_1.fq.gz",
        R2 = config['data'] + "/{sample}_2.fq.gz",
        idx = rules.index.output
        #idx = directory('index')
    output:
        #temp("{sample}.mapped.bam")
        bam = '{sample}.Aligned.sortedByCoord.out.bam',
        sj = '{sample}.SJ.out.tab',
        log = '{sample}.Log.final.out'
    params:
        out = '{sample}.'
    log:
        LOGS + '/{sample}.STAR.log'
    threads:
        24 # set the maximum number of available cores
    shell:
        'STAR --runThreadN {threads} '
            '--genomeDir {input.idx} '
            '--readFilesIn <(zcat {input.R1}) <(zcat {input.R2}) '
            '--outSAMtype BAM SortedByCoordinate ' 
            #'--readFilesCommand zcat ' 
            '--quantMode GeneCounts '
            '--outFileNamePrefix {params.out} >> {log} 2>&1'
        #'mv results/Aligned.sortedByCoord.out.bam {output.bam}'
        #'mv results/SJ.out.tab {output.sj}'
        #'mv results/Log.out {output.log}'


rule mark_duplicates:
    input:
        bam = rules.align_sort.output.bam,
    output:
        bam = '{sample}.Aligned.sortedByCoord.markDup.out.bam',
        met = '{sample}.marked_dup_metrics.txt'
    log:
        LOGS + '/{sample}.markDup.log'
    threads:
        48 # set the maximum number of available cores
    shell:
        """
        java -jar $PICARDHOME/picard.jar MarkDuplicates \
        --REMOVE_DUPLICATES \
        -I {input.bam} \
        -O {output.bam} \
        -M {output.met} >> {log} 2>&1
        """


rule fastqc_after:
    input:
        bam = rules.mark_duplicates.output.bam
    output:
        out = directory(RESULTS + "/fastqc/after_md/{sample}")
    log:
        LOGS + '/{sample}.fastqc_after_md.log'
    threads:
        48 # set the maximum number of available cores
    shell:
        """
        mkdir {output.out}
        fastqc {input.bam} -o {output.out} >> {log} 2>&1
        """


rule featureCounts:
    input:
        bam = expand("{sample}.Aligned.sortedByCoord.markDup.out.bam", sample=SAMPLES),
        gtf = config['gtf'] # provide your GTF file
    output:
        res1 = RESULTS + '/featureCounts_genes.txt',
        res2 = RESULTS + '/featureCounts_transcripts.txt'
    log:
        log1 = LOGS + '/featureCounts_genes.log',
        log2 = LOGS + '/featureCounts_transcripts.log'
    threads: 
        24
    shell:
        'featureCounts -a {input.gtf} '
        '-o {output.res1} '
        '-T {threads} {input.bam} >> {log.log1} 2>&1 \n'
        
        'featureCounts -a {input.gtf} '
        '-g transcript_id '
        '-o {output.res2} '
        '-T {threads} {input.bam} >> {log.log2} 2>&1'

