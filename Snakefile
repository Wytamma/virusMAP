# TODO: mapping to a combined genome https://groups.google.com/forum/#!topic/rna-star/TqOdXiEFYrI
# TODO: Magic-Blast
# TODO: minimap2?

configfile: "config.json"

if config["TESTING"]:
    SRA_IDS = ["SRR1553459"] # ebola
    VIRAL_GENBANK_IDS = ["NC_002549"] # ebola
else:
    SRA_IDS = config["SRA_IDS"].split()
    VIRAL_GENBANK_IDS = config["VIRAL_GENBANK_IDS"].split()

ALIGNERS = config["ALIGNERS"].split()

wildcard_constraints:
    #R="^[0-9]$",
    sra_id="^(SRR)[1-9]*$",
    #id="\d$"

localrules: all, multiqc

rule all:
    """
    Collect the main outputs of the workflow.
    """
    input:
        "results/multiqc.html",
        "results/jobgraph.png",
        'rulegraph.png'
    params:
        mem = '1gb'

rule get_genome_gb:
    """
    Retrieve the sequence in gb format from ncbi.
    """
    params:
        mem = '1gb'
    output:
        "data/viral_genome/{genbank_id}.gb"
    run:
        from Bio import Entrez
        from Bio import SeqIO
        Entrez.email = 'wytamma.wirth@me.com'
        with Entrez.efetch(db="nucleotide", id=wildcards.genbank_id, rettype="gb", retmode="text") as handle:
            record = SeqIO.read(handle, "genbank")
            SeqIO.write(record, output[0], "genbank")

rule convert_gb_to_fasta:
    """
    Converts gb format to fasta.
    """
    input:
        "data/viral_genome/{genbank_id}.gb"
    params:
        mem = '1gb'
    output:
        "data/viral_genome/{genbank_id}.fa"
    run:
        from Bio import SeqIO
        count = SeqIO.convert(input[0], "genbank", output[0], "fasta")

rule get_genome_gff:
    """
    Retrieve the sequence in gff format from ncbi.
    """
    params:
        mem = '1gb'
    output:
        "data/viral_genome/{genbank_id}.gff"
    shell:
        """
        curl -L --output {output} "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id={wildcards.genbank_id}" 
        """


rule get_reads_fastq:
    """
    Retrieve the reads in FASTQ format from the EBI SRA database.
    """
    output:
        "data/raw_reads/{srr_id}_{R}.fastq.gz"
    params:
        last_one=lambda wildcards : wildcards.srr_id[-1:],
        first_six=lambda wildcards : wildcards.srr_id[:6],
        mem = '1gb',
    shell:
        """
        curl -L \
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{params.first_six}/00{params.last_one}/{wildcards.srr_id}/{wildcards.srr_id}_{wildcards.R}.fastq.gz \
        -o {output}
        """


rule trim_galore:
    """
    Run Trim Galore! on a FASTQ file pair.

    30mins for 1gb of FASTQ with 4gb of RAM
    """
    input:
        "data/raw_reads/{srr_id}_1.fastq.gz",
        "data/raw_reads/{srr_id}_2.fastq.gz"
    params:
        mem = '4gb'
    threads: 4
    output:
        trimmed_reads = expand("data/trimmed_reads/{{srr_id}}_{R}_trimmed.fq.gz", R = [1,2]),
        trimming_report = expand("intermediate/trimming/{{srr_id}}_{R}.fastq.gz_trimming_report.txt", R = [1,2]),
    shell:
        """
        # Run trim_galore and save the output to the current directory
        trim_galore --paired {input} --output_dir .

        # Move the files which are used in the workflow
        mv {wildcards.srr_id}_1.fastq.gz_trimming_report.txt {output.trimming_report[0]}
        mv {wildcards.srr_id}_2.fastq.gz_trimming_report.txt {output.trimming_report[1]}
        mv {wildcards.srr_id}_1_val_1.fq.gz {output.trimmed_reads[0]}
        mv {wildcards.srr_id}_2_val_2.fq.gz {output.trimmed_reads[1]}
        """

rule fastqc:
    """
    Run FastQC on a trimmed FASTQ file.
    """
    input:
        "data/trimmed_reads/{id}_trimmed.fq.gz"
    params:
        mem = '4gb'
    output:
        "results/fastqc/{id}_trimmed_fastqc.html",
        "intermediate/fastqc/{id}_trimmed_fastqc.zip"
    shell:
        """
        # Run fastQC and save the output to the current directory
        fastqc {input} -q -o .

        # Move the files which are used in the workflow
        mv {wildcards.id}_trimmed_fastqc.html {output[0]}
        mv {wildcards.id}_trimmed_fastqc.zip {output[1]}
        """

rule genomeGenerate:
    """
    Generate genome for mapping.
    """
    input:
        fasta = "data/viral_genome/{genbank_id}.fa",
        gff = "data/viral_genome/{genbank_id}.gff",
    params:
        mem = '2gb'
    output:
        "data/STAR_genome/{genbank_id}/Genome"
    threads: 1
    run:
        from Bio import SeqIO
        from math import log2
        GenomeLength = len(next(SeqIO.parse(input.fasta, "fasta")))
        genomeSAindexNbases = round(min(14, log2(GenomeLength)/2 - 1))
        if 'exon' in open(input.gff).read():
            sjdbGTFfile = f"--sjdbGTFfile {input.gff}"
            sjdbGTFtagExonParentTranscript = "--sjdbGTFtagExonParentTranscript Parent"
        else:
            sjdbGTFfile = ''
            sjdbGTFtagExonParentTranscript = ''
        command = f"""
            STAR --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir data/STAR_genome/{wildcards.genbank_id}/ \
            --outFileNamePrefix data/STAR_genome/{wildcards.genbank_id}/ \
            --genomeFastaFiles {input.fasta} \
            {sjdbGTFfile} {sjdbGTFtagExonParentTranscript} \
            --genomeSAindexNbases {genomeSAindexNbases}"""
        shell(command)

        
rule STAR:
    """
    Map the trimmed reads to the genome with STAR.
     --outFileNamePrefix intermediate/STAR/
    """
    input:
        genome = "data/STAR_genome/{genbank_id}/Genome",
        trimmed_reads = expand("data/trimmed_reads/{{srr_id}}_{R}_trimmed.fq.gz", R = [1,2])
    params:
        mem = '8gb',
        datadirprefix = 'data/mapped_reads/{srr_id}.{genbank_id}.STAR',
        prefix = '{srr_id}.{genbank_id}.STAR'
    output:
        "data/mapped_reads/{srr_id}.{genbank_id}.STAR.bam",
        "intermediate/STAR/{srr_id}.{genbank_id}.Log.final.out"
    threads: 4
    shell:
        """
        STAR --runThreadN {threads} \
        --genomeDir data/STAR_genome/{wildcards.genbank_id} \
        --readFilesIn {input.trimmed_reads} \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "{params.datadirprefix}." \
        --readFilesCommand gunzip -c \
        --outReadsUnmapped Fastx

        # move log file 
        mv {params.datadirprefix}.Log.final.out {output[1]}

        # move unmapped 
        mkdir -p data/unmapped_reads/
        mv {params.datadirprefix}.Unmapped.out.mate1 data/unmapped_reads/{params.prefix}.Unmapped.mate1.fq
        mv {params.datadirprefix}.Unmapped.out.mate2 data/unmapped_reads/{params.prefix}.Unmapped.mate2.fq

        # rename 
        mv {params.datadirprefix}.Aligned.sortedByCoord.out.bam {params.datadirprefix}.bam

        # delete others 
        rm {params.datadirprefix}.Log.*
        rm {params.datadirprefix}.SJ.out.tab
        """

rule bwa_index:
    input:
        "data/viral_genome/{genbank_id}.fa"
    output:
        expand("data/bwa_index/{{genbank_id}}.fa.{extention}", extention = ['bwt', 'ann', 'amb', 'pac', 'sa'])
    params:
        mem = '8gb',
    shell:
        """
        bwa index {input}
        mv {input}.* data/bwa_index/
        """

rule bwa_map:
    input:
        expand("data/bwa_index/{{genbank_id}}.fa.{extention}", extention = ['bwt', 'ann', 'amb', 'pac', 'sa']),
        trimmed_reads = expand("data/trimmed_reads/{{srr_id}}_{R}_trimmed.fq.gz", R = [1,2])
    output:
        "data/mapped_reads/{srr_id}.{genbank_id}.bwa.bam"
    params:
        db_prefix = "data/bwa_index/{genbank_id}.fa",
        mem = '8gb',
    threads: 4
    shell:
        """
        bwa mem -t 8 {params.db_prefix} {input.trimmed_reads} | samtools sort - -@{threads} -o {output}
        """

rule minimap2:
    input:
        ref = "data/viral_genome/{genbank_id}.fa",
        trimmed_reads = expand("data/trimmed_reads/{{srr_id}}_{R}_trimmed.fq.gz", R = [1,2])
    output:
        "data/mapped_reads/{srr_id}.{genbank_id}.minimap2.bam"
    params:
        mem = '8gb',
    threads: 4
    shell:
        """
        # -a flag returns both mapped and unmapped reads.
        minimap2 -ax sr {input.ref} {input.trimmed_reads} | samtools sort - -@{threads} -o {output}
        """

rule qualimap:
    """
    Quality control of alignment
    """ 
    input:
        bam = "data/mapped_reads/{srr_id}.{genbank_id}.{aligner}.bam",
        gff = "data/viral_genome/{genbank_id}.gff",
    params:
        outdir = "intermediate/qualimap/{srr_id}.{genbank_id}.{aligner}",
        mem = '4gb',
        id = '{srr_id}.{genbank_id}.{aligner}'
    output:
        "intermediate/qualimap/{srr_id}.{genbank_id}.{aligner}/qualimapReport.html",
    shell:
        """
        qualimap bamqc \
        -bam {input.bam} \
        -outformat HTML \
        -outdir {params.outdir}/ \
        -gff {input.gff}

        #copy coverage results file
        mkdir -p results/qualimap/

        if [ -f {params.outdir}/images_qualimapReport/genome_coverage_across_reference.png ]; then
            cp {params.outdir}/images_qualimapReport/genome_coverage_across_reference.png \
            results/qualimap/{params.id}.coverage.png
        fi
        """


rule multiqc:
    """
    Aggregate all results into a MultiQC report.
    """
    input:
        expand("intermediate/fastqc/{srr_id}_{R}_trimmed_fastqc.zip", srr_id = SRA_IDS, R = [1,2]),
        expand("intermediate/trimming/{srr_id}_{R}.fastq.gz_trimming_report.txt", srr_id = SRA_IDS, R = [1,2]),
        expand("intermediate/qualimap/{srr_id}.{genbank_id}.{alinger}/qualimapReport.html", srr_id = SRA_IDS, genbank_id = VIRAL_GENBANK_IDS, alinger = ALIGNERS),
    params:
        mem = '2gb'
    output:
        html = "results/multiqc.html",
        stats = "intermediate/multiqc_general_stats.txt"
    wildcard_constraints:
        R="^[0-9]$",
    shell:
        """
        # Run multiQC and keep the html report
        multiqc -n multiqc.html intermediate/
        mv multiqc.html {output.html}
        mv multiqc_data/multiqc_general_stats.txt {output.stats}

        # Remove the other directory that multiQC creates
        rm -rf multiqc_data
        """

rule generate_rulegraph:
    """
    Generate a rulegraph for the workflow.
    """
    params:
        mem = '1gb'
    output:
        "rulegraph.png"
    shell:
        """
        rm -f {output}
        snakemake --rulegraph | dot -Tpng > {output}
        """

rule generate_jobgraph:
    """
    Generate a rulegraph for the workflow.
    """
    params:
        mem = '1gb'
    output:
        "results/jobgraph.png"
    shell:
        """
        rm -f {output}
        snakemake --dag | dot -Tpng > {output}
        """

