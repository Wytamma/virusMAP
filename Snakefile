# TODO: config file for samples

SRA_IDS = ["SRR7063616", "SRR7063617", "SRR7063617", "SRR7063618", "SRR7063619", "SRR7063620", "SRR7063621", "SRR7063622", "SRR7063623", "SRR7063624", "SRR7063625"] # SRR1553459 -ebola for testing
VIRAL_GENBANK_IDS = ["EU493091", "NC_011530"] # NC_002549 -ebola for testing

#SRA_IDS = ["SRR1553459"]
#VIRAL_GENBANK_IDS = ["NC_002549"]

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
        "results/jobgraph.pdf",
        'rulegraph.pdf'
    params:
        pmem = '1GB'

rule get_genome_gb:
    """
    Retrieve the sequence in gb format from ncbi.
    """
    params:
        pmem = '1GB'
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
        pmem = '1GB'
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
        pmem = '1GB'
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
        pmem = '1GB',
    shell:
        """
        curl -L \
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{params.first_six}/00{params.last_one}/{wildcards.srr_id}/{wildcards.srr_id}_{wildcards.R}.fastq.gz \
        -o {output}
        """


rule trim_galore:
    """
    Run Trim Galore! on a FASTQ file pair.
    """
    input:
        "data/raw_reads/{srr_id}_1.fastq.gz",
        "data/raw_reads/{srr_id}_2.fastq.gz"
    params:
        pmem = '4GB'
    output:
        trimmed_reads = expand("data/trimmed_reads/{{srr_id}}_{R}_trimmed.fq.gz", R = [1,2]),
        trimming_report = expand("intermediate/trimming/{{srr_id}}_{R}.fastq.gz_trimming_report.txt", R = [1,2]),
    threads: 4
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
        pmem = '2GB'
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
        pmem = '2GB'
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
        pmem = '8GB'
    output:
        "data/mapped_reads/{srr_id}.{genbank_id}.Aligned.sortedByCoord.out.bam",
        "intermediate/STAR/{srr_id}.{genbank_id}.Log.final.out"
    threads: 4
    shell:
        """
        STAR --runThreadN {threads} \
        --genomeDir data/STAR_genome/{wildcards.genbank_id} \
        --readFilesIn {input.trimmed_reads} \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "data/mapped_reads/{wildcards.srr_id}.{wildcards.genbank_id}." \
        --readFilesCommand gunzip -c 

        # move log file 
        mv data/mapped_reads/{wildcards.srr_id}.{wildcards.genbank_id}.Log.final.out {output[1]}

        # delete others 
        rm data/mapped_reads/{wildcards.srr_id}.{wildcards.genbank_id}.Log.*
        rm data/mapped_reads/{wildcards.srr_id}.{wildcards.genbank_id}.SJ.out.tab
        """

rule qualimap:
    """
    Quality control of alignment
    """ 
    input:
        bam = "data/mapped_reads/{srr_id}.{genbank_id}.Aligned.sortedByCoord.out.bam",
        gff = "data/viral_genome/{genbank_id}.gff",
    params:
        outdir = "intermediate/qualimap/{srr_id}.{genbank_id}",
        pmem = '4GB'
    output:
        "intermediate/qualimap/{srr_id}.{genbank_id}/qualimapReport.html",
    shell:
        """
        qualimap bamqc \
        -bam {input.bam} \
        -outformat HTML \
        -outdir {params.outdir}/ \
        -gff {input.gff}

        #copy coverage results file
        mkdir -p results/qualimap/
        cp {params.outdir}/images_qualimapReport/genome_coverage_across_reference.png results/qualimap/{wildcards.srr_id}.{wildcards.genbank_id}.coverage.png
        """


rule multiqc:
    """
    Aggregate all results into a MultiQC report.
    """
    input:
        expand("intermediate/STAR/{srr_id}.{genbank_id}.Log.final.out", srr_id = SRA_IDS, genbank_id = VIRAL_GENBANK_IDS),
        expand("intermediate/fastqc/{srr_id}_{R}_trimmed_fastqc.zip", srr_id = SRA_IDS, R = [1,2]),
        expand("intermediate/trimming/{srr_id}_{R}.fastq.gz_trimming_report.txt", srr_id = SRA_IDS, R = [1,2]),
        expand("intermediate/qualimap/{srr_id}.{genbank_id}/qualimapReport.html", srr_id = SRA_IDS, genbank_id = VIRAL_GENBANK_IDS),
    params:
        pmem = '2GB'
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
        pmem = '1GB'
    output:
        "rulegraph.pdf"
    shell:
        """
        rm -f {output}
        snakemake --rulegraph | dot -Tpdf > {output}
        """

rule generate_jobgraph:
    """
    Generate a rulegraph for the workflow.
    """
    params:
        pmem = '1GB'
    output:
        "results/jobgraph.pdf"
    shell:
        """
        rm -f {output}
        snakemake --dag | dot -Tpdf > {output}
        """

