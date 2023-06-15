import os
import glob
import pysam
import snakemake
import subprocess as sp




out_path = os.getcwd()
conditions = [os.path.basename(i) for i in glob.glob(os.path.join(out_path, "*"))]
samples = [os.path.basename(j) for i in conditions for j in glob.glob(os.path.join(out_path, i, "*"))]
data_dict = dict(zip(conditions, samples))

config = {"data":data_dict,
          "indir": out_path,
          "outdir": out_path,
	      "ref": out_path}


rule all:
	input:
		[expand(os.path.join("{condition}", "{sample}", "tumor.merged.realign.bam.bai"),
                condition=key, sample=value) for key, value in data_dict.items()]


rule add_RX: #not tested
    input:
        os.path.join("{condition}", "{sample}", "tumor.merged.bam")
    output:
        os.path.join("{condition}", "{sample}", "tumor.merged.RX.bam")
    run:
        bam = pysam.AlignmentFile(input, 'rb')
        out_bam = pysam.AlignmentFile(output, "wb", template=bam)
        iter = bam.fetch(until_eof=True, multiple_iterators=True)
        for rec in iter:
            fixed_umi = rec.qname.split(":")[-1]
            fixed_umi = fixed_umi.split("_")[1:]
            fixed_umi = '-'.join(fixed_umi)
            rec.set_tag("RX", fixed_umi)
            out_bam.write(rec)
        bam.close()
        out_bam.close()



rule sort_tag:
    input:
        rx_bam = os.path.join("{condition}", "{sample}", "tumor.merged.RX.bam")
    output:
        tagged = os.path.join("{condition}", "{sample}", "tumor.merged.fgbiotagged.bam")
    shell:
        """
        fgbio SortBam -s Queryname -i {input} -o /dev/stdout |\
        fgbio SetMateInformation -i /dev/stdin -o {output.tagged}
        """

rule umi_group:
    input:
        os.path.join("{condition}", "{sample}", "tumor.merged.fgbiotagged.bam")
    output:
        hist = os.path.join("{condition}", "{sample}", "size_hist"),
        grouped = os.path.join("{condition}", "{sample}", "tumor.merged.fgbiogrouped.bam")
    shell:
        """
        fgbio GroupReadsByUmi -i {input} -s adjacency \
        --family-size-histogram {output.hist} -o {output.grouped}
        """

rule make_consensus:
    input:
        os.path.join("{condition}", "{sample}", "tumor.merged.fgbiogrouped.bam")
    output:
        bam = os.path.join("{condition}", "{sample}", "tumor.merged.fgbioconsensus.bam")
    shell:
        """
        fgbio SortBam -s TemplateCoordinate -i {input} -o /dev/stdout | \
        fgbio CallMolecularConsensusReads --min-reads 3 --threads 6 \
        -i /dev/stdin --min-input-base-quality 20 -o {output.bam}
        """


rule filter_consensus:
    input:
        os.path.join("{condition}", "{sample}", "tumor.merged.fgbioconsensus.bam")
    output:
        R1_consensus = os.path.join("{condition}", "{sample}_R1_consensus.fastq.gz"),
        R2_consensus = os.path.join("{condition}", "{sample}_R2_consensus.fastq.gz")
    params:
        ref = "human_g1k_v37.fasta",
        tmp_dir ="/tmp/"
    threads:
        8
    shell:
        """
        fgbio FilterConsensusReads -i {input} -o /dev/stdout -r {params.ref} -M 4 -N 1 -E 0.010 \
        | samtools sort -@ 6 -n \
        | samtools fastq -1 {output.R1_consensus} -2 {output.R2_consensus} -
        """

rule align_consensus:
    input:
        R1_consensus = os.path.join("{condition}", "{sample}_R1_consensus.fastq.gz"),
        R2_consensus = os.path.join("{condition}", "{sample}_R2_consensus.fastq.gz")
    output:
        bam = os.path.join("{condition}", "{sample}", "tumor.merged.realign.bam"),
        index = os.path.join("{condition}", "{sample}", "tumor.merged.realign.bam.bai")
    threads:
        8
    params:
        ref = "human_g1k_v37.fasta"
    shell:
        """
        bwa mem -t 6 {params.ref} -M {input.R1_consensus} {input.R2_consensus}\
        | samtools sort -@6 -o {output.bam} -
        samtools index {output.bam}
        """
