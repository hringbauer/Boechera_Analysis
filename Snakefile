rule transform_cm:
    input:
        "Data/LinkageMap/7LGs_0916.fmap",
	"Data/IBD_raw/{sample}.ibd"
    output:
        "Data/IBD_cm/{sample}.ibd"
    shell:
        "python3 BlockMergeScript/transform_cm.py {input}"


rule merge_ibd:
    input:
        "Data/IBD_cm/{sample}.ibd"
    output:
        "Data/IBD_merged/{sample}.ibd"
    shell:
        "python3 BlockMergeScript/cm_remove_gap.py {input}"
