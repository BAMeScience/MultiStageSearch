rule genome2Proteome:
    input: 
        concat_strain_genomes = RESULT_DIR / "{sample}/FetchData/concat_strain_genomes.fasta"
    output:
        touch(RESULT_DIR / "{sample}/sixpack/genome2Proteome.done"),
        sixpack_temp = RESULT_DIR / "{sample}/sixpack/sixpack_temp.txt"
    params:
        sample_path = str(RESULT_DIR / "{sample}"),
        orfminsize = config["sixpack"]["orfminsize"],
        additional_params = config["sixpack"]["additional_parameters"],
    conda:
        "../envs/sixpack.yml"
    threads: 1
    script:
        "../scripts/genome2proteome.py"

rule concatProteomes:
    input:
        sixpack_done = RESULT_DIR / "{sample}/sixpack/genome2Proteome.done"
    output:
        concat_proteomes = RESULT_DIR / "{sample}/Database/proteomes.fasta"
    params:
        res_dir = RESULT_DIR,
        sample_name = "{sample}",
    conda:
        "../envs/base_python.yml"
    threads: 1
    script:
        "../scripts/concat_proteomes.py"

rule AddDecoysProteome:
    input: 
        proteomes = RESULT_DIR / "{sample}/Database/proteomes.fasta"
    output: 
        proteomes_decoy = RESULT_DIR / "{sample}/Database/proteomes_concatenated_target_decoy.fasta"
    log:
        stderr_log = RESULT_DIR / "logs/genome2proteome/proteomeDB/{sample}/stderr.log",
        stdout_log = RESULT_DIR / "logs/genome2proteome/proteomeDB/{sample}/stdout.log" 
    params:
        searchgui = config["SearchGUI"],
    conda:
        "../envs/java.yml"
    threads: 1
    shell: 
        "java -cp {params.searchgui} eu.isas.searchgui.cmd.FastaCLI -in {input.proteomes} -decoy > {log.stdout_log} 2> {log.stderr_log}"

rule FilterDuplicateORFs:
    input:
        proteomes_decoy = RESULT_DIR / "{sample}/Database/proteomes_concatenated_target_decoy.fasta"
    output:
        proteomes_decoy_unique = RESULT_DIR / "{sample}/Database/proteomes_concatenated_target_decoy_unique.fasta"
    log:
        stderr_log = RESULT_DIR / "logs/genome2proteome/FilterDuplicates/{sample}/stderr.log",
        stdout_log = RESULT_DIR / "logs/genome2proteome/FilterDuplicates/{sample}/stdout.log" 
    conda:
        "../envs/fetch_data.yml"
    threads: 1
    script:
        "../scripts/filter_duplicate_ORFs.py"