rule SearchSpectraAgainstProteomes:
    input: 
        proteomes_decoy_fasta = RESULT_DIR / "{sample}/Database/proteomes_concatenated_target_decoy.fasta",
        mgf = RESULT_DIR / "{sample}/RefFilter/Filtered_ref.mgf",
        par = PAR_FILE,
    output:  
        out_zip = RESULT_DIR / "{sample}/FinalSearch/proteomes_searchgui_out.zip"
    log:
        stderr_log=RESULT_DIR / "logs/FinalSearch/SearchSpectraAgainstProteomes/{sample}/stderr.log",
        stdout_log=RESULT_DIR / "logs/FinalSearch/SearchSpectraAgainstProteomes/{sample}/stdout.log"
    params:
        result_dir = str(RESULT_DIR / "{sample}/FinalSearch"),
        refname = "proteomes",
        search_engine = config["db_search"]["search_engine"],
        psm_fdr = config["db_search"]["psm_fdr"],
        peptide_fdr = config["db_search"]["peptide_fdr"],
        protein_fdr = config["db_search"]["protein_fdr"]
    # TODO: Conda variante zum Laufen bekommen
    conda:
        "../envs/java.yml"
    threads: workflow.cores / 2
    shell: 
        "java -cp /home/jpipart/project/SearchGUI-4.2.7/SearchGUI-4.2.7.jar eu.isas.searchgui.cmd.SearchCLI -spectrum_files {input.mgf} -fasta_file {input.proteomes_decoy_fasta} -output_folder {params.result_dir} -id_params {input.par} -output_default_name {params.refname}_searchgui_out -psm_fdr 1 -peptide_fdr 1 -protein_fdr 1 {params.search_engine} 1 -threads {threads} > {log.stdout_log} 2> {log.stderr_log}"


rule RunPeptideShakerProteomes:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/FinalSearch/proteomes_searchgui_out.zip",
        mgf = RESULT_DIR / "{sample}/RefFilter/Filtered_ref.mgf",
        proteomes_decoy_fasta = RESULT_DIR / "{sample}/Database/proteomes_concatenated_target_decoy.fasta"
    output: 
        peptide_shaker_psdb = RESULT_DIR / "{sample}/FinalSearch/proteomes.psdb"
    log:
        stderr_log=RESULT_DIR / "logs/FinalSearch/RunPeptideShakerProteomes/{sample}/stderr.log",
        stdout_log=RESULT_DIR / "logs/FinalSearch/RunPeptideShakerProteomes/{sample}/stdout.log",
        peptide_shaker_log=RESULT_DIR / "logs/FinalSearch/RunPeptideShakerProteomes/{sample}/PeptideShaker.log",
    params:
        refname = "proteomes",
    # TODO: Conda variante zum Laufen bekommen
    conda:
        "../envs/java.yml"
    threads: workflow.cores / 2
    shell:
        "java -cp /home/jpipart/project/PeptideShaker-2.2.22/PeptideShaker-2.2.22.jar eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference {params.refname} -fasta_file {input.proteomes_decoy_fasta} -identification_files {input.searchgui_zip} -spectrum_files {input.mgf} -out {output.peptide_shaker_psdb} -threads {threads} -log {log.peptide_shaker_log} > {log.stdout_log} 2> {log.stderr_log}" 


rule SimplePeptideListProteomes:
    input:
        peptide_shaker_psdb = RESULT_DIR / "{sample}/FinalSearch/proteomes.psdb"
    output: 
        #TODO hostname varialbe: ResultsDir+SampleName+'/FinalSearch/'+HostName+'.psdb'
        peptide_shaker_report = RESULT_DIR / "{sample}/FinalSearch/proteomes_Default_PSM_Report.txt",
    log:
        stderr_log=RESULT_DIR / "logs/FinalSearch/SimplePeptideListProteomes/{sample}/stderr.log",
        stdout_log=RESULT_DIR / "logs/FinalSearch/SimplePeptideListProteomes/{sample}/stdout.log"
    params:
        out_dir = str(RESULT_DIR / "{sample}/FinalSearch")
    # TODO: Conda variante zum Laufen bekommen
    conda:
        "../envs/java.yml"
    threads: 1
    shell:
        "java -cp /home/jpipart/project/PeptideShaker-2.2.22/PeptideShaker-2.2.22.jar eu.isas.peptideshaker.cmd.ReportCLI -in {input.peptide_shaker_psdb} -out_reports {params.out_dir} -reports 3 > {log.stdout_log} 2> {log.stderr_log}" 
   
# rule extractSearchGuiResults:
#     input:
#         searchgui_zip = RESULT_DIR / "{sample}/FirstSearch/ref_searchgui_out.zip",
#     output:
#         out_file = RESULT_DIR / "{sample}/FirstSearch/Filtered_host.t.xml.gz"
#     log:
#         stderr_log=RESULT_DIR / "logs/FirstSearch/extractSearchGuiResults/{sample}/stderr.log",
#         stdout_log=RESULT_DIR / "logs/FirstSearch/extractSearchGuiResults/{sample}/stdout.log"  
#     params:
#         out_dir = str(RESULT_DIR / "{sample}/FirstSearch")
#     threads: 1
#     shell:
#         "unzip -u {input.searchgui_zip} -d {params.out_dir} > {log.stdout_log}"


# rule createMS2RescoreConfig:
#     input: 
#         xTandem_out = RESULT_DIR / "{sample}/FirstSearch/Filtered_host.t.xml.gz"
#     output:
#         out_dir = RESULT_DIR / "{sample}/FirstSearch/MS2Rescore/config.json"
#     params:
#         pipeline = config["ms2rescore"]["RescorePipeline"],
#         features = config["ms2rescore"]["RescoreFeatures"],
#         percolator = config["ms2rescore"]["RunPercolator"],
#         fragmodel = config["ms2rescore"]["FragModel"]
#     threads: 1
#     script:
#         "../scripts/ms2rescore_config.py"

# rule RunMS2Rescore:
#     input:
#         mgf = SearchDB.get_input_MGF()["mgf"],
#         tandem_xml = RESULT_DIR / "{sample}/FirstSearch/Filtered_host.t.xml.gz"
#     output: 
#         #TODO hostname varialbe: ResultsDir+SampleName+'/FirstSearch/'+HostName+'.psdb'
#         ms2rescore_out = RESULT_DIR / "{sample}/FirstSearch/ref_Default_PSM_Report.txt",
#     log:
#         stderr_log=RESULT_DIR / "logs/FirstSearch/RunMS2Rescore/{sample}/stderr.log",
#         stdout_log=RESULT_DIR / "logs/FirstSearch/RunMS2Rescore/{sample}/stdout.log"
#     params:
#         out_dir = str(RESULT_DIR / "{sample}/FirstSearch/MS2Rescore")
#     conda:
#         "../envs/ms2rescore.yml"
#     threads: 1
#     shell:
#         "peptide-shaker eu.isas.peptideshaker.cmd.ReportCLI -in {input.peptide_shaker_psdb} -out_reports {params.out_dir} -reports 3 > {log.stdout_log} 2> {log.stderr_log}"     




# rule getTargets:
#     input:
#         PSM_Report = RESULT_DIR / "{sample}/" + config["db_search"]["ref"] + "_Default_PSM_Report.txt",
#         ResourcesDir + TaxidMapping + 'accessions_hashed.npy',
#         ResourcesDir + TaxidMapping + 'taxids.txt'
#     params:
#         query=ResultsDir + SampleName + '/{DBname}_query_accessions.txt',
#         samplename=SampleName,
#         hostname=HostName,
#         DBname=ReferenceDBName
#     conda: 'envs/graphenv.yml'
#     output:
#         ResultsDir + SampleName + '/{DBname}_mapped_taxids.txt',ResultsDir + SampleName + '/{DBname}_mapped_taxids_weights.csv'
#     shell: "python3 workflow/scripts/getTargets.py -rq {input[0]} -q {params.query} -d {input[1]} -t {input[2]} -r {output[0]} "
