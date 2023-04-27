class General:
    @staticmethod
    def get_input_all(wildcards):
        input_list = []
        #input_list += expand(RESULT_DIR / "{sample}/FirstSearch/MS2Rescore/config.json", sample=list(SAMPLES.index))
        #input_list += expand(RESULT_DIR / "{sample}/taxids/mapped_taxids.tsv", sample=list(SAMPLES.index))
        #input_list += expand(RESULT_DIR / "{sample}/RefFilter/Filtered_ref.mgf", sample=list(SAMPLES.index))
        #input_list += expand(RESULT_DIR / "{sample}/Database/concat_proteomes.fasta", sample=list(SAMPLES.index))
        #input_list += expand(RESULT_DIR / "{sample}/Database/proteomes_concatenated_target_decoy.fasta", sample=list(SAMPLES.index))
        #input_list += expand(RESULT_DIR / "{sample}/FinalSearch/proteomes_Default_PSM_Report.txt", sample=list(SAMPLES.index))
        #input_list += expand(RESULT_DIR / "{sample}/PepGM/config/config.yaml", sample=list(SAMPLES.index))
        input_list += expand(RESULT_DIR / "{sample}/refSeqViral_Default_PSM_Report.txt", sample=list(SAMPLES.index))
        #input_list += expand(RESULT_DIR / "{sample}/FetchData/strain_accessions.tsv", sample=list(SAMPLES.index))
        #input_list += expand(RESULT_DIR / "{sample}/taxids/strain_name_counts.tsv", sample=list(SAMPLES.index))

        # Plots
        input_list += expand(RESULT_DIR / "{sample}/Plots/{sample}_strain_bar_plot.png", sample=list(SAMPLES.index))
        input_list += expand(RESULT_DIR / "{sample}/Plots/{sample}_taxIdScores_bar_plot.png", sample=list(SAMPLES.index))
        input_list += expand(RESULT_DIR / "{sample}/Plots/{sample}_proportions_pie_chart.png", sample=list(SAMPLES.index))
        input_list += expand(RESULT_DIR / "{sample}/Plots/{sample}_first_search_confidence_histogram.png", sample=list(SAMPLES.index))
        input_list += expand(RESULT_DIR / "{sample}/Plots/{sample}_final_search_confidence_histogram.png", sample=list(SAMPLES.index))
        
        return input_list

class SearchDB:
    @staticmethod
    def get_input_MGF():
        if not SKIP_HOST_FILTERING:
            data_dict = {
                "mgf": RESULT_DIR / "{sample}/SpectraFilter/Filtered_host.mgf"
                }
        else:
            data_dict = {
                "mgf": MGF_FILE,
                }
        return data_dict
    
    @staticmethod
    def get_output_AddDecoysRef():
        ref_path = config["db_search"]["ref"]
        ref_path = ref_path.split(".")
        ref_name = ""
        num_name_segments = len(ref_path[:-1])
        for i in range(num_name_segments):
            if i == num_name_segments - 1:
                ref_name += ref_path[i]
            else:
                ref_name += ref_path[i] + "."

        ref_name += "_concatenated_target_decoy." + ref_path[-1]

        data_dict = {
            "ref_decoy_fasta": ref_name 
        }
        return data_dict

class Mapping:
    @staticmethod
    def get_input_getTargets():
        ref_path = str(config["db_search"]["ref"])
        input_list = [
            RESULT_DIR / "{sample}/{ref_path}_Default_PSM_Report.txt",
            "resources/taxidMapping/accessions_hashed.npy",
            "resources/taxidMapping/taxids.txt"
        ]
        return input_list

class Filtering:
    @staticmethod
    def get_input_filter_undefined():
        if not SKIP_HOST_FILTERING:
            data_dict = {
                "mgf": RESULT_DIR / "{sample}/SpectraFilter/Filtered_host.mgf"
                }
        else:
            data_dict = {
                "mgf": MGF_FILE,
                }
        return data_dict
 
        