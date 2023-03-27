class General:
    @staticmethod
    def get_input_all(wildcards):
        input_list = []
        input_list += expand(RESULT_DIR / "{sample}/FirstSearch/MS2Rescore/config.json", sample=list(SAMPLES.index))
        input_list += expand(RESULT_DIR / "{sample}/taxids/mapped_taxids.tsv", sample=list(SAMPLES.index))
        input_list += expand(RESULT_DIR / "{sample}/RefFilter/Filtered_ref.mgf", sample=list(SAMPLES.index))
        #input_list += expand(RESULT_DIR / "{sample}/Database/concat_proteomes.fasta", sample=list(SAMPLES.index))
        input_list += expand(RESULT_DIR / "{sample}/Database/proteomes_concatenated_target_decoy.fasta", sample=list(SAMPLES.index))
        input_list += expand(RESULT_DIR / "{sample}/FinalSearch/proteomes_Default_PSM_Report.txt", sample=list(SAMPLES.index))
        
        
        #input_list.append(RESULT_DIR / "taxidMapping/accessions_hashed.npy")

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
 
        