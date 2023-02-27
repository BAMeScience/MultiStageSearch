class General:
    @staticmethod
    def get_input_all(wildcards):
        input_list = []

        # host filtering
        # if not SKIP_HOST_FILTERING:
        #     input_list += expand(RESULT_DIR / "{sample}/SpectraFilter/Filtered_host.mgf", sample=list(SAMPLES.index))
        
        input_list += expand(RESULT_DIR / "{sample}/FirstSearch/ref_Default_PSM_Report.txt", sample=list(SAMPLES.index))

        return input_list

# TODO: Not used yet
class SearchDB:
    @staticmethod
    def get_input_SearchDB():
        if not SKIP_HOST_FILTERING:
            data_dict = {
                "mfg": RESULT_DIR / "{sample}/SpectraFilter/Filtered_host.mgf"
                }
        else:
            data_dict = {
                "mfg": MGF_FILE,
                }
        return data_dict