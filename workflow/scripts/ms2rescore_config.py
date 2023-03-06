out_file = snakemake.output[0]

RescorePipeline = snakemake.params[0]
RescoreFeatures = snakemake.params[0]
RunPercolator = snakemake.params[0]
FragModel = snakemake.params[0]

with open(out_file, "w") as f:
    
    lines = [' {\"$schema\":\"./config_schema.json\",','\"general\":{']
    lines.append(f'"pipeline":{RescorePipeline},')
    lines.append(f'"feature_sets":{RescoreFeatures},')
    lines.append(f'"run_percolator":{RunPercolator},')
    lines.append('"id_decoy_pattern": "DECOY",')
    lines.append('"num_cpu":'+'1'+',')
    lines.append('"config_file":null,')
    lines.append('"tmp_path":null,')
    lines.append('"mgf_path":null,')
    lines.append('"output_filename":null,')
    lines.append('"log_level":\"info\",')
    lines.append('"plotting":false')
    lines.append('},')
    lines.append('\"ms2pip\":{')
    lines.append(f'\"model\":{FragModel},')
    lines.append('\"modifications\":[')
    # if mods:
    #     for mod in mods:
    #     lines.append(InputMod(mod[0],mod[1],mod[2],mod[3],mod[4],mod[5]))
    #lines[-1] = lines[-1][:-1]

    lines.append(']')
    lines.append('},')
    lines.append('\"maxquant_to_rescore\":{},')
    lines.append('\"percolator\":{}')
    lines.append('}')

    f.writelines([line + "\n" for line in lines])