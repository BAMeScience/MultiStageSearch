rule createPlots:
    input:
        strain_mappings = RESULT_DIR / "{sample}/taxids/strain_name_counts.tsv",
        taxID_scores = RESULT_DIR / "{sample}/taxids/taxid_scores.tsv",
        frist_search = RESULT_DIR / "{sample}/FirstSearch/ref_Default_PSM_Report.txt",
        final_search = RESULT_DIR / "{sample}/FinalSearch/proteomes_Default_PSM_Report.txt",
    output:
        strain_bar_plot = RESULT_DIR / "{sample}/Plots/{sample}_strain_bar_plot.png",
        taxIdScores_bar_plot = RESULT_DIR / "{sample}/Plots/{sample}_taxIdScores_bar_plot.png",
        proportions_pie_chart = RESULT_DIR / "{sample}/Plots/{sample}_proportions_pie_chart.png",
        first_search_confidence_histogram = RESULT_DIR / "{sample}/Plots/{sample}_first_search_confidence_histogram.png",
        final_search_confidence_histogram = RESULT_DIR / "{sample}/Plots/{sample}_final_search_confidence_histogram.png",
    log:
        stderr_log=RESULT_DIR / "logs/Plots/{sample}/stderr.log",
        stdout_log=RESULT_DIR / "logs/Plots/{sample}/stdout.log"
    conda:
        "../envs/plots.yml"
    threads: 1
    script:
        "../scripts/create_plots.py"