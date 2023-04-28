rule catDBResults:
     input: 
          firstDBSearch = RESULT_DIR / "{sample}/FirstSearch/ref_Default_PSM_Report.txt",
          finalDBSearch = RESULT_DIR / "{sample}/FinalSearch/proteomes_Default_PSM_Report.txt",
     output: 
          concatResults = RESULT_DIR / "{sample}/concat_PSM_Report.txt"
     log:
        stderr_log=RESULT_DIR / "logs/concatResults/{sample}/stderr.log",
        stdout_log=RESULT_DIR / "logs/concatResults/{sample}/stdout.log"
     conda:
        "../envs/pandas.yml"
     threads: 1
     script: 
          "../scripts/concat_results.py"