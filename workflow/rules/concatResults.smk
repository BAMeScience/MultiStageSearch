# TODO: FDR beachten
rule catDBResults:
     input: 
          firstDBSearch = RESULT_DIR / "{sample}/FirstSearch/ref_Default_PSM_Report.txt",
          finalDBSearch = RESULT_DIR / "{sample}/FinalSearch/proteomes_Default_PSM_Report.txt",
     output: 
          concatResults = RESULT_DIR / "{sample}/refSeqViral_Default_PSM_Report.txt" # TODO refseq als variable
     shell: 
          "cat {input.firstDBSearch} {input.finalDBSearch} > {output.concatResults}"