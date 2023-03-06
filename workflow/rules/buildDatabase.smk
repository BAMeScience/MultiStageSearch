rule catDatabase:
     input: 
          protacc2taxids00 = "resources/Database/taxidMapping/protacc2taxids00",
          protacc2taxids01 = "resources/Database/taxidMapping/protacc2taxids01",
          protacc2taxids02 = "resources/Database/taxidMapping/protacc2taxids02"
     output: 
          protacc2taxids_virus = RESULT_DIR / "taxidMapping/protacc2taxids_virus.txt"
     shell: 
          "cat {input.protacc2taxids00} {input.protacc2taxids01} {input.protacc2taxids02} > {output.protacc2taxids_virus}"


rule splitToAccessions:
     input: 
          protacc2taxids_virus = RESULT_DIR / "taxidMapping/protacc2taxids_virus.txt"
     output: 
          accessions = RESULT_DIR / "taxidMapping/accessions.txt"
     shell:
        "awk '{{print $1}}' {input.protacc2taxids_virus} > {output.accessions}"


rule splitToTaxids:
     input:
          protacc2taxids_virus = RESULT_DIR / "taxidMapping/protacc2taxids_virus.txt"
     output:
          taxids = RESULT_DIR / "taxidMapping/taxids.txt"
     shell:
          "awk '{{print $2}}' {input.protacc2taxids_virus} > {output.taxids}"


rule hashDatabase:
     input:
          accessions = RESULT_DIR / "taxidMapping/accessions.txt"
     output:
          accessions_hashed = RESULT_DIR / "taxidMapping/accessions_hashed.npy"
     conda: 
          "../envs/hashing.yml"
     script:
          "../scripts/hashDatabase.py"