# MultiStageSearch

## Goal of this Pipeline
In Proteomic analysis, Viruses often yield low identification rates, especially on strain level.
The aim is the usage of a proteogenomic approach in a multistage search approach to overcome this problem.

It consists roughly of the following steps:

1. Host and Contaminant Filtering (optional)
2. Search against a reference proteome database provided by the user
3. Mapping and scoring of PSMs to find interesting species
4. Download of genomic data by one of the following approaches:
    1. Search for sequences linked to the species taxon
    2. Infer descendants of the given Taxon and search for sequences according to the descendant Taxon name
5. Translation of the nucleotide sequences to proteomes using a 6-frame translation approach
6. Concatenation of the proteomes for the sample
7. Reduction of the search space by combining duplicate sequences
8. Search agains the created proteome database of translates nucleotide sequences
9. Combining of the results of both search steps for potential downstream analysis
10. Creation of different plots to evaluate the results

Often, strain level genomes are only linked to the species taxon.
Therefore a workaround is required to fetch sequences for strain of the given taxon.
One way is just download sequences linked to the species taxon and infer the strain from the metadata of the sequence records.
An other way is to search sequences with similar names to the strain name. In Theory, this approach should yield a more precise result.
But sometimes also the names of strains found in the taxonomy are not equal to sequence name records in the NCBI Database.
If using the descendants approach, when no sequence is found, the workflow automatically falls back to the other approach for the affected taxon id.
The user should test both approaches for the given sample.

The workflow was developed in a way that allows the analysis of multiple samples at the same time.

## Input:
- samples.tsv: A tab separated list with the following information for each sample:
    - sample: Name of the sample
    - mgf: Path to the mgf file
    - par: Path to the par file
    - host_fasta: Path to the fasta file containing the host database. Please make sure that the file extension is "fasta". 
    - host_name: Scientific Name of the Host
- Additionally a fasta file containing common contaminants is required for the optional host filtering step.
- Please make sure to adjust the config.yaml according to your paths and the NCBI account data (Mail and API key). 

## Getting started:
Please make sure that 
[SearchGUI](http://compomics.github.io/projects/searchgui) 
and 
[PeptideShaker](http://compomics.github.io/projects/peptide-shaker) 
are installed and the paths are configured in the config file.\
Other dependencies (other than the paths in the config) should be handled by the workflow automatically.\
The workflow was tested using Ubuntu 22.04LTS.

Make sure that snakemake is installed by using conda. A tutorial for that can be found [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).\
The workflow can be started with the following command while being in the MultistageSearch directory:
```
snakemake --use-conda -c <n_cores> 
```

## Output:
- Plots showing statistics of the results
- a combined Peptide-shaker report of both search iterations