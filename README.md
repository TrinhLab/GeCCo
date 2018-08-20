# Gene expression classifier
# Installation
1. Clone the repository by running the following command (from [Git Bash in Windows](https://git-scm.com/downloads) or the terminal in Linux):
    ~~~
    git clone https://github.com/trinhlab/gec.git
    ~~~
2. Install the python package (for Windows the [Anaconda Prompt](https://www.anaconda.com/distribution/) is recommended):
    ~~~
    pip install gec
    ~~~
 
# Usage
1. Create a _problem directory_ with the following files in an input subdirectory:
    - __tpm.csv__ (mandatory): Must contain the following headers: `Gene|WT_t1_rep1|WT_t1_rep2|WT_t1_rep3|WT_t2_rep1|WT_t2_rep2|MT_t1_rep1|MT_t1_rep2|MT_t2_rep1|MT_t2_rep2` Any number of replicates is acceptable. 
    - __gene_features.csv__(optional): Headers are `Gene|Feature`. Examples of features are mutated_gene or transcription_factor
    - __coexpression_tpm.csv__ (optional): Used to construct co-expression network, First column must be labeled 'Gene', every subsequent column will be treated as a sample. 
     
    For example a problem directory named "p1" should have the following structure: 
    ~~~
    \---p1
        |   README.md
        \---input
                coexpression_tpm.csv
                gene_features.csv
                tpm.csv
    ~~~
    
    :warning: The transcript per million (tpm) data provided in any of the input files must not be log transformed. 
2. Open the command line(for Windows the Anacondas Prompt is recommended) and execute:  
    ~~~
    gec <problem_directory>
    ~~~
    where `<problem_directory>` is the path to the problem you would like to run. For example, to run problem p1 execute:
    ~~~
    gec gec\problems\p1
    ~~~
    Then your output figures and cytoscape input will be saved in your problem directory. 
    
    To explore additional options execute
    ~~~
    gec --help
    ~~~
