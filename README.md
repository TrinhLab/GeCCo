# Gene expression classifier (GEC)

## What is GEC
(TODO:expand on this explanation)

GEC serves to analyze omics data between two conditions or time points. Such comparison is referred to as case vs control. GEC will determine genes or proteins that fall under different categories:
- highly-expressed
- lowly-expressed
- up-regulated
- down-regulated
- changed-regulation
- no-change
Additionally, GEC will perform network co-expression analysis and produce several files ready for analysis in [Cytoscape](https://cytoscape.org/).

## Preliminary notes
GEC is a python program with a command line interface so that you can get things done efficiently. The commands below are intended for a unix-like OS (MacOS, *BSD, GNU/Linux), if you are on Windows you can use a variety of emulation options (e.g., Windows Subsystem for Linux, Cygwin, Virtualbox) to reproduce a linux command line environment or simply use the Windows command line, which may require some adjustments (If you choose to use Windows directly check out how to run python programs, one option is the [Anaconda Prompt](https://www.anaconda.com/distribution/))

## Installation
1. Clone the repository or download the zip file.
2. Install the python package: `pip install -e gec`

## Usage
1. Create a _problem directory_ with the following files in an input subdirectory:
    - __tpm.csv__ (mandatory): Must contain the following headers: `Gene|WT_t1_rep1|WT_t1_rep2|WT_t1_rep3|WT_t2_rep1|WT_t2_rep2|MT_t1_rep1|MT_t1_rep2|MT_t2_rep1|MT_t2_rep2` Any number of replicates is acceptable. Note that WT corresponds to _case_ and MT to _control_.
    - __gene_features.csv__(optional): Headers are `Gene|Feature`. Any arbitrary feature name may be used, for example  mutated_gene or transcription_factor.
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

    :warning: The transcript per million (tpm) data provided in any of the input files must not be log transformed. :warning:

    The file `header_map.csv` keeps tracks of the original data headers but it is currently not used by GEC.

2. Open the command line(for Windows the Anacondas Prompt is recommended) and execute:
    ~~~
    gec <problem_directory>
    ~~~
    where `<problem_directory>` is the path to the problem you would like to run. For example, to run problem p1 execute:
    ~~~
    gec gec/problems/p1
    ~~~
    Then your output figures and cytoscape input will be saved in your problem directory.

    To explore additional options execute
    ~~~
    gec --help
    ~~~
