# sRNA-jbrowse
 A set of processing scripts for setting up jbrowse with ShortStack sRNA-data


### Motivation

Jbrowse (1) is an extremely powerful and modifiable genome browser. However, it can also be quite opaque how to set it up - especially with unique data (like sRNA seq clusters and alignment).

This repo is meant to be a set of directions to process sRNA-seq alignments, annotations, and genomic annotations to produce a high-quality jbrowse environment.


### Prerequisites
The following software must be runable from command line.


```python3```

```gt``` (genometools, http://genometools.org/) - this can also be installed easily through package managers (brew install genometools for mac)

```bgzip``` and ```tabix``` - both available through the htslib package http://www.htslib.org/download/

```ShortStack``` - available from https://github.com/MikeAxtell/ShortStack, with additional required packages.


### General pipeline



#### Reformatting ShortStack gff

The general output of all annotations from a ShortStack run is a gff3 file. Several modifications need to be performed to make it work for jbrowse. We first add a header giving ```sequence-region``` identifiers. We also remove gff description keys that start with a capital letter.

```python cleanup_gff.py -g ShortStack_All.gff3 ```

Subsequent steps sort, block-zip, and index this gff.

```
gt gff3 -retainids -sortlines -tidy ShortStack_All.clean.gff3 > ShortStack_All.clean.sorted.gff3
bgzip ShortStack_All.clean.sorted.gff3 
tabix -f -p gff ShortStack_All.clean.sorted.gff3.gz 

```

All of these steps should hopefully function using the single python script:

```sh prepare_gff.sh  ```



