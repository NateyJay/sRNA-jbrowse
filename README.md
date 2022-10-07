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






