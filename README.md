# MTBDiff

<img align="right" src=https://github.com/dmnfarrell/mtbdiff/blob/master/img/logo.png width=250px>

Mycobacterium tuberculosis is the causative agent of human TB. Not only are there are multiple closely related lineages of the bacteria present in human populations, there are also animal adapted organisms such as M. bovis that infect a broad range of animal species beyond their most prominent host in cattle. The entire group is referred to as the M. tuberculosis complex (MTC). The animal adapted species are characterized by high sequence similarity to MTB but with key sequence polymorphisms called regions of difference in which deletions or insertions are present. Some of these may be important for adaption to the chosen niche of the species.

With the increased use of whole genome sequencing it is important to be able to compare new and existing isolates systematically and rapidly. MTBDiff performs quick analysis of structural and SNP changes compare to the laboratory reference species MTB H37Rv. It outputs a table indicating of previously known regions (RDs) are present in the species. Multiple genomes can be compared at once. It also indicates if any novel structural changes are present.

This tool uses MUMmer for aligning genomic sequences and NucDiff which parses the results. NucDiff locates and categorizes differences between two closely related nucleotide sequences. It is able to deal with very fragmented genomes, structural rearrangements and various local differences. It is also suitable for assemblies.

## Installation

MUMmer is required. If using and Ubuntu based distro use apt to install MUMmer:
```sudo apt install mummer```

Then:
```pip install mtbdiff```

## Usage

Put your genome files in fasta format inside a folder and then run:

```mtbdiff -i <folder>```

See the wiki for more details.

## References

* Brites, Daniela et al. “A New Phylogenetic Framework for the Animal-Adapted Mycobacterium tuberculosis Complex.” Frontiers in microbiology vol. 9 2820. 27 Nov. 2018, doi:10.3389/fmicb.2018.02820
* Faksri, Kiatichai et al. “In silico region of difference (RD) analysis of Mycobacterium tuberculosis complex from sequence reads using RD-Analyzer.” BMC genomics vol. 17,1 847. 2 Nov. 2016, doi:10.1186/s12864-016-3213-1
* R. Brosch, S. V. Gordon, M. Marmiesse, P. Brodin, C. Buchrieser, K. Eiglmeier, T. Garnier, C. Gutierrez, G. Hewinson, K. Kremer, L. M. Parsons, A. S. Pym, S. Samper, D. van Soolingen, S. T. Cole. A new evolutionary scenario for the Mycobacterium tuberculosis complex. Proceedings of the National Academy of Sciences Mar 2002, 99 (6) 3684-3689; DOI: 10.1073/pnas.052548299
* Mostowy S, Inwald J, Gordon S, et al. Revisiting the evolution of Mycobacterium bovis. J Bacteriol. 2005;187(18):6386–6395. doi:10.1128/JB.187.18.6386-6395.2005