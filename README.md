## README

### Introduction

This software is used for calculating copy number of genes, and get top 5% genes which selected by humans.

###  Dependencies

Python modules:  
&ensp;&ensp;&ensp;&ensp;numpy  
&ensp;&ensp;&ensp;&ensp;scipy.stats  
Software:  
&ensp;&ensp;&ensp;&ensp;mosdepth  
&ensp;&ensp;&ensp;&ensp;samtools  

### Installation

```sh
cd /path/to/install
git clone https://github.com/sc-zhang/popCNV.git
cd popCNV/bin
chmod +x *.py
echo 'export PATH=/path/to/install/popCNV/bin:$PATH' >> ~/.bash_profile
source ~/.bash_profile
```

### Usage

```sh
popCNV.py -g genome.fasta -s 1000 -r read_depth/ -b bam_files/ -l gene.list -w wrk_dir --group group.list --sample sample.group --wild 0
```

**-g, --genome** fasta file of genome

**-s** window size, must same to the parameter which used in mosdepth

**-r** read_depth, a folder which contains results generated by mosdepth

**-b** bam_file , a folder which contains all bam files and index files or a file contains two columns: sample_name sequence_depth

**Notice: Usage of mosdepth**

```sh
mkdir read_depth
cd read_depth
mosdepth -b 1000 sample_name /path/to/bam/sample_name.bam
```

**Notice: the bam file must be named start with sample_name and a dot, the prefix used in mosdepth must be sample_name**

**-l** gene.list, a list file contains 4 columns: chromosome_name start_position end_position gene_name

**-w** wrk_dir, a directory for working

**--group** group.list, a file contains two columns: group_id group_name

**--sample** sample.list, a file contains two columns: sample_name group_id

**Notice: there must be no spaces in names of groups and names of samples**

**--wild** 0, the group_id of the group which was used as wild group

**-t, --threads** threads

**--rerun** if this parameter was used, the software will restart all steps

### Results

**06.genes.round.cn** Copy number of genes in all groups

**07.RFD**

RFD is relative frequency difference which was calculated based on wild group

a) a set of list files named with group_name, contain a list of genes with top 5% of RFD based on p<0.05

b) Total.list contains all RFD and p-value of genes in all groups