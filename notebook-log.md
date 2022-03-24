# 2/17
## object of study: fliY gene of bacillus family
## database: UniProt

### Data already have:
1) Bacillus beveridgei
2) bacillus horikoshii
3) bacillus altitudines
4) bacillus stratosphericus
5) bartpnella rochalimae

### QC software: FastQC

Problem: still confuse of how to run FastQC

# 2/18
## try to use clustalw
### code:
grep ">" fily.fasta
~/desktop/software/clustalw2 -ALIGN -INFILE=fily.fasta -OUTFILE=fily-aligned.fasta -OUTPUT=PHYLIP

### result:
CLUSTAL 2.1 Multiple Sequence Alignments


Sequence format is Pearson
Sequence 1: tr|I2C561|I2C561_BACAM           378 aa
Sequence 2: tr|A0A1D7QVZ1|A0A1D7QVZ1_9BACI   418 aa
Sequence 3: tr|E6YKZ9|E6YKZ9_9HYPH           170 aa
Sequence 4: tr|A0A1Y0CMP9|A0A1Y0CMP9_9BACI   410 aa
Sequence 5: tr|A0A5K1NDT8|A0A5K1NDT8_9BACI   373 aa
Sequence 6: tr|A0A4Y9K066|A0A4Y9K066_BACIT   373 aa
Start of Pairwise alignments
Aligning...

Sequences (1:2) Aligned. Score:  57
Sequences (1:3) Aligned. Score:  27
Sequences (1:4) Aligned. Score:  59
Sequences (1:5) Aligned. Score:  72
Sequences (1:6) Aligned. Score:  73
Sequences (2:3) Aligned. Score:  28
Sequences (2:4) Aligned. Score:  53
Sequences (2:5) Aligned. Score:  59
Sequences (2:6) Aligned. Score:  59
Sequences (3:4) Aligned. Score:  21
Sequences (3:5) Aligned. Score:  24
Sequences (3:6) Aligned. Score:  24
Sequences (4:5) Aligned. Score:  60
Sequences (4:6) Aligned. Score:  60
Sequences (5:6) Aligned. Score:  94
Guide tree file created:   [fily.dnd]

There are 5 groups
Start of Multiple Alignment

Aligning...
Group 1: Sequences:   2      Score:7777
Group 2: Sequences:   3      Score:6933
Group 3: Sequences:   4      Score:6428
Group 4: Sequences:   5      Score:6421
Group 5:                     Delayed
Alignment Score 14891

[fily-aligned.fasta]

# 2/23
## Try to use MUSCLE

### first try
code:
$ pwd
$ cd desktop/data/fily
$ ~/desktop/software/muscle5.1.win64.exe -in fily.fasta -out fily-muscle-aligned.fasta

### result
false: unknown option in

# 3/2
## HW Aligning your own data
### Aligning method/software choice: clustalw
### Assumptions:
1) sequences are down-weighted compared to how closely related they are to other sequences
2) the program varies the gap penalties (GP) for sequences and positions
### Limitations:
1) The guide tree has a big impact on alignments
2) Errors made early in the process persist since subsequent mergers never change the alignments they are merging together
Process and code: see 2/18

# 3/4
## Make a tree
### change the phylip file into a fasta form
code:
$ grep ">" filY.fasta
$ ~/desktop/software/clustalw2 -ALIGN -INFILE=fily.fasta -OUTFILE=fily-aligned.fasta -OUTPUT=FASTA
### use R to make a tree
code:
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep = TRUE)
``` (Install necaccary packages)

```{r}
library(ape)
library(adegenet)
library(phangorn)

seq <- read.FASTA(file="../filY/fily-aligned.fasta",type = "AA")
``` (Read the sequence as AA)

```{r}
D <-  dist.ml(seq)
tre <- nj(D)
tre <-ladderize (tre)

plot (tre, cex =.6)
title ("A simple NJ tree")
``` (Make a simple evolution tree)
```

# 3/9
## Align the formal data (FliY protein sequences of 52 subtype of Bacillus)
### Software: Git Bash & Clustalw
code:
$ pwd
$ cd desktop/data/fily
$ grep ">" filY.fasta
$ ~/desktop/software/clustalw2 -ALIGN -INFILE=fily.fasta -OUTFILE=filY-aligned.fasta -OUTPUT=FASTA

## Make a brief evolution tree
### Software: R
code:
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep = TRUE)
``` (Install necaccary packages)

```{r}
library(ape)
library(adegenet)
library(phangorn)

seq <- read.FASTA(file="../fliY/fliY-aligned.fasta",type = "AA")
``` (Read the sequence as AA)

```{r}
D <-  dist.ml(seq)
tre <- nj(D)
tre <-ladderize (tre)

plot (tre, cex =.6)
title ("A simple NJ tree")
``` (Make a simple evolution tree)
```
# 3/23
Do nothing this week :(
So many midterms
