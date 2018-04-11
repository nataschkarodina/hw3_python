# Homework 3

This skript finds all k-mers with set length, creates and plots the kmers spectrum and calculates the approximate genome length

# Usinfg the skript 

You will need to download the skipt hw3_python_Rodina_test.py

## Input file
You will need a fastq sequence containing contigs
If you use your own sequence file don't forget to change the name of file in the skript!

Change the current file name in #110 to that of yours. Remember, the file shoud be in fastq format
```
kmers = Kmer_spectrum('test_kmer.fastq', kmer_size, quality_for_nucleotide)
```
## Setting parameters

You will need to set up the parameters:
kmer_size is the size of the kmers 
quality_for_nucleotide is the value of the nucleotide reading quality
cut is the cutoff parameter

Just change it in the skript:
```
kmer_size = 11 #the size of kmer
quality_for_nucleotide = 30 #quality of single nucleotide reading setup
cut = 20 #set the cutoff parameter
```
You will also need to change the parameters of the plot:
x_lim and y_lim are limits for the size of the plot
You should also set the name for the plot - the plot will be saved in png format

```
x_lim = 200
y_lim = 40000
kmers.plot_spectra('q30_3', cut, x_lim, y_lim)
```
## The skript functions:
The skript contains the class Kmer_spectrum that takes for input the name of the fastq file, 
kmersize, quality for nucleotide reading
```
kmers = Kmer_spectrum('test_kmer.fastq', kmer_size, quality_for_nucleotide)
```

### Function kmer_dictionary:
creates the dictionary of all kmers with set length kmer_size and quality higher then quality_for_nucleotide:

```
kmers.kmer_dictionary()
```
### Function kmer_spectra:
creates a dictionary of kmer distribution from frequenct and prints it as a table
```
kmers_list = kmers2.kmer_spectra()
```
### Function plot_spectra:
plots the kmer count from multiplicity as a bar plot - takes as input the name of the plot for saving, the cutoff value, 
and the limits for the x and y axis
```
x_lim = 300
y_lim = 40000
kmers.plot_spectra('q0_3', cut, x_lim, y_lim)
```
### Function genome_length 
Calculates the approximate genome length - takes as input the cut value
```
genome_size = kmers.genome_length(cut)
```
## What is the output?

The skript will provide you a table of the kmers distributaion. And plot the kmers spectra.
According to the set cuf off parameter it will calculate the approximate genome length

# Report
The current skipt will run the class two times with different parameters:
## WITH QUALITY OF NUCLEOTIDE READING CUTOFF

The set parameters:
```
kmer_size = 11 #the size of kmer
quality_for_nucleotide = 30 #quality of single nucleotide reading setup
cut = 20 #set the cutoff parameter
x_lim = 200
y_lim = 40000
```
* Hope you liked the skript
* And found it useful
* Good luck!
