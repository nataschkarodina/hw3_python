# Homework 3
# created by Natalia Rodina
# 08.04.2018

#packages
from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd

#The class Kmer_spectrum
class Kmer_spectrum:
    k = 0 #kmer size
    q = 0 #quality
    kmer_dict = {} #dictionary of kmers
    kmer_list = [] #list of kmers
    dict_spectra = {} #dictionary of kmer distribution
    genome_size = 0 #the size of genome

    def __init__(self, filename, k, q):
        self.k = k
        self.q = q
        self.filename = filename

    #2.1 creates the dictionary of all kmers with set length kmer_size and quality higher then quality_for_nucleotide
    def kmer_dictionary (self):
        sequences = SeqIO.parse(self.filename, 'fastq')
        for record in sequences:
            seq_len = len(record.seq)
            for index in range(seq_len - self.k + 1):
                current_kmer = record.seq[index:(index + self.k)] #finds current kmer
                current_quality = record.letter_annotations["phred_quality"][index:(index + self.k)] #finds the current kmer quality

                #checks if the quality (v) in current_quality is not lower then q
                if any(v < self.q for v in current_quality):
                    pass
                # checks if the kmer is in dictionary (and has the right quality) -
                else:
                    if current_kmer in self.kmer_dict:
                        self.kmer_dict[current_kmer]+=1
                    else:
                        self.kmer_dict[current_kmer]=1
        return (self.kmer_dict)

    # 2.2 creates a dictionary of kmer distribution from frequenct and prints it as a table
    def kmer_spectra(self):
        kmer_list_y = []
        kmer_list_x = []
        #checks if the current kmer is already in dict_spectra and writes the frequency
        for key in self.kmer_dict:
            if self.kmer_dict[key] in self.dict_spectra:
                self.dict_spectra[self.kmer_dict[key]] += 1
            else:
                self.dict_spectra[self.kmer_dict[key]] = 1

         #writes the kmers and their frequency to a list
        for key in self.dict_spectra:
            self.kmer_list += [[key, self.dict_spectra[key]]] #a list of kmer count and multiplicity
            kmer_list_x.append(key) #list of kmer multiplicity
            kmer_list_y.append(self.dict_spectra[key]) #list of kmer count
        output = pd.DataFrame(data={'kmer multiplicity': kmer_list_x, 'kmer count': kmer_list_y}) #dataframe of kmer spectrum
        print(output) #prints the dataframe
        return (self.kmer_list)

    # 2.3 plots the kmer count from multiplicity as a bar plot - takes as input the cut value
    def plot_spectra (self, name, cutoff, xlim, ylim):
        for key in self.dict_spectra:
            x_x = key
            y_y = self.dict_spectra[key]
            plt.bar(x_x, y_y, color = 'green')
        #max_y = max(self.dict_spectra.values()) + 1
        #max_x = max(self.dict_spectra.keys()) + 1
        plt.bar(cutoff, ylim, color = 'red') #plots the cutoff
        plt.xlabel('Diversity of found kmers')
        plt.ylabel('The frequency of kmers')
        plt.title('The spectrum of kmers')
        plt.axis([0, xlim, 0, ylim])
        plt.savefig(name)
        plt.show()
        plt.close()

    # Calculates the approximate genome length - takes as input the cut value
    def genome_length(self, cut):
        self.genome_size = 0

        #determine the peak position
        x_max = 0 #maximum by x
        y_max = 0 #maximum by y
        for key,value in self.dict_spectra.items():
            if (value > y_max) and (key > cut):
                x_max = key
                y_max = value
        print ('x_max', x_max)

        #determine the genome size by summarizing all compositions of key*value in dict_spectra and deviding it by the peak position
        for key, value in self.dict_spectra.items():
                if key > cut:
                    self.genome_size += key * value
        self.genome_size = round(self.genome_size/x_max)
        return (self.genome_size)

# RUN THE SKRIPT
#NOTE: the fastq file should be in the same directory as the skript!

#WITH QUALITY OF NUCLEOTIDE READING CUTOFF
#SET PARAMETERS
kmer_size = 11 #the size of kmer
quality_for_nucleotide = 30 #quality of single nucleotide reading setup
cut = 20 #set the cutoff parameter

kmers = Kmer_spectrum('test_kmer.fastq', kmer_size, quality_for_nucleotide)
#USES THE FUNCTION KMER_DICTIONARY -
kmers.kmer_dictionary()
#USES THE FUNCTION KMER_SPECTRA-
kmers_list = kmers.kmer_spectra()
#USES THE FUNCTION PLOT_SPECTRA -
x_lim = 200
y_lim = 40000
kmers.plot_spectra('q30_3', cut, x_lim, y_lim)
#USES FUNCTION GENOME_LENGTH -
genome_size = kmers.genome_length(cut)
print('Genome size,'
      'kmer size =', kmer_size,
      'and quality of nucleotide reading', quality_for_nucleotide)
print(genome_size)

#WITHOUT QUALITY OF NUCLEOTIDE READING CUTOFF
#SET PARAMETERS
kmer_size = 11 #the size of kmer
quality_for_nucleotide = 0 #quality of single nucleotide reading setup
cut = 50 #set the cutoff parameter

kmers2 = Kmer_spectrum('test_kmer.fastq', kmer_size, quality_for_nucleotide)
#USES THE FUNCTION KMER_DICTIONARY -
kmers2.kmer_dictionary()
#USES THE FUNCTION KMER_SPECTRA-
kmers_list2 = kmers2.kmer_spectra()
#USES THE FUNCTION PLOT_SPECTRA -
x_lim = 300
y_lim = 40000
kmers2.plot_spectra('q0_3', cut, x_lim, y_lim)
#USES FUNCTION GENOME_LENGTH -
genome_size_2 = kmers2.genome_length(cut)
print('Genome size,'
      'kmer size =', kmer_size,
      'and quality of nucleotide reading', quality_for_nucleotide)
print(genome_size_2)





