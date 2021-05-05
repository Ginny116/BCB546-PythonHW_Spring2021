#!/usr/bin/env python
# coding: utf-8

# # Python_Assignment

# In[2]:


from Bio import SeqIO

from Bio.Data import CodonTable

import pandas as pd


# ## Function

# 1. Document Dr. X's function with comments and with markdown text in your Jupyter notebook

# In[3]:


def get_sequences_from_file(fasta_fn): # define a function to get sequences from fasta files
    sequence_data_dict = {} # create a dictionary
    for record in SeqIO.parse(fasta_fn, "fasta"): # start a for-loop over the records in the fasta file 
        description = record.description.split() # extract the unique elements
        species_name = description[1] + " " + description[2] # create species_name by extracting element 1 and space and then element 2 in discription
        sequence_data_dict[species_name] = record.seq # the dictionary contains "species_name" as key and sequence as value 
    return(sequence_data_dict) # reture to the sequense data dictionary


# In[5]:


seq = get_sequences_from_file("bears_cytb.fasta")


# In[6]:


from Bio.Data import CodonTable # import Codon table from Bio.Data
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"] # get the codon table and name it mito_table
print(mito_table)


# 2. Write a function that translates a string of nucleotides to amino acids based on Dr. X's pseudo-code suggestion.

# In[7]:


def translate_function(string_nucleotides):
    mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"] # get codon table and name mito_table
    aa_seq_string = "" # set an empty string-type of "aa_seq_string"
    i = 0 # set value of element i to 0
    while i < len(string_nucleotides): # A while-loop goes on until the value of "i" is equal the length of the sequence defined by string_nucleotides
        codon = string_nucleotides[i:i+3] # Create "codon" containing string_nucleotides. A codon contains 3 nucleotides (from i to i+3)
        if (codon == 'TAG' or codon == 'TAA' or codon == 'AGA' or codon == 'AGG'):
            break  # stop the loop if a codon matches to one of the known "stop" codons
        aa = mito_table.forward_table[codon] # get the aa in mito_table that matches its codon
        aa_seq_string += aa # increase the aa_seq_string by adding aa
        i += 3 # value of element i is inceased by 3
    return(aa_seq_string)


# In[11]:


translate_function("ATGACCAACATCCGAAAA") # test the translate_function


# 3. Write an alternative translation function

# In[52]:



def translate_function_alt(nucleotides): # define translate_funtion_alt
    from Bio.Seq import Seq # import seq from Bio.Seq package
    seq = Seq(nucleotides)  # create seq sequence 
    aa_seq = seq.translate(table=2, to_stop=True) # translate seq by using transcription table and stop at stop codon
    return(str(aa_seq))


# In[53]:


translate_function_alt("ATGACCAACATCCGAAAA") # test translation_function_alt


# 4. Write a function that calculates the molecular weight of each amino acid sequence.

# In[38]:


from Bio.SeqUtils.ProtParam import ProteinAnalysis
def cal_molecular_weight(aa_seq):
    aa_seq_Analysed = ProteinAnalysis(aa_seq) # analyses sequence
    molecular_weight = aa_seq_Analysed.molecular_weight() # using molecular weight attribute to calculate weight
    return(molecular_weight)  # returning the molecular weight


# In[39]:


cal_molecular_weight("MTNIRK") # test cal_molecular_weight function


# 5. Write a function that computes the GC-content of each DNA sequence.

# In[54]:


from Bio.SeqUtils import GC # import the GC function from Bio.Sequtils
def gc_content(nucleotides): # define gc_content function
    gc_content = GC(nucleotides) # calculate the GC content
    return(gc_content) 


# In[55]:


gc_content("ATGACCAACATCCGAAAA") #  test gc_content function


# ## MAIN part of the script

# In[42]:


cytb_seqs = get_sequences_from_file("bears_cytb.fasta") 
bears_df = pd.read_csv("bears_mass.csv")
species_list = list(bears_df.species)


# 6. Add two new columns to the bears DataFrame: (1) molecular weight and (2) GC content

# In[56]:


import numpy as np
import pandas as pd
bears_df["molecular_weight"]='Non' # in bears_df data frame add columns "molecular_weight" and set value as Non
bears_df["GC_content"]='Non' # in bears_df data frame add columns "GC_content" and set value as Non
bears_df 


# 7. Write a for-loop that translates each sequence and also gets molecular weight and computes the GC content of each translated sequence and adds those data to DataFrame

# In[45]:


molecular_weight = []
gc_content = []
    
for key, value in cytb_seqs.items(): #A for-loop through all elements of the cytb sequence
    DNA_seq = (str(value)) # get the DNA_seq 
    aa_seq = translate_function(DNA_seq) # translate DNA_seq and name as aa_seq
    
    molecular_weight.append(cal_molecular_weight(aa_seq)) # calculate the molecular weight 
    gc_content.append(gc_content1(DNA_seq)) # calculate the GC_contents

bears_df['molecular_weight'] = molecular_weight # add the molecular_weight into the molecular weight column in the bears_df data frame
bears_df['GC_content'] = gc_content  # add the GC_content into the molecular weight column in the bears_df data frame

print(bears_df)


# 8. Plot a bar-chart of the mass with the x-axes labeled with species names.

# In[48]:


import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

plot_dims = (16, 5)
sns.set(font_scale = 1.5)
fig, ax = plt.subplots(figsize=plot_dims)
sns.barplot(x = 'species', y = 'mass', data = bears_df)
ax.set(xlabel = 'Species ID', ylabel = 'Mass')
plt.xticks(rotation=90)


# *Q1* What is the largest bear species? 
# Ursus Spelaeus
# 
# *Q2* What else is interesting about this species?
# Ursus Spelaeus are found mainly in caves, therefore, they are called cave bears

# 9. Plot a visualization of the molecular weight (y-axis) as a function of GC-content (x-axis).

# In[50]:


import matplotlib.pyplot as plt
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')
sns.lmplot( 'GC_content','molecular_weight', data = bears_df, fit_reg=False, height=7, aspect=2.0, hue="species")


# 10. Save the new DataFrame to a file called "bears_mass_cytb.csv"

# In[51]:


bears_df.to_csv("bears_mass_cytb.csv", index = False)


# In[ ]:




