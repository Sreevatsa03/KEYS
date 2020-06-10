# What can we learn from a DNA sequence?

**Overview of the hands-on exercise**

This hands-on workshop focuses on an unknown DNA sequence collected in a surface seawater metagenomic sample. The sequence is available in contig_unknown.fasta file.

**Discussion: What can we learn from an unknown DNA sequence?**

## 1- What kind of organism this sequence belongs to?
The first question that comes in mind when dealing with an unknown sequence is trying to identify the "who"?

### 1.A Run a BlastN query against the NCBI GenBank database

Step1: Go to the [NCBI Blast page](https://blast.ncbi.nlm.nih.gov/Blast.cgi) and select the Nucleotide Blast tool.

Step2: Paste the sequence in the "Enter Query Sequence" imput box and choose the "Nucleotide collection (nr/nt)" datasbase to query. 

Step3: Run the algorithm to identify "Highly similar sequences (megablast)".

**Discussion: Look at the results. What can we say about the sequence? What type of organism is it?**

**To go further: What is a phage? What kind of organism does it infects?**

## 2- What biological functions does this sequence encode ?

Next we want to get an idea of what information we can retrieve on the functions carried in this DNA sequence. 

**Discussion: What is a gene? How can we identify a gene in a DNA sequence?**

### 2.A Unsing GeneMark to identify ORFs

Step1: To retrieve the ORFs from our sequence, we'll use the heuristic [GeneMark](http://opal.biology.gatech.edu/GeneMark/heuristic_gmhmmp.cgi).

Step2: Paste the sequence in the input box. Select "Protein sequence" as output. Leave the other parameters as default. 

Step3: Run GeneMark and open the protein sequence output.

**Discussion: How many ORFs did the algoritm identify in our sequence?**

## 2.B Explore putative protein functions with BlastP

Step1: Go to [BlastP](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins) at NCBI.

Step2: Paste our protein sequences (sequence available in unknown_proteins.fasta file) in the input query box, and choose to run BlastP against the "Reference proteins (RefSeq_protein)" database. Leave the other parameters as default

Step3: Browse the results for each of the putative proteins.

**Discussion: What type of functions did we identify? Is it surprising?**

Step4: Let's look closely at the protein identified by GenMark called 'gene_4'. 

Step5: Select the Top 10 hits for this protein and look at the result as a "Distance tree".

**Discussion: In what type of organism do we find this protein? What do you think happened here?**

## 3. Where can we find this sequence in ecosystems ?

Step 1: Finally, we want to explore the distribution and abundance of this gene in ecosystems. Go to the [Ocean Gene Atlas](http://tara-oceans.mio.osupytheas.fr/ocean-gene-atlas/) webpage.

Step 2: Using the 'gene_4|GeneMark.hmm|341_aa|+|9466|10491' as input, query the "OM-RGCv1 - Tara Ocean Microbiome" database. 

Step3: Set the "Expect threshold" to 1e-200.

**Discussion: Where this gene can be found in the ocean? Rather at the surface or deep in the sea ? Does that make sense from what we know about this gene ?**





