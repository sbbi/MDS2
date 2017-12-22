MDS2 is a tool for sequence motif identification on the short nucleotide sequences. MDS2 integrates the graph algorithms and large-scale sampling techniques to detect the sequence motif on the short sequences. It balances the specificity and sensitivity of motif and utilizes the statistical evaluation to minimize the risk of false identification. It has been proven that MDS2 is capable of predicting the motifs as short as 4bp-long, which are usually overlooked by other motif finding tools, such as MEME and DREME.
Online version can be found here: sbbi-panda.unl.edu/MDS2

Requirement:
1. MDS2 can only run on Linux system.
2. The following python packages need to be installed:
1. bitarray
2. graphviz
3. joblib
4. scipy
5. pickle
6. numpy
7. meme
8. gs-921-linux-x86_64
3. Background of RNA and miRNA data is required. They can be found here:

How to run:
1. Make links of uniprobe2meme and tomtom using cmd: 
    "ln -s [MEME_folder]/bin/uniprobe2meme ."
    "ln -s [MEME_folder]/bin/tomtom ."
2. Copy gs-921-linux-x86_64 to MDS2 folder.
3. MiRNA background data need to be renamed as "allHumanMiRNA(HSA).fa" and put in data folder.
4. RNA background data need to be renamed as "human_CDS.fasta" and put in data folder.
5. type in cmd "python GraphBasedMotifFinding.py [inputFolder] [maxMotifLength] [alphabet] [samplingFrequence] [enableMiRNASamplingOrNot] [useNegtiveInput] [clusteringMethod] [covPrecThres]"
inputFolder: The path of the folder which includes input fasta data. The fasta data need to be renamed as "userInput.fa".
maxMotifLength: The maximum of motif length. The default value is 6.
alphabet: The alphabet of the input. The default value is ACGU.
samplingFrequency: Frequency of sampling in background. The default value is 100000.
enableMiRNASamplingOrNot: Whether or not do sampling on miRNA background. The default value is False.
useNegtiveInput: Whether or not consider negtive data. The default value is False.
clusteringMethod: The method for get motif from different clusters. The default value is 2.
covPercThres: The threshold of coverage percentage considered as valid kmer. The default value is 0.1.

Run a test:
A test case is included in folder "output/testU". 
Use the following cmd to run a test:
python GraphBasedMotifFinding.py output/testU 6 ACGU 100 True False 2 0.1

Analysis of the result:
The final motif reported can be found in [inputFolder]/finalRlt.txt.
Logo files can be found in [inputFolder]/png/*