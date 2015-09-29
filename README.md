# MetaCRAM
Lossless compression tool for metagenomic reads

##Requirements and Dependencies

MetaCram is suitable for all unix-like operating systems with perl installation with perl modules 1) `File::Slurp` 2) `Array::Utils` 3) `Capture::Tiny`.

MetaCram uses the following software packages:
* Kraken - available from https://ccb.jhu.edu/software/kraken/
  - download minikraken database
* Bowtie2 - available from http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
* BLAST - available from http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download 
  - download nt database.
* IDBA - available from http://i.cs.hku.hk/~alse/hkubrg/projects/idba/
* CRAM - 3 jar files packaged together with the script 
  - ExtendedGolombCRAM.jar
  - GolombCRAM.jar 
  - HuffmanCRAM.jar
* SAMtools - available from samtools.sourceforge.net
* MFCompress - available from http://bioinformatics.ua.pt/software/mfcompress/

You will need to download and install these packages before running MetaCRAM.


MetaCram also makes use of the bacterial genome library from NCBI, which can be downloaded from 
ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/all.fna.tar.gz. 
After decompressing the file, we need to combine multiple strands into one folder. 

Alternatively, we have a genus database in the specific format needed by MetaCRAM, called "genus_db_final.tgz",  which can be downloaded here: https://uofi.box.com/s/sy92jv098eb0e4s277vuiod72m2jjf7c



##Commands to run MetaCRAM

**Compression**

`perl MetaCram.pl --compress --output <output directory> --paired <path to reads> --<exGolomb, huffman, golomb>`

Example:

`[shared3]$ perl MetaCram.pl --compress --output /shared3/MetaCRAM_SRR359032_Huffman --paired /shared3/SRR359032_1.fasta /shared3/ SRR359032_2.fasta --huffman &>MetaCramLOG_SRR359032_Huffman.txt`

**Decompression**

`perl MetaDeCram.pl --input <path to folder containing the Round1 and Round2 folders>`

Example:

`[shared3]$ perl MetaDeCram.pl --input /shared3/MetaCRAM_processedSRR359032_Huffman/MetaCRAM &>decompressorLogSRR359032_Huffman.txt`

(*--paired is optional)

(*<> indicates a choice)

(* to log, append “&> <log file>”)

##Options
For more options, run `perl MetaCram.pl --help`.

##Contact
Minji Kim (mkim158@illinois.edu) and Xiejia Zhang (xzhan121@illinois.edu)
