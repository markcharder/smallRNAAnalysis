////////////////////////////////////////////////////////////////////////////////////// HHMMiR Classifier//// Written By: Sabah Kadri//// Programs://	1. ExtractHairpins.jar//	2. HHMMiR.jar//	3. fold_regions.pl//// Reference:// S. Kadri, V. Hinman, P.V. Benos, 
// "HHMMiR: Efficient de novo prediction of microRNAs using hierarchical hidden Markov models."// BMC Bioinformatics (Proc APBC), (2009) 10(Suppl 1):S35//// PDF: http://genomebiology.com/content/pdf/gb-2009-10-2-r18.pdf//
//
//
// VERSIONS
// ========
// Started: March 2008
// Version 1 (released): January 2009 
//// Version 1.2: April 2009 
// Comment: correct a bug (incorectly repeated base in the overlap)//
//
//
//// Copyright 2008 Sabah Kadri//// This file is part of HHMMiR.//// HHMMiR is free software.  It is distributed in the hope that// it will be useful, but WITHOUT ANY WARRANTY; without even the// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR// PURPOSE.//////////////////////////////////////////////////////////////////////////////////////HHMMiR Classifier

About HHMMiR: This program can be used to classify a given genomic hairpin as a microRNA precursor or a random hairpin (without a microRNA gene). For this, we require the secondary stucture of the RNA sequence. The input to the program is a file containing the query single-loop hairpins and the output is a comma-separated file that contains the hairpin name, its sequence, the folded sequence, the log likelihood from the positive model, the log likelihood from the negative model, the ratio of the log likelihood scores and finally, based on the threshold specified, TRUE(microRNA precursor) or FALSE(random hairpin).
Ref[1] has the details.

Pre-requisites:
You will need Java (at least 5.0) installed to run the jar file.You will also need a Perl compiler.

Follow these steps to make a prediction:
1. Make a secondary structure prediction using RNAfold (the perl file will do this)[2]
2. Parse and extract the single-loop hairpins from the above file as input to HHMMiR
3. Use HHMMiR to classify the hairpins

For example:

1. If you have a section of the genome, as in "seq_region.fa" or simply, a list of sequences in Fasta format. Run the perl script as follows:

perl fold_regions.pl <Complete path for RNAfold> <fastafile>

i.e.

perl fold_regions.pl /home/ViennaPackage/ViennaRNA-1.6.5/Progs/RNAfold seq_region.fa

This will take the input sequence and fold windows of 1kb with overlap of 150nts between consecutive windows using RNAfold[2]. The output file is "folded.fa"

NOTE 1 : This will generate two other files "ip.txt" and "rna.ps". You do not need these further so you can delete these.
NOTE 2 : RNAfold takes some time 

2. Now parse "folded.fa" and extract all single loop hairpins. Note you have to specify:

Minimum hairpin length 
Maximum Bulge length
Minimum base pairs in the hairpin for this purpose
If not specified, these are set by default to: 55,19,15 respectively.

java -jar ExtractHairpins.jar <inputFile> <outputFile> >> removed.txt


For example:

java -jar ExtractHairpins.jar folded.fa hairpins.fa >> removed.txt

OR

java -jar ExtractHairpins.jar folded.fa hairpins.fa 50 25 10 >> removed.txt

NOTE: "removed.txt" will contain the hairpins that did not meet the criteria and the reason for each removal.

3. Now use the single looped hairpins returned above, as input to HHMMiR.jar to make predictions.

java -jar HHMMiR.jar <positive_model> <negative_model> <input_file> <output_csv_file> <threshold>

where <positive_model> is the parameter file for the model trained on positive dataset
<negative_model> is the parameter file for the model trained on negative dataset
<input_file> is the hairpin file from above
<threshold> = optimal (0.71) for BaumWelch trained model
					  (0.98) for MLE trained model 
					  Refer to [1] for details

For example, for above,

java -jar HHMMiR.jar posModel.txt negModel.txt hairpins.fa output.csv 0.71

The output comma separated file has the following columns:
Name_of_hairpin,Hairpin_seq,Hairpin_fold,positive_likelihood,negative_likelihood,ratio_of_likelihoods,TRUE/FALSE

****************************
In summary (for the example_files provided to you,

perl fold_regions.pl /home/ViennaPackage/ViennaRNA-1.6.5/Progs/RNAfold seq_region.fa

java -jar ExtractHairpins.jar folded.fa hairpins.fa 50 25 10 >> removed.txt

java -jar HHMMiR.jar posModel.txt negModel.txt hairpins.fa output.csv 0.71

****************************
References:
1. Sabah Kadri, Veronica Hinman and Panayiotis V Benos: HHMMiR: efficient de novo prediction of microRNAs using hierarchical hidden Markov models; BMC Bioinformatics 2009, 10(Suppl 1):S35doi:10.1186/1471-2105-10-S1-S35
2. Hofacker IL:  Vienna RNA secondary structure server. Nucleic Acids Res 2003, 31(13):3429-3431.


*Please contact Sabah Kadri at sabah.kadri@gmail.com for any further questions or datasets*