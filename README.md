# MUNIn
MUNIn (Multiple tissue UNifying long-range chromatin Interaction detector): a statistical framework for identifying long-range chromatin interactions from multiple tissues

MUNIn adopts a hierarchical hidden Markov random field (H-HMRF) model for identifying long-range chromatin interactions from multiple tissues, which is an extension of our previous HMRF peak caller HMRFBayesHiC (Xu et al., Bioinformatics, 2016). Comparing to the existing HiC peak calling methods, MUNin simultaneously account for spatial dependency within the same tissue, as well as dependency among different tissues. Specially, in MUNIn, the status of each pair of interacting chromatin loci (peak or non-peak) depends not only on the status of loci pairs in its neighborhood region, but also on the status of the same loci pair in other tissues.

MUNIn is maintained by Yuchen Yang [yyuchen@email.unc.edu]

## News and Updates
Sep 20, 2019
* Version 1.0 released
  + First offical release

## Installation
After downloading the MUNIn_1.0.tar.gz into a chosen local folder "local_path",
1. Unzip the file MUNIn_1.0.tar.gz, you will get a C++ executable program MUNIn, a folder Example including toy data, and a folder MUNIn_outputs including output results of MUNIn.
2. Copy the C++ executable program MUNIn to any folder and it is ready to use.

## MUNIn Examples
In this example, we use the TAD in chromosome 1 from 6595000 bp to 7965000 bp (denoted as "TAD_6595000_7965000") from two tissues, the cortical and subcortical plate (CP) and the germinal zone (GZ), at 10 KB resolution (Won et al. Nature, 2016). 

1. For each tissue, we start from HiC contact matrix, and calculate expected frequency using a modified version of Fit-Hi-C, which can be downloaded from here. The command interface of our utility software is exactly the same as Fit-Hi-C. Please refer to Fit-Hi-C for more details at https://noble.gs.washington.edu/proj/fit-hi-c/. 

2. We first perform peak calling in each tissue. To conduct peak calling, users need to prepare HiC data file for HiC_HMRF_Bayes_Files to load, which is a text file with 5 columns, separated by the table delimiter, respectively as tissue index, middle point of fragment 1, middle point of fragment 2, observed frequency and expected frequency. 

The 8 required command parameters are:

-I, HiC input data file, which is a text file, with 5 columns respectively as tissue index, middle point of fragment 1, middle point of fragment 2, observed frequency and expected frequency. The example file for CP is CP_1_6595000_7965000.txt.

-NP, size of HiC contact matrix.

-Tune, .

-NG, number of Gibbs sample.

-Bininitial, the middle point of the first fragment 1.

-Binsize, fragment length.

-SEED, seed of the random number generator. Setting the seed to a fixed value can make the results reproducible.

-O, output folder, which contains the output files of inferred peak status and parameters in the HMRF peak calling model. The example file is CP_output.

To run HMRF tissue by tissue, use 
./HMRF -I Example/CP_1_6595000_7965000.txt -NP 138 -Tune 100 -NG 10000 -Bininitial 6595000 -Binsize 10000 -SEED 123 -O Example/CP_output/ 

3. With the peak calling results from each tissue, we lable the tissues with different indices, i.e. 0, 1, 2..., and concatenate the long format output files together as the input file for MUNIn, which contains 6 columns respectively as tissue index, middle point of fragment 1, middle point of fragment 2, observed frequency, expected frequency and peak status. Users also need to prepare four files respectively containing the estimated parameters of theta, phi, gamma and psi of each tissue. 

The 8 required command parameters are
-I, input data file for MUNIn, which is a text file with 6 columns respectively as tissue index, middle point of fragment 1, middle point of fragment 2, observed frequency, expected frequency and peak status. The example file is CP_GZ_Record_long_format.txt.

-NP, size of HiC contact matrix.

-NT, number of tissues.

-NG, number of Gibbs sample.

-Bininitial, the middle point of the first fragment 1.

-Binsize, fragment length.

-Theta, theta input file, which is text file of one column listing the estimated theta from each tissue.

-Phi, phi input file, which is text file of one column listing the estimated phi from each tissue.

-Gamma, gamma input file, which is text file of one column listing the estimated gamma from each tissue.

-Psi, psi input file, which is text file of one column listing the estimated psi from each tissue.

-Alpha, tissue dependency input file. When there are two tissues, it contains 5 columns respectively as order index, peak status in tissue 1, peak statues in tissue 2, heterogeneity of peak status in the two tissues (0, shared background; 1, tissue-specific peak; 2, shared peak) and proportion of each status in all the fragment pairs. The example file is alpha_CP_GZ_1_6595000_7965000.txt.

-SEED, seed of the random number generator. Setting the seed to a fixed value can make the results reproducible.

-O, output folder, which contains the output files of inferred peak status and parameters in the HMRF peak calling model. The example file is MUNIn_outputs.

To run MUNIn, use 
./MUNIn -I Example/CP_GZ_Record_long_format.txt -NP 138 -NT 2 -NG 10000 -Bininitial 6595000 -Binsize 10000 -Theta Example/theta.txt -Phi Example/phi.txt -Gamma Example/gamma.txt -Psi Example/psi.txt -Alpha Example/alpha_CP_GZ_1_6595000_7965000.txt -SEED 1 -O MUNIn_output/ 




