# MUNIn
MUNIn (Multiple tissue UNifying long-range chromatin Interaction detector): a statistical framework for identifying long-range chromatin interactions from multiple tissues/samples.

MUNIn adopts a hierarchical hidden Markov random field (H-HMRF) model for identifying long-range chromatin interactions from multiple samples, which is an extension of our previous HMRF peak caller HMRFBayesHiC (Xu et al., Bioinformatics, 2016). Comparing to the existing HiC peak calling methods, MUNin simultaneously account for spatial dependency within the same sample, as well as dependency among different samples. Specially, in MUNIn, the status of each pair of interacting chromatin loci (peak or non-peak) depends not only on the status of loci pairs in its neighborhood region, but also on the status of the same loci pair in other samples.

MUNIn is maintained by Yuchen Yang [yyuchen@email.unc.edu], Ming Hu [hum@ccf.org] and Yun Li [yun_li@med.unc.edu].

## News and Updates
Sep 20, 2019
* Version 1.0 released
  + First offical release

## Installation
After downloading the MUNIn_1.0.tar.gz into a chosen local folder "local_path",
1. Unzip the file MUNIn_1.0.tar.gz, you will get a C++ executable program MUNIn, a folder Example including toy data, and a folder MUNIn_outputs including output results of MUNIn.
2. Copy the C++ executable program MUNIn to any folder and it is ready to use.

## MUNIn Examples
In this example, we use the TAD in chromosome 1 from 50875000 bp to 51725000 bp (denoted as "TAD_50875000_51725000") from two cell lines, GM12878 and IMR90, at 10 KB resolution (Rao et al. Nature, 2016). 

### Calculate expected contact frequency using Fit-Hi-C
For each sample, we start from HiC contact matrix, and calculate expected contact frequency using a modified version of Fit-Hi-C, which can be downloaded from here. The command interface of our utility software is exactly the same as Fit-Hi-C. Please refer to Fit-Hi-C for more details at https://noble.gs.washington.edu/proj/fit-hi-c/. 

### Call peaks for each sample separately
We first perform peak calling in each sample separately using ***H-HMRF*** method. To conduct peak calling, users need to prepare HiC data file for HiC_HMRF_Bayes_Files to load, which is a text file with 5 columns, separated by the table delimiter, respectively as sample index, middle point of fragment 1, middle point of fragment 2, observed frequency and expected frequency and p-value estimated by Fit-Hi-C. For example, the first several lines of GM_1_50875000_51725000.txt are 

50875000        50885000        820     511.407636      2.803035e-36 <br>
50875000        50895000        383     264.041244      4.051192e-12 <br>
50875000        50905000        272     184.994981      1.315629e-09 <br>
50875000        50915000        186     149.215834      2.027127e-03 <br>
50875000        50925000        124     126.085980      5.855781e-01 <br>
... <br>

***The 8 required command parameters*** by H-HMRF include:

-I, HiC input data file, which is a text file, with 5 columns respectively as sample index, middle point of fragment 1, middle point of fragment 2, observed frequency and expected frequency. The example file for GM12878 is GM_1_50875000_51725000.txt.

-NP, size of HiC contact matrix.

-Tune, .

-NG, number of Gibbs sample.

-Bininitial, the middle point of the first fragment 1.

-Binsize, fragment length.

-SEED, seed of the random number generator. Setting the seed to a fixed value can make the results reproducible.

-O, output folder, which contains the output files of inferred peak status and parameters in the HMRF peak calling model. The example file is GM_output.

The command line for executing ***H-HMRF method*** in each sample is <br>
./HMRF -I Example/GM_1_50875000_51725000.txt -NP 138 -Tune 100 -NG 10000 -Bininitial 50875000 -Binsize 10000 -SEED 123 -O Example/GM_output/ 

### Call peaks across samples using MUNIn
With the peak calling results from each sample, we lable the sample with different indices, i.e. 0, 1, 2..., and concatenate the long format output files together as the input file for MUNIn, which contains 6 columns respectively as sample index, middle point of fragment 1, middle point of fragment 2, observed frequency, expected frequency and peak status. For example, the first several lines of GM_IMR90_Record_long_format.txt are

// sample_index	frag1	frag2	Oij Eij	peak_status <br>
0 50875000  50885000 820 511.407636  -1 <br>
0 50875000  50895000 383 264.041244  -1 <br>
0 50875000  50905000 272 184.994981  -1 <br>
0 50875000  50915000 186 149.215834  -1 <br>
0 50875000  50925000 124 126.085980  -1 <br>
... <br>
1 50875000  50885000  248 156.12634000  -1 <br>
1 50875000  50895000  89  77.77461700 -1 <br>
1 50875000  50905000  71  54.13943900 -1 <br>
1 50875000  50915000  58  43.29327000 -1 <br>
1 50875000  50925000  26  35.84396800 -1 <br>
...

MUNIn also requires an alpha file, which show the dependency level between different samples (e.g. tissues or cells lines). Alpha file is a text file with 5 columns, when there are two samples, respectively as order index, peak status in sample 1, peak statues in sample 2, heterogeneity of peak status in the two samples (0, shared background; 1, sample-specific peak; 2, shared peak) and proportion of each status in all the fragment pairs. Here is an example for alpha_GM_IMR90_1_50875000_51725000.txt

// order_index	peak_status_sample1	peak_status_sample2	proportion <br>
0	0	0	0 0.51573187 <br>
1	1	0	1	0.18276334 <br>
2	0	1	1	0.02927497 <br>
3	1	1	2	0.27222982 <br>

Users also need to files for each of the four parameters, theta, phi, gamma and psi, according to the estimation results from H-HMRF method. The file of each parameter is a text file of one column listing the estimated parameter, for example phi, from each sample. Here is an example for phi.txt

11.9797 <br>
9.1174 

***The 8 required command parameters*** by MUNIn include: 
-I, input data file for MUNIn, which is a text file with 6 columns respectively as sample index, middle point of fragment 1, middle point of fragment 2, observed frequency, expected frequency and peak status. The example file is GM_IMR90_Record_long_format.txt.

-NP, size of HiC contact matrix.

-NT, number of samples.

-NG, number of Gibbs sample.

-Bininitial, the middle point of the first fragment 1.

-Binsize, fragment length.

-Theta, theta input file, which is text file of one column listing the estimated theta from each sample.

-Phi, phi input file, which is text file of one column listing the estimated phi from each sample.

-Gamma, gamma input file, which is text file of one column listing the estimated gamma from each sample.

-Psi, psi input file, which is text file of one column listing the estimated psi from each sample.

-Alpha, sample dependency input file. When there are two samples, it contains 5 columns respectively as order index, peak status in sample 1, peak statues in sample 2, heterogeneity of peak status in the two samples (0, shared background; 1, sample-specific peak; 2, shared peak) and proportion of each status in all the fragment pairs. The example file is alpha_GM_IMR90_1_50875000_51725000.txt.

-SEED, seed of the random number generator. Setting the seed to a fixed value can make the results reproducible.

-O, output folder, which contains the output files of inferred peak status and parameters in the HMRF peak calling model. The example file is MUNIn_outputs.

The command line for ***MUNIn***, use 
./MUNIn -I Example/GM_IMR90_Record_long_format.txt -NP 86 -NT 2 -NG 10000 -Bininitial 50875000 -Binsize 10000 -Theta Example/theta.txt -Phi Example/phi.txt -Gamma Example/gamma.txt -Psi Example/psi.txt -Alpha Example/alpha_GM_IMR90_1_50875000_51725000.txt -SEED 1 -O MUNIn_output/

## Output formats of MUNIn
MUNIn outputs multiple files, the majority of which are files recoding peak status and parameters for each sample.
The output Hi-C peak recode file is a text file, with 5 columns respectively as middle point of fragment 1, middle point of fragment 2, observed frequency, expected frequency, peak status and original Fit-Hi-C pvalue in sample 1 and sample 2, respectively, and posterior probablitis for shared peaks, sample 1-specific peak and sample 2-specific peak. For example, the first ten lines of PP.txt are

// frag1	frag2	Oij_sample1	Eij_sample1	peak_status_sample1 original_Fit-Hi-C_pvalue_sample1  Oij_sample2	Eij_sample2	peak_status_sample2 original_Fit-Hi-C_pvalue_sample2 <br>
50875000  50885000  820 511.4076  1 0.000000e+00  248 156.1263  -1  0.000000e+00  0.32240863  0.39526601  0.12683203 <br>
50875000  50895000  383 264.0412  1 0.000000e+00  89  77.7746 -1  1.134804e-01  0.00105054  0.09128604  0.01032674 <br>
50875000  50905000  272 184.9950  -1  0.000000e+00  71  54.1394 -1  1.596326e-02  0.00008979  0.10087206  0.00079958 <br>
50875000  50915000  186 149.2158  -1  2.027130e-03  58  43.2933 -1  1.881117e-02  0.00000319  0.00319080  0.00099670 <br>
50875000  50925000  124 126.0860  -1  5.855781e-01  26  35.8440 -1  9.635358e-01  0.00000008  0.00083995  0.00009658 <br>
...

The output parameter recode file is a text file, with maximum likelihood value, estimated parameters of theta, phi, gamma and psi, and the number of Gibbs sample and seed used for MUNIn outputted. For example, the texts of the file Record_Para.txt are

Best LogLike = 1503260.6415 <br>
Tissue = 0 <br>
Best Theta = 0.7971 <br>
Best Phi = 9.8234 <br>
Best Gamma = -0.0185 <br>
Best Psi = 0.6000 <br>
Tissue = 1 <br>
Best Theta = 0.7693 <br>
Best Phi = 8.2968 <br>
Best Gamma = -0.0294 <br>
Best Psi = 0.5982 <br>
NumGibbs = 10000 <br>
SEED = 1

## Citation
Liu, W., Abnousi, A., Zhang, Q., Li, Y., Hu, M., Yang, Y. (2021+) MUNIn (Multiple cell-type UNifying long-range chromatin Interaction detector): a statistical framework for identifying long-range chromatin interactions from multiple cell types. *biorxiv*, [DOI: 10.1101/2020.11.12.380782](https://www.biorxiv.org/content/10.1101/2020.11.12.380782v1)


