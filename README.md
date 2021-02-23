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
We first perform peak calling in each sample separately using ***H-HMRF*** method. To conduct peak calling, users need to prepare HiC data file for HiC_HMRF_Bayes_Files to load, which is a text file with 5 columns, separated by the table delimiter, respectively as sample index, middle point of fragment 1, middle point of fragment 2, observed frequency and expected frequency and p-value estimated by Fit-Hi-C. 

```
zcat Fithic_output/GM12878_chr1_10kb.spline_pass2.significances.txt.gz | awk '$2>=50875000 && $2<=51725000 && $4>=50875000 && $4<=51725000 && $4<$2' | awk '{print "0\t"$4"\t"$2"\t"$5"\t"$8"\t"$6}' | sort -nk 2 >GM12878_1_50875000_51725000.txt
```

Here are the first several lines of GM12878_1_50875000_51725000.txt

```
50875000        50885000        820     511.407636      2.803035e-36
50875000        50895000        383     264.041244      4.051192e-12
50875000        50905000        272     184.994981      1.315629e-09
50875000        50915000        186     149.215834      2.027127e-03
50875000        50925000        124     126.085980      5.855781e-01
... 
```

***The 8 required command parameters*** by H-HMRF include:

-I, HiC input data file, which is a text file, with 5 columns respectively as sample index, middle point of fragment 1, middle point of fragment 2, observed frequency and expected frequency. The example file for GM12878 is GM_1_50875000_51725000.txt.

-NP, size of HiC contact matrix.

-Tune, .

-NG, number of Gibbs sample.

-Bininitial, the middle point of the first fragment 1.

-Binsize, fragment length.

-SEED, seed of the random number generator. Setting the seed to a fixed value can make the results reproducible.

-O, output folder, which contains the output files of inferred peak status and parameters in the HMRF peak calling model. The example file is GM_output.

The command lines for executing ***H-HMRF method*** in each of the two samples, GM12878 and IMR90, are <br>

```
mkdir GM12878_output
mkdir IMR90_output
./HMRF -I GM12878_1_50875000_51725000.txt -NP 138 -Tune 100 -NG 10000 -Bininitial 50875000 -Binsize 10000 -SEED 123 -O GM12878_output/
./HMRF -I IMR90_1_50875000_51725000.txt -NP 138 -Tune 100 -NG 10000 -Bininitial 50875000 -Binsize 10000 -SEED 123 -O IMR90_output/

mkdir HMRF_output
mv GM12878_output HMRF_output
mv IMR90_output HMRF_output
```

### Call peaks across samples using MUNIn
With the peak calling results from each sample, we lable the sample with different indices, i.e. 0, 1, 2..., and concatenate the long format output files together as the input file for MUNIn. For example,

```
awk '{print "0\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' GM12878_output/Record_long_format.txt >GM12878_output/GM12878_1_50875000_51725000_long_format.txt
awk '{print "0\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' IMR90_output/Record_long_format.txt >IMR90_output/IMR90_1_50875000_51725000_long_format.txt
cat GM12878_output/GM12878_1_50875000_51725000_long_format.txt IMR90_output/IMR90_1_50875000_51725000_long_format.txt >GM12878_IMR90_Record_long_format.txt
```

The input file contains 6 columns respectively as sample index, middle point of fragment 1, middle point of fragment 2, observed frequency, expected frequency and peak status. For example, the first several lines of GM12878_IMR90_Record_long_format.txt are

```
// sample_index	frag1	frag2	Oij Eij	peak_status
0 50875000  50885000 820 511.407636  -1
0 50875000  50895000 383 264.041244  -1
0 50875000  50905000 272 184.994981  -1
0 50875000  50915000 186 149.215834  -1
0 50875000  50925000 124 126.085980  -1
...
1 50875000  50885000  248 156.12634000  -1
1 50875000  50895000  89  77.77461700 -1
1 50875000  50905000  71  54.13943900 -1
1 50875000  50915000  58  43.29327000 -1
1 50875000  50925000  26  35.84396800 -1
...
```

MUNIn requires an alpha file, which contains the dependency level between different samples (e.g. tissues or cells lines). Sample dependency can be estimated using ***estAlpha*** function via the folloing command line. 

```
./estAlpha -I GM_IMR90_Record_long_format.txt -NP 138 -NT 2 -Bininitial $begin -Binsize 10000 -O ./
mv alpha.txt alpha_GM12878_IMR90_1_50875000_51725000.txt
```

Alpha file is a text file with 5 columns, when there are two samples, respectively as order index, peak status in sample 1, peak statues in sample 2, heterogeneity of peak status in the two samples (0, shared background; 1, sample-specific peak; 2, shared peak) and proportion of each status in all the fragment pairs. Here is an example for alpha_GM12878_IMR90_1_50875000_51725000.txt

```
// order_index	peak_status_sample1	peak_status_sample2	proportion
0	0	0	0 0.51573187
1	1	0	1	0.18276334
2	0	1	1	0.02927497
3	1	1	2	0.27222982
```

MUNIn also requires files for each of the four parameters, theta, phi, gamma and psi, according to the estimation results from H-HMRF method. For example, users can generate the file from the H-HMRF output via following commands:

```
phi1=`head -n 2 HMRF_outdir/GM12878/Record_Para.txt | tail -n 1 | cut -f 4 -d ' '`
phi2=`head -n 2 HMRF_outdir/IMR90/Record_Para.txt | tail -n 1 | cut -f 4 -d ' '`
echo $phi1 > phi.txt
echo $phi2 >> phi.txt
```

The file of each parameter is a text file of one column listing the estimated parameter, for example phi, from each sample. Here is an example for phi.txt

```
11.9797
9.1174 
```

***The 8 required command parameters*** by MUNIn include: 
-I, input data file for MUNIn, which is a text file with 6 columns respectively as sample index, middle point of fragment 1, middle point of fragment 2, observed frequency, expected frequency and peak status. The example file is GM12878_IMR90_Record_long_format.txt.

-NP, size of HiC contact matrix.

-NT, number of samples.

-NG, number of Gibbs sample.

-Bininitial, the middle point of the first fragment 1.

-Binsize, fragment length.

-Theta, theta input file, which is text file of one column listing the estimated theta from each sample.

-Phi, phi input file, which is text file of one column listing the estimated phi from each sample.

-Gamma, gamma input file, which is text file of one column listing the estimated gamma from each sample.

-Psi, psi input file, which is text file of one column listing the estimated psi from each sample.

-Alpha, sample dependency input file. When there are two samples, it contains 5 columns respectively as order index, peak status in sample 1, peak statues in sample 2, heterogeneity of peak status in the two samples (0, shared background; 1, sample-specific peak; 2, shared peak) and proportion of each status in all the fragment pairs. The example file is alpha_GM12878_IMR90_1_50875000_51725000.txt.

-SEED, seed of the random number generator. Setting the seed to a fixed value can make the results reproducible.

-O, output folder, which contains the output files of inferred peak status and parameters in the HMRF peak calling model. The example file is MUNIn_outputs.

The command line for ***MUNIn***, use 
```
./MUNIn -I GM_IMR90_Record_long_format.txt -NP 86 -NT 2 -NG 10000 -Bininitial 50875000 -Binsize 10000 -Theta theta.txt -Phi phi.txt -Gamma gamma.txt -Psi psi.txt -Alpha alpha_GM_IMR90_1_50875000_51725000.txt -SEED 1 -O MUNIn_output/
```

## Output formats of MUNIn
MUNIn outputs multiple files, the majority of which are files recoding peak status and parameters for each sample.
The output Hi-C peak recode file is a text file, with 5 columns respectively as middle point of fragment 1, middle point of fragment 2, observed frequency, expected frequency, peak status and original Fit-Hi-C pvalue in sample 1 and sample 2, respectively, and posterior probablitis for shared peaks, sample 1-specific peak and sample 2-specific peak. For example, the first ten lines of PP.txt are

```
// frag1	frag2	Oij_sample1	Eij_sample1	peak_status_sample1 original_Fit-Hi-C_pvalue_sample1  Oij_sample2	Eij_sample2	peak_status_sample2 original_Fit-Hi-C_pvalue_sample2
50875000  50885000  820 511.4076  1 0.000000e+00  248 156.1263  -1  0.000000e+00  0.32240863  0.39526601  0.12683203
50875000  50895000  383 264.0412  1 0.000000e+00  89  77.7746 -1  1.134804e-01  0.00105054  0.09128604  0.01032674
50875000  50905000  272 184.9950  -1  0.000000e+00  71  54.1394 -1  1.596326e-02  0.00008979  0.10087206  0.00079958
50875000  50915000  186 149.2158  -1  2.027130e-03  58  43.2933 -1  1.881117e-02  0.00000319  0.00319080  0.00099670
50875000  50925000  124 126.0860  -1  5.855781e-01  26  35.8440 -1  9.635358e-01  0.00000008  0.00083995  0.00009658
...
```

The output parameter recode file is a text file, with maximum likelihood value, estimated parameters of theta, phi, gamma and psi, and the number of Gibbs sample and seed used for MUNIn outputted. For example, the texts of the file Record_Para.txt are

```
Best LogLike = 1503260.6415
Tissue = 0
Best Theta = 0.7971
Best Phi = 9.8234
Best Gamma = -0.0185
Best Psi = 0.6000
Tissue = 1
Best Theta = 0.7693
Best Phi = 8.2968
Best Gamma = -0.0294
Best Psi = 0.5982
NumGibbs = 10000
SEED = 1
```

## Citation
Liu, W., Abnousi, A., Zhang, Q., Li, Y., Hu, M., Yang, Y. (2021+) MUNIn (Multiple cell-type UNifying long-range chromatin Interaction detector): a statistical framework for identifying long-range chromatin interactions from multiple cell types. *biorxiv*, [DOI: 10.1101/2020.11.12.380782](https://www.biorxiv.org/content/10.1101/2020.11.12.380782v1)


