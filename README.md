Calculate Runs of Homozigosity (in Mb) using genotype likelihoods.
This was done as part of my master thesis

Input: Beagle File, more information on how to get this data at http://www.popgen.dk/angsd/index.php/Genotype_Likelihoods.
Output: Length of each ROH detected and total sum of ROHs

To calculate ROH call function calcROH(beagle, PhomMin, PhomMax), where beagle = beagle file; PhomMin = minimum Likelihood to be considered within a ROH if the adjacent sites have Likelihood > PhomMax; PhomMax = The likelihood requires for a site to be considered part of ROH.
eg: calcRoh(beagle = file.beagle, PhomMin = 0.6, PhomMax = 0.9)

To calculate the total sum of ROHs, call function calcSNroh(roh, minSize), where roh = list of roh length (output from calcROH) and minSize = the minimum size in Mb to consider when calculating the total sum of ROHs.
eg: calcSROH(calcRoh(file.beagle), minSize = 3)
![alt text](https://github.com/FCoroado/genolikeRoh/blob/master/genocalc.JPG?raw=true)

Author: Francisco Coroado Santos
