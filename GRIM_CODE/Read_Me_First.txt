The codes show how to evaluate the relative importance of random variables by using GRIM. The codes are based on the paper "Kim T and Song J 2018" published at Reliability Engineering and System Safety. The paper’s full reference information is as follows:

Kim, T., and J. Song (2018). Generalized reliability importance measure by Gaussian mixture.Reliability Engineering and System Safety. Vol. 173, 105-115.
DOI: https://doi.org/10.1016/j.ress.2018.01.005

To estimate GRIM, cross-entropy based adaptive importance sampling (CE-AIS-GM) is employed. For this purpose, this code employs the developed CE-AIS-GM Matlab code and a tool box of Finite Element Reliability Using Matlab (FERUM), directly. For more information, please visit http://systemreliability.wordpress.com and http://www.ce.berkeley.edu/projects/ferum/.

How to use the GRIM Matlab code
1. Run a CEAISGM input file (the CEAISGM input file is similar to FERUM's, but not the same)
2. Run CEAISGM on the command window of Matlab (type 'ceaisgm' on the command window)
3. The analysis results are available in the Matlab workspace (output and sensitivity).

Remark:
- The performance function is defined in input file (e.g. 'CEAISGM_inputfile_component.m'), one could modify this function to solve other problems of interest. Please note that the input file is DIFFERENT from the FERUM input file.

Acknowledgement:
The authors of GRIM would like to thank Ryan H.Y. Wong and Prof. Terje Haukaas who is the code developer of CE-AIS-GM and FERUM, respectively.