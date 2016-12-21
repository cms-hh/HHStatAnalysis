In this repo we have the tools for constructing shapes in the kt X kl plane analitically.

We had calculated analytical formulas for signal efficiencies in bins. 
The bins are defined as:

binsx = [250.,270.,300.,330.,360.,390., 420.,450.,500.,550.,600.,700.,800.,1000.]
binsy  = [ -1., -0.55,0.55,1.  ]

The coefficients calculated bin by bin are in "coefficientsByBin_klkt.txt"

In the  "function_klkt.py" we have a template how to use the above file to calculate event-by-event weights. 

    ==> from lines 55-96 we read the coefficients
    ==> 216-223 we calculate the weights event by event
    ==> from lines 201-212 we read events from some files generated to test the formula

Of course you will need to to adapt the script and input events to your analysis workflow.
Note: 
- in root the bin counting in histograms start from 1 (and not 0). 
- The counting of benchmarks start from 1  

============================================================================================
==> The events for this template are in txt format (lines 105 and 106), and that is why it takes long to read and store
============================================================================================
By now you can test the formula works, the input events are assumed as the 12 benchmarks of the clustering JHEP version:

# It tests with simulation to: 
#  kl	kt	c2	cg	c2g		
#  1.0	1.0	0	0	0	: python function_klkt.py  --LHC 13 --kl 1 --kt 1 --c2 0 --cg 0 --c2g 0 --v1 False
# -10.	0.5	0	0	0	: python function_klkt.py  --LHC 13 --kl -10 --kt 0.5 --c2 0 --cg 0 --c2g 0  --v1 False
#  0.0001	2.25	0	0	0 : python function_klkt.py  --LHC 13 --kl 0.0001 --kt 2.25 --c2 0 --cg 0 --c2g 0  --v1 False

If you ask for a point that is not one of those will only draw to you the shape calculated by the reweighting, 
If you ask for one of those points it will superimpose it with an actual MC simulation

==> We stil need to adapt the denominator of the weight calculation to the v1 version (that we have fullsim)

===========================================================================================
The file Distros_5p_SM3M_sumBenchJHEP_13TeV.root is not yet explicitly used
 
