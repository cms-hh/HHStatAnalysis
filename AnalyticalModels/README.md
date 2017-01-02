In this repo we have the tools for constructing shapes in the kt X kl plane analitically (this framework is extensible to c2, cg, c2g and more when necessary).

===========================================================================================

In the file "python nonResonant_test.py" there are functions to apply analytical formulas for signal efficiencies in bins. 

The bins are defined in the gen-level variables as:

binGenMHH = [250.,270.,300.,330.,360.,390., 420.,450.,500.,550.,600.,700.,800.,1000.] 
binGenCostS  = [ -1., -0.55,0.55,1.  ] 

The coefficients calculated bin by bin are in "coefficientsByBin_klkt.txt"

   ==> The lines are respectivelly: nbins GenMhh GenCostStar NenventsSM NenventsSumV1 A1 A3 A7 errorA1 errorA3 error A7. 

   Where: NenventsSM and NenventsSumV1 is the number of events (in units of 10k events) in that bin,  
   as calculated respectivelly by simulations of  3M events for SM and the nevents calculated from the 12 benchmarks defined in 1507.02245v4 (JHEP version) each one with 100k events
   The numbering of the coefficients follow  arXiv:1608.06578

============================================================================================

In the file "test/nonResonant_test.py" we have a template how to use the above file to calculate event-by-event weights. 

The functions main() and plotting() are merelly templates on how to apply the above mentionend functions
We calculate the weights event by event with the tree bellow lines (in the main() function):

By now you can test the formula works, the input events are assumed as the 12 benchmarks of the clustering JHEP version:

It tests with simulation to the following points: 

              kl	kt			
             1.0	1.0	: python nonResonant_test.py --LHC 13 --kl 1 --kt 1 --v1 0
             -10.	0.5	: python nonResonant_test.py --LHC 13 --kl -10 --kt 0.5   --v1 0
            0.0001	2.25	: python nonResonant_test.py --LHC 13 --kl 0.0001 --kt 2.25   --v1 0
            2.5		1.0	: python nonResonant_test.py --LHC 13 --kl 0.0001 --kt 2.25   --v1 0

If you ask for a point that is not one of those will only draw to you the shape calculated by the reweighting, 
If you ask for one of those points it will superimpose it with an actual MC simulation

==> We stil need to adapt the denominator of the weight calculation to the v1 version (matching to use the events we have from fullsim)

===========================================================================================
 
==> The events for this template are in txt format in AFS
