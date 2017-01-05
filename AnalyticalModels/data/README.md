The coefficients calculated bin by bin are in "coefficientsByBin_klkt.txt"

The bins are defined in the gen-level variables as:

binGenMHH = [250.,270.,300.,330.,360.,390., 420.,450.,500.,550.,600.,700.,800.,1000.] 
binGenCostS  = [ -1., -0.55,0.55,1.  ] 

   ==> The lines are respectivelly: nbins GenMhh GenCostStar NenventsSM NenventsSumV1 A1 A3 A7 errorA1 errorA3 error A7. 

   Where: NenventsSM and NenventsSumV1 is the number of events (in units of 10k events) in that bin,  
   as calculated respectivelly by simulations of  3M events for SM and the nevents calculated from the 12 benchmarks defined in 1507.02245v4 (JHEP version) each one with 100k events
   The numbering of the coefficients follow  arXiv:1608.06578

==================================================================================================================
In "Hist2DSum_V0_SM_box.root" you find the distribution for both 
     - the SM point (simulated from 100k events) - H0binX  
     - sum of the events simulated to Mriond 2016 (SM + box + the V0 benchmarks from 2-13 (simulated from ~300 events/each) - H1binX  

where X = 1 means the binning ( 90,0.,1800.,10,-1,1. ) and X = 1 means the same binning used for the analytical formula (see above)

==================================================================================================================
In "Distros_5p_SM3M_sumBenchJHEP_13TeV.root" you find the distribution for both 
     - the SM point (simulated from 3M events) - H0binX  
     - sum of the events from simulation of the 12 benchmarks defined in 1507.02245v4 (JHEP version) - H1binX 

===================================================================================================================
In "Distros_5p_500000ev_12sam_13TeV_JHEP_500K.root" you find the 2D histograms for the 12 benchmarks defined in 1507.02245v4 (JHEP version) ( 90,0.,1800.,10,-1,1. )

