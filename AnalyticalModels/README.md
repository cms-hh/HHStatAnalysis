In this repo we have the tools for constructing shapes in the kt X kl plane analitically (this framework is extensible to c2, cg, c2g and more when necessary).

cmsrel CMSSW_7_4_7
cd CMSSW_7_4_7/src
cmsenv
tar xzvf HHStatAnalysis.tar.gz
scram b
cd HHStatAnalysis/AnalyticalModels/test
python nonResonant_bbbb.py --kl 1 --kt 1


The input elements are described inside the folder "data"

===========================================================================================

In the file "test/nonResonant_test_v(0/JHEP).py" we have a template how to use the above file to calculate event-by-event weights. 
In suma, To calculate the event weight one need the following three lines

      mhhcost= [tree.Genmhh(),tree.GenHHCost()] # to store [mhh , cost] of that event
      effSum = sumHAnalyticalBin.GetBinContent(bmhh,bcost)  # quantity of simulated events in that bin (without cuts)
      weight = model.getScaleFactor(mhhcost,kl, kt,v1,model.effSM,effSum,model.MHH,model.COSTS,model.A1,model.A3,model.A7) 

The sum of weights without cuts is automatically equal to 1. (with ~ 5% precision), so the sum of weights after cuts is already the signal efficiency. 

=============================================================================================

To test this template the events are in a public space AFS

  ==> The events for the 12 benchmarks defined in 1507.02245v4 (JHEP version) each one with 100k event are in txt format (a reader is implemented in the python class)
      ==> We sum the benchmarks from 1-12 

  ==> The events for V0 (the same of the fullsim version of Moriond 2016) are in root format. Those are common to all CMS final states
      ==> We sum SM + box + the V0 benchmarks from 2-13

The functions main() and plotting() are merelly templates on how to apply the above mentionend functions
We calculate the weights event by event with the tree bellow lines (in the main() function):

============================================================================================

It tests with simulation (events stored in txt format in AFS) against the reweigthed distribution to the following points: 

              kl	kt			
             1.0	1.0	: python nonResonant_test.py --LHC 13 --kl 1 --kt 1 --v1 0
             -10.	0.5	: python nonResonant_test.py --LHC 13 --kl -10 --kt 0.5   --v1 0
            0.0001	2.25	: python nonResonant_test.py --LHC 13 --kl 0.0001 --kt 2.25   --v1 0
            2.5		1.0	: python nonResonant_test.py --LHC 13 --kl 0.0001 --kt 2.25   --v1 0

If you ask for a point that is not one of those will only draw to you the shape calculated by the reweighting, 
If you ask for one of those points it will superimpose it with an actual MC simulation

===========================================================================================


