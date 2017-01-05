#! /usr/bin/env python
# Analytical reweighting implementation for H->4b
# This file is part of https://github.com/cms-hh/HHStatAnalysis.
# python nonResonant_test_JHEP.py  --kl 1 --kt 1
# compiling
from optparse import OptionParser
import ROOT
import numpy as np
from HHStatAnalysis.AnalyticalModels.NonResonantModel import NonResonantModel

parser = OptionParser()
parser.add_option("--kl", type="float", dest="kll", help="Multiplicative factor in the H trilinear wrt to SM")
parser.add_option("--kt", type="float", dest="ktt", help="Multiplicative factor in the H top Yukawa wrt to SM")

parser.add_option("--c2", type="float", dest="c22", help="ttHH with triangle loop", default=0)
parser.add_option("--cg", type="float", dest="cgg", help="HGG contact", default=0)
parser.add_option("--c2g", type="float", dest="c2gg", help="HHGG contact", default=0)

(options, args) = parser.parse_args()
print " "
kl = options.kll
kt = options.ktt
c2 = options.c22
cg = options.cgg
c2g = options.c2gg

print "Weights calculated from the 12 benchmarks defined in 1507.02245v4 (JHEP version) each one with 100k events "
if c2 != 0 or cg != 0 or c2g != 0 :  print "The analytical function is not yet implemented"

###########################################################
# read events and apply weight
###########################################################
def main():
  model = NonResonantModel()
  # obtaining BSM/SM coeficients
  dumb = model.ReadCoefficients("../data/coefficientsByBin_klkt.txt",model.effSM,model.effSum,model.MHH,model.COSTS,model.A1,model.A3,model.A7) 
  # We sum SM + box + the benchmarks from 2-13 
  # read the 2D histo referent to the sum of events
  fileHH=ROOT.TFile("../data/Distros_5p_SM3M_sumBenchJHEP_13TeV.root")
  sumJHEPAnalyticalBin = fileHH.Get("H1bin2")
  SMAnalyticalBin = fileHH.Get("H0bin2")
  xaxis = sumJHEPAnalyticalBin.GetXaxis()
  yaxis = sumJHEPAnalyticalBin.GetYaxis()
  print "Sum hist ",sumJHEPAnalyticalBin.GetNbinsX(),sumJHEPAnalyticalBin.GetNbinsY(),sumJHEPAnalyticalBin.Integral(),sumJHEPAnalyticalBin.GetXaxis().GetBinLowEdge(1),sumJHEPAnalyticalBin.GetXaxis().GetBinUpEdge(xaxis.GetNbins())
  #print SMAnalyticalBin.GetBinContent(4,4)
  # now loop over events, calculate weights using the coeffitients and  plot histograms
  # events to reweights, in text format (for testing only)
  pathBenchEvents="/afs/cern.ch/work/a/acarvalh/public/toAnamika/GF_HH_BSM/" # events to reweight
  # declare the histograms 
  CalcMhh = np.zeros((1200000))
  CalcCost = np.zeros((1200000))
  CalcWeight = np.zeros((1200000))
  ##########################################
  # initialize tables of coefficients by bins
  # calculate mhh and cost* and find the bin
  # initialize events reading 
  countline=0
  countevent=0
  counteventSM=0
  # read events as text files for events to test 
  # particuliarity of the text file with events = each 2 lines are one event there
  # save the information of the two Higgses
  Px = np.zeros((2)) 
  Py = np.zeros((2)) 
  Pz = np.zeros((2)) 
  En = np.zeros((2))  
  for sam in  range(1,13): # read events from the list of 12 benchmarks
       filne = pathBenchEvents+"GF_HH_"+str(sam)+".lhe.decayed"    # 0 is SM = here it does not enter in the events to be reweighted in this version
       f = open(filne, 'r+')
       lines = f.readlines() # get all lines as a list (array)
       countline = 0 # particuliarity of the text file with events = each 2 lines are one event there
       for line in  lines:
          model.ReadLine(line, countline,Px,Py,Pz,En)
          #print countline
          countline+=1
          mhhcost= [0,0] # to store [mhh , cost] of that event
          if countline==2 : # if read 2 lines 
            model.CalculateMhhCost(mhhcost,countline,Px,Py,Pz,En) # ==> adapt to your input 
            bmhh = sumJHEPAnalyticalBin.GetXaxis().FindBin(mhhcost[0])
            bcost = sumJHEPAnalyticalBin.GetYaxis().FindBin(mhhcost[1])
            effSumV0 = sumJHEPAnalyticalBin.GetBinContent(bmhh,bcost)  # quantity of simulated events in that bin (without cuts)
            #weight = model.getScaleFactor(mhhcost,kl, kt,0,model.effSM,model.effSum,model.MHH,model.COSTS,model.A1,model.A3,model.A7,0)  
            #print mhhcost[1],bcost
            #print effSumV0
            weight = model.getScaleFactor(mhhcost,kl, kt,model.effSM,model.MHH,model.COSTS,model.A1,model.A3,model.A7,effSumV0)  
            countline=0
            #############################################
            # fill histograms
            #############################################
            if weight > 0: 
               #print countevent
               CalcMhh[countevent] = float(mhhcost[0]) 
               CalcCost[countevent] = float(mhhcost[1]) 
               CalcWeight[countevent] = weight 
               countevent+=1
       f.close()
  print "plotted hostogram reweighted from ",countevent," events, ", float(100*(1200000-countevent)/1200000)," % of the events was lost in empty bins in SM simulation"
  ############################################################################################################################
  # Draw test histos
  ###############################################
  drawtest =0 
  nevtest=50000
  if kl == 1 and kt == 1 and c2 ==0 and cg == 0 and c2g ==0 : 
     filne = pathBenchEvents+"GF_HH_0.lhe.decayed"    # 0 is SM
     nevtest = 100000
     drawtest = 1 
  # BSM events
  pathBSMtest="/afs/cern.ch/work/a/acarvalh/public/toAnamika/GF_HH_toRecursive/" # events of file to superimpose a test
  # see the translation of coefficients for this last on: If you make this script smarter (to only read files we ask to test) you can implement more
  # https://github.com/acarvalh/generateHH/blob/master/fit_GF_HH_lhe/tableToFitA3andA7.txt
  if kl == -10 and kt == 0.5 and c2 ==0 and cg == 0 and c2g ==0 :
     drawtest =1
     filne = pathBSMtest+"GF_HH_42.lhe.decayed"
  if kl == 0.0001 and kt == 2.25 and c2 ==0 and cg == 0 and c2g ==0  :
     drawtest =1
     filne = pathBSMtest+"GF_HH_9.lhe.decayed"
  if kl == 2.5 and kt == 1.0 and c2 ==0 and cg == 0 and c2g ==0  :
     drawtest =1
     filne = pathBSMtest+"GF_HH_60.lhe.decayed"
  ############################################################################################################################
  CalcMhhTest = np.zeros((nevtest))
  CalcCostTest = np.zeros((nevtest))
  if drawtest ==1 :
     print "draw plain histogram to test"
     model.LoadTestEvents(CalcMhhTest,CalcCostTest,filne)  
  model.plotting(kl,kt,CalcMhh,CalcCost,CalcWeight,CalcMhhTest,CalcCostTest,drawtest)
###############################################################################################################################
if __name__ == "__main__":  
   main()


#print len(MHH),A1[0][0]

#options.kll, options.ktt)


# obtaining BSM/SM scale factors

#canvas = ROOT.TCanvas("")
#scaleFactors.Draw('colz')
#canvas.SaveAs("{}.pdf".format(scaleFactors.GetName()))
