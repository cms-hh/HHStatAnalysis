#!/usr/bin/env python
# RUN: python teste_klkt.py  --LHC 13 --kl 1 --kt 1 --v1 False
##########################################
# It tests with simulation to: 
#  1.0		1.0	: python function_klkt.py  --LHC 13 --kl 1 --kt 1 --c2 0 --cg 0 --c2g 0 --v1 False
# -10.		0.5	: python function_klkt.py  --LHC 13 --kl -10 --kt 0.5 --c2 0 --cg 0 --c2g 0  --v1 False
#  0.0001	2.25	 : python function_klkt.py  --LHC 13 --kl 0.0001 --kt 2.25 --c2 0 --cg 0 --c2g 0  --v1 False

from array import array
import math
import bisect
from optparse import OptionParser
import numpy as np
import ROOT
import matplotlib
import matplotlib.pyplot as plt

# functions to analytical formula
#from function_klkt import readCoefficients , functionGF , findBin

parser = OptionParser()
parser.add_option("--kl", type="float", dest="kll", help="Multiplicative factor in the H trilinear wrt to SM")
parser.add_option("--kt", type="float", dest="ktt", help="Multiplicative factor in the H top Yukawa wrt to SM")

parser.add_option("--LHC", type="int", dest="lhc", help="LHC CM energy in TeV", default=13)
parser.add_option("--v1", type="int", dest="vv1", help="The version of clustering where the input events are based")
parser.add_option("--c2", type="float", dest="c22", help="ttHH with triangle loop", default=0)
parser.add_option("--cg", type="float", dest="cgg", help="HGG contact", default=0)
parser.add_option("--c2g", type="float", dest="c2gg", help="HHGG contact", default=0)

(options, args) = parser.parse_args()
print "LHC @ %s TeV" % options.lhc
print " "
CM =options.lhc 
kl = options.kll
kt = options.ktt
c2 = options.c22
cg = options.cgg
c2g = options.c2gg
v1 = options.vv1

if v1 == 1 : print "Only the  calculated from the 12 benchmarks defined in 1507.02245v4 (JHEP version) each one with 100k events (v1 == 0)"

if c2 != 0 or cg != 0 or c2g != 0 :  print "The analytical function is not yet implemented"


######################################################
# Declare the function
####################################################
def functionGF(kl,kt,c2,cg,c2g,A):
    return A[0]*kt**4 + (A[1]*kt**2 )*kl**2 + A[2]*kt*kl*kt**2 
#############################################
# read the GF coeficients in 15 2D matrices
############################################# 
def readCoefficients():
  global A13tev
  A13tev = [2.09078, 0.282307, -1.37309] # coefficients for total cross section arxiv:1608.06578
  # tables for store coefficients by bin
  global effSM  
  global effSum  
  global MHH  
  global COSTS  
  global A1 
  global A3  
  global A7  
  effSM = np.zeros((3,13)) 
  effSum = np.zeros((3,13)) 
  MHH = np.zeros((3,13)) 
  COSTS = np.zeros((3,13)) 
  A1 = np.zeros((3,13)) 
  A3 = np.zeros((3,13)) 
  A7 = np.zeros((3,13)) 
  filne = "coefficientsByBin_klkt.txt"
  f = open(filne, 'r+')
  lines = f.readlines() # get all lines as a list (array)
  # Read coefficients by bin
  countercost=0
  countermhh=0
  for line in  lines:
    l = []
    tokens = line.split()
    for token in tokens:
        num = ""
        num_char = "."
        num2 = "e"
        num3 = "-"
        for char in token: 
            if (char.isdigit() or (char in num_char) or (char in num2) or (char in num3)): num = num + char
        try: 
           l.append(float(num))
        except ValueError: pass
    MHH[countercost][countermhh] = l[1] 
    COSTS[countercost][countermhh] = l[2] 
    effSM[countercost][countermhh] = l[3] #/10000. # in units of 10k events
    effSum[countercost][countermhh] = l[4] #/10000. # in units of 10k events # 12 JHEP benchmarks 
    # Just for testing purposes the above contains the number of events by bin from an ensenble of events 
    # calculated from the 12 benchmarks defined in 1507.02245v4 (JHEP version) each one with 100k events
    A1[countercost][countermhh] = l[5]
    A3[countercost][countermhh] = l[6]
    A7[countercost][countermhh] = l[7]
    countercost+=1
    if countercost == 3 :
          countercost=0
          countermhh+=1
  f.close()
################################################
# distribute the calculated GenMHH and CostS in the bins numbering following "coefficientsByBin_klkt.txt"
################################################
def findBin(mhhroot,costhetast):   
   # Note that in root the bin counting in histograms start from 1 (and not 0).
   binGenMHH = [250.,270.,300.,330.,360.,390., 420.,450.,500.,550.,600.,700.,800.,1000.]
   binGenCostS  = [ -1., -0.55,0.55,1.  ]
   global binmhh
   global bincost
   for ii in range (0,13) : 
     if mhhroot > binGenMHH[12-ii] : 
        binmhh = 12-ii 
        break
   for ii in range (0,3) : 
     if costhetast > binGenCostS[2-ii] : 
        bincost = 2-ii
        break
################################################
################################################
# begin of template on how to apply the above functions
################################################
###################################################
###################################################
# Draw the histograms
#####################################################
def plotting(CalcMhh,CalcCost,CalcWeight,CalcMhhTest,CalcCostTest):
  bin_size = 30; min_edge = 250; max_edge = 1000
  N = (max_edge-min_edge)/bin_size; Nplus1 = N + 1
  bin_list = np.linspace(min_edge, max_edge, Nplus1)
  plt.xlim(min_edge, max_edge)
  #plt.yscale('log', nonposy='clip')
  #if kl == 1 and kt == 1 and c2 ==0 and cg == 0 and c2g ==0 : 
  plt.hist(CalcMhhTest, bin_list , normed=1,  histtype='bar', label='simulated', color= 'r', edgecolor='r', lw=3)
  #if kl == -10 and kt == 0.5 and c2 ==0 and cg == 0 and c2g ==0 : plt.hist(CalcMhhTest, bin_list , normed=1,  histtype='bar', label='simulated', color= 'r', edgecolor='r', lw=3)
  #if kl == 0.0001 and kt == 2.25 and c2 ==0 and cg == 0 and c2g ==0 : plt.hist(CalcMhhTest, bin_list , normed=1,  histtype='bar', label='simulated', color= 'r', edgecolor='r', lw=3)
  plt.hist(CalcMhh, bin_list, weights=CalcWeight , normed=1, histtype='bar', label='reweigted', fill=False, color= 'k', edgecolor='k')
  plt.legend(loc='upper right')
  plt.title(" In  kl =="+str(kl)+", kt =="+str(kt)+", c2 =="+str(c2)+", cg =="+str(cg)+", c2g ==" +str(c2g) )
  plt.xlabel("Mhh (GeV)")
  plt.ylabel("events")
  plt.savefig("Mhh_kl_"+str(kl)+"_kt_"+str(kt)+"_c2_"+str(c2)+"_cg_"+str(cg)+"_c2g_" +str(c2g)+".pdf")
  plt.cla()   # Clear axis
  plt.clf()   # Clear figure
  plt.close() # Close a figure window
  bin_size = 0.1; min_edge = -1; max_edge = 1
  N = (max_edge-min_edge)/bin_size; Nplus1 = N + 1
  bin_list = np.linspace(min_edge, max_edge, Nplus1)
  plt.xlim(min_edge, max_edge)
  #plt.yscale('log', nonposy='clip')
  #if kl == 1 and kt == 1 and c2 ==0 and cg == 0 and c2g ==0 : 
  plt.hist(CalcCostTest, bin_list , normed=1,  histtype='bar', label='simulated', color= 'r', edgecolor='r', lw=3)
  #if kl == -10 and kt == 0.5 and c2 ==0 and cg == 0 and c2g ==0 : plt.hist(CalcCostTest, bin_list , normed=1,  histtype='bar', label='simulated', color= 'r', edgecolor='r', lw=3)
  #if kl == 0.0001 and kt == 2.25 and c2 ==0 and cg == 0 and c2g ==0 : plt.hist(CalcCostTest, bin_list , normed=1,  histtype='bar', label='simulated', color= 'r', edgecolor='r', lw=3)
  plt.hist(CalcCost, bin_list, weights=CalcWeight , normed=1, histtype='bar', label='reweigted', fill=False, color= 'k', edgecolor='k', lw=3)
  plt.title(" In  kl =="+str(kl)+"and kt =="+str(kt)+"and c2 =="+str(c2)+"and cg =="+str(cg)+"and c2g ==" +str(c2g) )
  plt.legend(loc='upper right')
  plt.xlabel("Cost*")
  plt.ylabel("events")
  plt.savefig("CostS_kl_"+str(kl)+"_kt_"+str(kt)+"_c2_"+str(c2)+"_cg_"+str(cg)+"_c2g_" +str(c2g)+".pdf")
  plt.cla()   # Clear axis
  plt.clf()   # Clear figure
  plt.close() # Close a figure window
###########################################################
# read events and apply weight
###########################################################
def main():
  # events for testing, in text format 
  pathBSMtest="/afs/cern.ch/work/a/acarvalh/public/toAnamika/GF_HH_toRecursive/" # events of file to superimpose a test
  # see the translation of coefficients for this last on: If you make this script smarter (to only read files we ask to test) you can implement more
  # https://github.com/acarvalh/generateHH/blob/master/fit_GF_HH_lhe/tableToFitA3andA7.txt
  if v1 == 0 : 
    # events to reweights, in text format (for testing only)
    pathBenchEvents="/afs/cern.ch/work/a/acarvalh/public/toAnamika/GF_HH_BSM/" # events to reweight
    # Save numpys with the events from benchmarks 
    CalcMhh = np.zeros((1200000))
    CalcCost = np.zeros((1200000))
    CalcWeight = np.zeros((1200000))
  ##########################################
  # initialize tables of coefficients by bins
  readCoefficients()
  ###################################################
  # Read events from MC simulation in three points to compare 
  if kl == 1 and kt == 1 and c2 ==0 and cg == 0 and c2g ==0 : ntest = 100000
  else : ntest = 50000
  CalcMhhTest = np.zeros((100000))
  CalcCostTest = np.zeros((100000)) # the BSM points have 50000
  CalcMhhTest2 = np.zeros((50000))
  CalcCostTest2 = np.zeros((50000)) # the BSM points have 50000
  CalcMhhTest3 = np.zeros((50000))
  CalcCostTest3 = np.zeros((50000)) # the BSM points have 50000
  ##########################################
  # read the events and fill histograms with weights
  # calculate mhh and cost* and find the bin
  ########################################## 
  countline=0
  countevent=0
  counteventSM=0
  # the numbers > 12 are the 
  countevent13=0 
  countevent14=0
  # read events as text files for events to test 
  Px = np.zeros((2)) 
  Py = np.zeros((2)) 
  Pz = np.zeros((2)) 
  En = np.zeros((2))  
  ##################################################
  for sam in  range(0,15):
   fail =0
   if sam ==0 : filne = pathBenchEvents+"GF_HH_"+str(sam)+".lhe.decayed"    # 0 is SM = here it does not enter in the events to be reweighted
   elif sam < 13 : filne = pathBenchEvents+"GF_HH_"+str(sam)+".lhe.decayed"
   elif sam ==13 : filne = pathBSMtest+"GF_HH_42.lhe.decayed"
   elif sam ==14 : filne = pathBSMtest+"GF_HH_9.lhe.decayed"
   # follow the numbering of: https://github.com/acarvalh/generateHH/blob/master/fit_GF_HH_lhe/tableToFitA3andA7.txt
   f = open(filne, 'r+')
   lines = f.readlines() # get all lines as a list (array)
   for line in  lines:
    l = []
    tokens = line.split()
    for token in tokens:
        num = ""
        num_char = "."
        num2 = "e"
        num3 = "-"
        for char in token: 
            if (char.isdigit() or (char in num_char) or (char in num2) or (char in num3)): num = num + char
        try: 
           l.append(float(num))
        except ValueError: pass
    if countline < 2 :
       Px[countline] = l[1] 
       Py[countline] = l[2]
       Pz[countline] = l[3]
       En[countline] = l[4]
    countline+=1
    if countline==2 :
       # calculate reweigthing 
       if abs(Px[0])!= abs(Px[1]) : print "error parsing ascii file"
       countline=0
       P1 = ROOT.TLorentzVector()
       P1.SetPxPyPzE(Px[0],Py[0],Pz[0],En[0])
       P2 = ROOT.TLorentzVector()
       P1.SetPxPyPzE(Px[1],Py[1],Pz[1],En[1])
       SUM = ROOT.TLorentzVector()
       SUM.SetPxPyPzE(Px[0]+Px[1],Py[0]+Py[1],Pz[0]+Pz[1],En[0]+En[1])
       mhhroot=SUM.M()
       P1boost = P1
       P1boost.Boost(-SUM.BoostVector())
       costhetast = float(P1boost.CosTheta())
       ####################################
       findBin(mhhroot,costhetast)
       #################################
       # Read the events and do histos for the test MC if they are there ( you can implement more ) 
       #################################
       if sam ==0 and kl == 1 and kt == 1 and c2 ==0 and cg == 0 and c2g ==0 :
             CalcMhhTest[counteventSM] = float(mhhroot)
             CalcCostTest[counteventSM] = float(costhetast)
             counteventSM+=1
       if sam ==13 and kl == -10 and kt == 0.5 and c2 ==0 and cg == 0 and c2g ==0  :
             CalcMhhTest2[countevent13] = float(mhhroot)
             CalcCostTest2[countevent13] = float(costhetast)
             countevent13+=1
       if sam ==14 and kl == 0.0001 and kt == 2.25 and c2 ==0 and cg == 0 and c2g ==0  :
             CalcMhhTest3[countevent14] = float(mhhroot)
             CalcCostTest3[countevent14] = float(costhetast)
             countevent14+=1
       #############################################
       # the input events and calculate the weights ===> to adapt to you
       #############################################
       if sam >0 and sam < 13 :
           if effSum[bincost][binmhh] > 0 and A1[bincost][binmhh] > 0: 
             CalcMhh[countevent] = float(mhhroot) # just to draw an histogram
             CalcCost[countevent] = float(costhetast) # just to draw an histogram
             A = [A1[bincost][binmhh],A3[bincost][binmhh],A7[bincost][binmhh]]
             effBSM = float(effSM[bincost][binmhh]*functionGF(kl,kt,c2,cg,c2g,A)/functionGF(kl,kt,c2,cg,c2g,A13tev))
             CalcWeight[countevent] = effBSM/float(effSum[bincost][binmhh]) # ==> JHEP sum
             countevent+=1
   f.close()
  if  kl == 1 and kt == 1 and c2 ==0 and cg == 0 and c2g ==0 :  plotting(CalcMhh,CalcCost,CalcWeight,CalcMhhTest,CalcCostTest)
  if  kl == -10 and kt == 0.5 and c2 ==0 and cg == 0 and c2g ==0 :  plotting(CalcMhh,CalcCost,CalcWeight,CalcMhhTest2,CalcCostTest2)
  if  kl == 0.0001 and kt == 2.25 and c2 ==0 and cg == 0 and c2g ==0 :  plotting(CalcMhh,CalcCost,CalcWeight,CalcMhhTest3,CalcCostTest3)
  print countevent
##########################################
if __name__ == "__main__":  
  main()


