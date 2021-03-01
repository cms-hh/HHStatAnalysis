# Provides scale factors for event reweighting based on an analytical model.
# This file is part of https://github.com/cms-hh/HHStatAnalysis.
# compiling

import ROOT
ROOT.gROOT.SetBatch(True)
import numpy as np
from array import array
import matplotlib
import matplotlib.pyplot as plt
import math

class NonResonantModelNLO:
    def __init__(self):
        # read coefficients from the input file here
        # to store coefficients use self.
        self.NCostHHbin=4
        self.binGenCostS  = [ 0.0, 0.4, 0.6, 0.8, 1.0 ] #

        self.NMHHbin=36
        self.binGenMHH  = [ 250.,   270.,  290.,  310.,  330.,
                            350.,   370.,  390.,  410.,  430.,
                            450.,   470.,  490.,  510.,  530.,
                            550.,   570.,  590.,  610.,  630.,
                            650.,   670.,  700.,  750.,  800.,
                            850.,   900.,  950., 1000., 1100.,
                            1200., 1300., 1400., 1500., 1750., 2000., 5000.]

        self.NCoef=23
        self.A =  [[[0 for MHHbin in range(self.NMHHbin)] for CostHHbin in range(self.NCostHHbin)] for coef in range(self.NCoef)]
        self.A13tev = [62.5088, 345.604, 9.63451, 4.34841, 39.0143, -268.644, -44.2924, 96.5595, 53.515, -155.793, -23.678, 54.5601, 12.2273, -26.8654, -19.3723, -0.0904439, 0.321092, 0.452381, -0.0190758, -0.607163, 1.27408, 0.364487, -0.499263] # from https://github.com/pmandrik/VSEVA/blob/master/HHWWgg/reweight/reweight_HH.C#L117
        self.Cnorm=0
        print "initialize"

    # Declare the function
    def functionGF(self, kl,kt,c2,cg,c2g,A): 
        return A[0] * kt**4 + \
        A[1] * c2**2 + \
        A[2] * kt**2 * kl**2 + \
        A[3] * cg**2 * kl**2 + \
        A[4] * c2g**2 + \
        A[5] * c2 * kt**2 + \
        A[6] * kl * kt**3 + \
        A[7] * kt * kl * c2 + \
        A[8] * cg * kl * c2 + \
        A[9] * c2 * c2g + \
        A[10] * cg * kl * kt**2 + \
        A[11] * c2g * kt**2 + \
        A[12] * kl**2 * cg * kt + \
        A[13] * c2g * kt * kl + \
        A[14] * cg * c2g * kl + \
        A[15] * kt**3 * cg + \
        A[16] * kt * c2 * cg + \
        A[17] * kt * cg**2 * kl + \
        A[18] * cg * kt * c2g + \
        A[19] * kt**2 * cg**2 + \
        A[20] * c2 * cg**2 + \
        A[21] * cg**3 * kl + \
        A[22] * cg**2 * c2g

    def FindBin(self,mhh,cost):
        mhhbin=-1
        if mhh < 250.:
            print "[ERROR]: a value mhh<250 was provided"
        elif mhh >= self.binGenMHH[self.NMHHbin]:
            mhhbin=self.NMHHbin
        else:
            for ibin in range (0,self.NMHHbin):
                if mhh>=self.binGenMHH[ibin] and mhh<self.binGenMHH[ibin+1]:
                    mhhbin=ibin
                    break

        costbin=-1
        if cost<0. or cost>1.:
            print "[ERROR]: a value costheta<0. or costheta>1. was provided"
        else:
            for ibin in range (0,self.NCostHHbin):
                if cost>=self.binGenCostS[ibin] and cost<self.binGenCostS[ibin+1]:
                    costbin=ibin
                    break

        return mhhbin, costbin        

    def ReadCoefficients(self,inputFileName) : #,effSM,MHH,COSTS,A1,A3,A7):
        f = open(inputFileName, 'r')
        lines = f.readlines() # get all lines as a list (array)
        # Read coefficients by bin
        for line in  lines:
            #print "reading line", line
            l = []
            tokens = line.split(",")
            #print "tokens are", tokens
            if tokens[0]=='""': #nominal coefficients have no label
                #print "this is nominal" 
                MHHmin=float(tokens[1])
                costhetamin=float(tokens[3])
                #print "MHHmin, costhetamin = ", MHHmin, costhetamin
                MHHbin, costhetabin = self.FindBin(MHHmin,costhetamin)
                #print "the bin is ",MHHbin, costhetabin
                for coef in range (0,self.NCoef):
                    self.A[coef][costhetabin][MHHbin] = float(tokens[5+coef])
        f.close()
        print "Stored coefficients by bin"

    def getTotalXS(self, kl, kt, c2, cg, c2g):
        return self.functionGF(kl,kt,c2,cg,c2g,self.A13tev)

    def getDifferentialXS2D(self, mhh, cost, kl, kt, c2, cg, c2g):
        #assume XS~0 for mHH > 5000 GeV
        if mhh > self.binGenMHH[ self.NMHHbin ]:
            return 0.
        binmhh, bincost = self.FindBin(mhh,cost)
        dmhh = self.binGenMHH[binmhh+1] - self.binGenMHH[binmhh]
        dcost = self.binGenCostS[bincost+1] - self.binGenCostS[bincost]
        Acoeff = [ self.A[coef][bincost][binmhh] for coef in range(0,self.NCoef) ] 
        return self.functionGF(kl,kt,c2,cg,c2g,Acoeff)/1000. * dcost

    def getDifferentialXSmHH(self, mhh, kl, kt, c2, cg, c2g):
        # Integrate the 2D differental XS over costheta
        dXS=0.
        for bincost in range(0,self.NCostHHbin):
            costhetamin = self.binGenCostS[bincost] 
            costhetamax = self.binGenCostS[bincost+1]
            dcostheta = costhetamax - costhetamin
            costheta = 0.5*(costhetamax+costhetamin)
            dXS += self.getDifferentialXS2D(mhh,costheta,kl,kt,c2,cg,c2g)
        return dXS

    def getBenchmark(self, BM):
        # load benchmarks
        BMcouplings={}
        BMcouplings["kl"]=[1.0,  7.5,  1.0,  1.0,  -3.5, 1.0, 2.4, 5.0, 15.0, 1.0, 10.0, 2.4, 15.0]
        BMcouplings["kt"]=[1.0,  1.0,  1.0,  1.0,  1.5,  1.0, 1.0, 1.0, 1.0,  1.0, 1.5,  1.0, 1.0]
        BMcouplings["c2"]=[0.0,  -1.0, 0.5, -1.5, -3.0,  0.0, 0.0, 0.0, 0.0,  1.0, -1.0, 0.0, 1.0]
        BMcouplings["cg"]=[0.0,  0.0, -0.8,  0.0, 0.0,   0.8, 0.2, 0.2, -1.0, -0.6, 0.0, 1.0, 0.0]
        BMcouplings["c2g"]=[0.0, 0.0, 0.6, -0.8, 0.0, -1.0, -0.2,-0.2,  1.0,  0.6, 0.0, -1.0, 0.0] 
        return [ BMcouplings[param][BM] for param in ["kl", "kt", "c2", "cg", "c2g"] ]


    def CalculateMhhCost(self,mhhcost,countline,Px,Py,Pz,En) :
        # calculate reweigthing 
        if abs(Px[0])!= abs(Px[1]) : print "error parsing ascii file"
        P1 = ROOT.TLorentzVector()
        P1.SetPxPyPzE(Px[0],Py[0],Pz[0],En[0])    
        P2 = ROOT.TLorentzVector()
        P1.SetPxPyPzE(Px[1],Py[1],Pz[1],En[1])
        SUM = ROOT.TLorentzVector()
        SUM.SetPxPyPzE(Px[0]+Px[1],Py[0]+Py[1],Pz[0]+Pz[1],En[0]+En[1])
        mhhcost[0]=SUM.M()
        P1boost = P1
        P1boost.Boost(-SUM.BoostVector())
        mhhcost[1] = float(P1boost.CosTheta())
        mhhcost[2] = float(P1.Pt())
        mhhcost[3] = float(SUM.Pt())
