# 
# This class provides a tool to access the differential XS(mHH, costhetaHH) 
# and total XS of the HH @ NLO through an analytical formula.
# The coefficients were derived by P. Mandrik and originally stored in the 
# repository https://github.com/pmandrik/VSEVA/tree/master/HHWWgg/reweight  
#
# For an example on how to use the reweight see the script usage_example.py
# in the same folder of this file

import ROOT
ROOT.gROOT.SetBatch(True)
import numpy as np
from array import array
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
            raise Exception("[ERROR]: a value mhh<250 was provided")
        elif mhh >= self.binGenMHH[self.NMHHbin]:
            mhhbin=self.NMHHbin-1
        else:
            for ibin in range (0,self.NMHHbin):
                if mhh>=self.binGenMHH[ibin] and mhh<self.binGenMHH[ibin+1]:
                    mhhbin=ibin
                    break

        costbin=-1
        if cost<0. or cost>1.:
            raise Exception("costheta<0 or costheta>1 was provided")
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
            l = []
            tokens = line.split(",")
            if tokens[0]=='""': #nominal coefficients have no label
                MHHmin=float(tokens[1])
                costhetamin=float(tokens[3])
                MHHbin, costhetabin = self.FindBin(MHHmin,costhetamin)
                for coef in range (0,self.NCoef):
                    self.A[coef][costhetabin][MHHbin] = float(tokens[5+coef])
        f.close()
        print "Stored coefficients by bin"

    def getTotalXS(self, kl, kt, c2, cg, c2g):
        return self.functionGF(kl,kt,c2,cg,c2g,self.A13tev)

    def getmHHbinwidth(self,mHHvalue):
        binmhh, bincost = self.FindBin(mHHvalue,0.)
        dmhh = self.binGenMHH[binmhh+1] - self.binGenMHH[binmhh]
        return dmhh

    def getcosthetabinwidth(self,costhetavalue):
        binmhh, bincostheta = self.FindBin(500,costhetavalue)
        dcostheta = self.binGenCostS[bincostheta+1] - self.binGenCostS[bincostheta]
        return dcostheta

    def getDifferentialXS2D(self, mhh, cost, kl, kt, c2, cg, c2g):
        #assume XS~0 for mHH > 5000 GeV
        if mhh >= self.binGenMHH[ self.NMHHbin ]:
            return 0.
        binmhh, bincost = self.FindBin(mhh,cost)
        Acoeff = [ self.A[coef][bincost][binmhh] for coef in range(0,self.NCoef) ] 
        return self.functionGF(kl,kt,c2,cg,c2g,Acoeff)/1000.

    def getDifferentialXSmHH(self, mhh, kl, kt, c2, cg, c2g):
        # Integrate the 2D differental XS over costheta
        dXS=0.
        for bincost in range(0,self.NCostHHbin):
            costhetamin = self.binGenCostS[bincost] 
            costhetamax = self.binGenCostS[bincost+1]
            dcostheta = costhetamax - costhetamin
            costheta = 0.5*(costhetamax+costhetamin)
            dXS += self.getDifferentialXS2D(mhh,costheta,kl,kt,c2,cg,c2g) * dcostheta
        return dXS

    def getBenchmark(self, BM): #benchmark points from JHEP04(2016)126
        # load benchmarks
        BMcouplings={}
        BMcouplings["kl"]=[1.0,  7.5,  1.0,  1.0,  -3.5, 1.0, 2.4, 5.0, 15.0, 1.0, 10.0, 2.4, 15.0]
        BMcouplings["kt"]=[1.0,  1.0,  1.0,  1.0,  1.5,  1.0, 1.0, 1.0, 1.0,  1.0, 1.5,  1.0, 1.0]
        BMcouplings["c2"]=[0.0,  -1.0, 0.5, -1.5, -3.0,  0.0, 0.0, 0.0, 0.0,  1.0, -1.0, 0.0, 1.0]
        BMcouplings["cg"]=[0.0,  0.0, -0.8,  0.0, 0.0,   0.8, 0.2, 0.2, -1.0, -0.6, 0.0, 1.0, 0.0]
        BMcouplings["c2g"]=[0.0, 0.0, 0.6, -0.8, 0.0, -1.0, -0.2,-0.2,  1.0,  0.6, 0.0, -1.0, 0.0] 
        return [ BMcouplings[param][BM] for param in ["kl", "kt", "c2", "cg", "c2g"] ]

    def getBenchmark8a(self): # A new benchmark is defined in JHEP09(2018)057 to re-produce the expected features of benchmark 8 
        return [1.0, 1.0, 0.5, 0.8/3, 0.0]

    def getBenchmark2020(self, BM): #new benchmarks from JHEP03(2020)091
        BMcouplings={}
        BMcouplings["kl"]=[   3.94,   6.84,    2.21,   2.79,   3.95,  5.68, -0.10]
        BMcouplings["kt"]=[   0.94,   0.61,    1.05,   0.61,   1.17,  0.83,  0.94]
        BMcouplings["c2"]=[ -1./3.,  1./3.,  -1./3.,  1./3., -1./3., 1./3.,  1.]
        BMcouplings["cg"]=[    0.5,    0.0,     0.5,   -0.5,  1./6.,  -0.5,  1./6.]
        BMcouplings["c2g"]=[ 1./3., -1./3.,    0.5 ,  1./6.,  -0.5 , 1./3., -1./6.]
        #a redefinition of the cg and c2g parameters was done in the paper
        for ivalue in range(len(BMcouplings["cg"])):
            BMcouplings["cg"][ivalue] *= 1.5
        for ivalue in range(len(BMcouplings["c2g"])):
            BMcouplings["c2g"][ivalue] *= (-3.)

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

