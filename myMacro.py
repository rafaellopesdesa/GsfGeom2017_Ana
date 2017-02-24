from ROOT import *
gROOT.SetBatch(True)
from array import *
import numpy as np

import sys
import CMS_lumi, tdrstyle

iPos    = 0
iPeriod = 0
H_ref = 600; 
W_ref = 800; 

etaCategories = [0.0,1.0,1.444,1.556,2.0,2.5]
colorSeven = [TColor.GetColor(244,109,67),
              TColor.GetColor(253,174,97),
              TColor.GetColor(254,224,144),
              TColor.GetColor(224,243,248),
              TColor.GetColor(171,217,233),
              TColor.GetColor(116,173,209),
              TColor.GetColor(69,117,180)]
colorTwo = [TColor.GetColor(255,255,204), TColor.GetColor(12,44,132)]
seedNames = ['initialStepSeeds',
             'highPtTripletStepSeeds',
             'mixedTripletStepSeeds',
             'pixelLessStepSeeds',
             'tripletElectronSeeds',
             'pixelPairElectronSeeds',
             'stripPairElectronSeeds']
seedTypes = ['ECALdriven', 'TRKdriven']

class EffHistograms:
    
    def __init__(self, name, binning, categories):
        self.name = name
        self.categories = categories
        self.allHist = TH1D('all'+name, '', len(binning)-1, array('d', binning))

        self.catHists = []
        for cat in categories:
            self.catHists.append(TH1D(cat+name,'', len(binning)-1, array('d', binning)))

    def Fill(self, var, eta, cat, wgt = 1):
        self.allHist.Fill(var)
        self.catHists[cat].Fill(var, wgt)

    def MakeLegend(self, x1, y1, x2, y2, colors):
        legend = TLegend(x1, y1, x2, y2)
        if colors is not None:
            legend.SetFillColor(colors[0])
        for ih, hist in enumerate(self.catHists):
            legend.AddEntry(hist, '%s [%.1f%%]' % (self.categories[ih], 100*hist.Integral()/self.allHist.Integral()), 'f')
        return legend
    
    def MakeStack(self, title, colors, etaCategory = None):
        stack = THStack('stack'+self.name, title)
        for color, hist in zip(colors, self.catHists):
            hist.SetLineColor(kBlack)
            hist.SetFillColor(color)
            hist.Divide(self.allHist)
            stack.Add(hist, "hist")
        return stack

def getBin(val, cats):
    cat = int(np.digitize([val],cats)[0])-1
    if cat > len(cats)-2: cat = len(cats)-2
    return cat

# common configuration for texts
def setupGraphics():

    tdrstyle.setTDRStyle()
    CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"
    CMS_lumi.lumi_8TeV = "18.3 fb^{-1}"
    CMS_lumi.lumi_13TeV = "36.3 fb^{-1}"
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText  = "(Simulation)"
    CMS_lumi.extraText2 = "Z#rightarrow ee"
    CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)    
    if( iPos==0 ): CMS_lumi.relPosX = 0.12

# common configuration for hists
def setupHists(h):
    h.Draw()
    h.GetXaxis().SetTitleSize(0.75*h.GetXaxis().GetTitleSize())
    h.GetYaxis().SetTitleSize(0.75*h.GetYaxis().GetTitleSize())
    h.GetXaxis().SetTitleOffset(1.3*h.GetXaxis().GetTitleOffset())
     
# common configuration for canvases
def setupCanvas(c):
    W = W_ref
    H  = H_ref
    T = 0.08*H_ref
    B = 0.12*H_ref 
    L = 0.12*W_ref
    R = 0.04*W_ref
    c.SetGridx(1)
    c.SetGridy(1)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetLeftMargin( L/W )
    c.SetRightMargin( R/W )
    c.SetTopMargin( T/H )
    c.SetBottomMargin( B/H )
    c.SetTickx(0)
    c.SetTicky(0)

#common configuration for legends
def setupLegend(l):
    l.SetBorderSize(0)
    l.SetTextFont(42)
    l.SetTextSize(1.2*l.GetTextSize())

def printHist(hist,legend = None):

    canvas = TCanvas('canvas_%s'%hist.GetName(),'Canvas %s'%hist.GetTitle(),50,50,W_ref,H_ref)
    canvas.cd()
    setupCanvas(canvas)
    setupHists(hist)
    hist.Draw()
    CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)
    canvas.cd()
    canvas.Update()
    canvas.RedrawAxis()    
    if legend is not None:
        setupLegend(legend)
        legend.Draw()
    return canvas

def analyzeSeeds(chain):

    # Define the binning to be logarithmic uniform
    ptBinning = [0.5]
    etaBinning = [-2.5]
    for i in xrange(1,51):
        ptBinning.append(TMath.Power(10,TMath.Log10(5)-1+i*(TMath.Log10(2)+2.-TMath.Log10(5)+1.)/50))
        etaBinning.append(-2.5+(i*5.)/50.)
        
    # Add histograms
    ptECALdrivenSeedHist = EffHistograms('ptECALdrivenSeedHist', ptBinning, seedNames)
    etaECALdrivenSeedHist = EffHistograms('etaECALdrivenSeedHist', etaBinning, seedNames)
    ptALLSeedHist = EffHistograms('ptALLSeedHist', ptBinning, seedTypes)
    etaALLSeedHist = EffHistograms('etaALLSeedHist', etaBinning, seedTypes)

    # In eta categories
    ptALLSeedHistEtaCategories = [EffHistograms('ptALLSeedHistEta1', ptBinning, seedTypes),
                                  EffHistograms('ptALLSeedHistEta2', ptBinning, seedTypes),
                                  EffHistograms('ptALLSeedHistEta3', ptBinning, seedTypes),
                                  EffHistograms('ptALLSeedHistEta4', ptBinning, seedTypes),
                                  EffHistograms('ptALLSeedHistEta5', ptBinning, seedTypes)]

    nentries = chain.GetEntries()
    for ientry, events in enumerate(chain):
        if ientry % 1000 == 0:
            print '%d / %d' % (ientry, nentries)
        usedSimIdx = [] # Apparently we have an irritating duplication of tracks in electronGsfTracks.
                        # I should fix that in the next version 
        for seedIdx, simIdxs in zip(events.trk_seedIdx, events.trk_simTrkIdx):
            if len(simIdxs) > 0:                                                     # sim-matched 
                if usedSimIdx.count(simIdxs[0]) > 0:                                 # discard duplicated entries
                    continue
                if TMath.Abs(events.sim_pdgId[simIdxs[0]]) == 11:                    # matched to electrons
                    usedSimIdx.append(simIdxs[0])
                    simVec = Math.PxPyPzMVector(events.sim_px[simIdxs[0]],
                                                events.sim_py[simIdxs[0]],
                                                events.sim_pz[simIdxs[0]],
                                                0.000511)
                    
                    etaCategory = getBin(TMath.Abs(simVec.eta()),etaCategories)
                    if events.see_ecalDriven[seedIdx]:                               
                        ptALLSeedHist.Fill(simVec.pt(), simVec.eta(), 0)
                        ptALLSeedHistEtaCategories[etaCategory].Fill(simVec.pt(), simVec.eta(), 0)
                        if simVec.pt() > 10.: etaALLSeedHist.Fill(simVec.eta(), simVec.eta(), 0)
                        algo = events.see_algoOriginal[seedIdx]
                        ptECALdrivenSeedHist.Fill(simVec.pt(), simVec.eta(), algo)
                        if simVec.pt() > 10.: etaECALdrivenSeedHist.Fill(simVec.eta(), simVec.eta(), algo)
                    else:
                        ptALLSeedHist.Fill(simVec.pt(), simVec.eta(), 1)
                        ptALLSeedHistEtaCategories[etaCategory].Fill(simVec.pt(), simVec.eta(), 1)
                        if simVec.pt() > 10.: etaALLSeedHist.Fill(simVec.eta(), simVec.eta(), 1)

    ptECALdrivenLegend = ptECALdrivenSeedHist.MakeLegend(0.55, 0.15, 0.85, 0.55, colorSeven)
    etaECALdrivenLegend = etaECALdrivenSeedHist.MakeLegend(0.55, 0.15, 0.85, 0.55, colorSeven)

    ptECALdrivenStack  = ptECALdrivenSeedHist.MakeStack(';Simulated Transverse Momentum [GeV];Fraction of Reconstructed GsfElectrons', colorSeven)
    etaECALdrivenStack = etaECALdrivenSeedHist.MakeStack(';Simulated Pseudorapidity;Fraction of Reconstructed GsfElectrons', colorSeven)
        
    ptECALdrivenStack.SetMaximum(1)
    etaECALdrivenStack.SetMaximum(1)

    ptALLLegend = ptALLSeedHist.MakeLegend(0.55, 0.15, 0.85, 0.35, colorTwo)
    etaALLLegend = etaALLSeedHist.MakeLegend(0.55, 0.15, 0.85, 0.35, colorTwo)

    ptALLStack  = ptALLSeedHist.MakeStack(';Simulated Transverse Momentum [GeV];Fraction of Reconstructed GsfElectrons', colorTwo)
    etaALLStack = etaALLSeedHist.MakeStack(';Simulated Pseudorapidity;Fraction of Reconstructed GsfElectrons', colorTwo)
        
    ptALLStack.SetMaximum(1)
    etaALLStack.SetMaximum(1)

    ptALLLegendEta0 = ptALLSeedHistEtaCategories[0].MakeLegend(0.55, 0.15, 0.85, 0.35, colorTwo)
    ptALLStackEta0  = ptALLSeedHistEtaCategories[0].MakeStack(';Simulated Transverse Momentum [GeV];Fraction of Reconstructed GsfElectrons', colorTwo)
    ptALLStackEta0.SetMaximum(1)

    ptALLLegendEta1 = ptALLSeedHistEtaCategories[1].MakeLegend(0.55, 0.15, 0.85, 0.35, colorTwo)
    ptALLStackEta1  = ptALLSeedHistEtaCategories[1].MakeStack(';Simulated Transverse Momentum [GeV];Fraction of Reconstructed GsfElectrons', colorTwo)
    ptALLStackEta1.SetMaximum(1)

    ptALLLegendEta2 = ptALLSeedHistEtaCategories[2].MakeLegend(0.55, 0.15, 0.85, 0.35, colorTwo)
    ptALLStackEta2  = ptALLSeedHistEtaCategories[2].MakeStack(';Simulated Transverse Momentum [GeV];Fraction of Reconstructed GsfElectrons', colorTwo)
    ptALLStackEta2.SetMaximum(1)

    ptALLLegendEta3 = ptALLSeedHistEtaCategories[3].MakeLegend(0.55, 0.15, 0.85, 0.35, colorTwo)
    ptALLStackEta3  = ptALLSeedHistEtaCategories[3].MakeStack(';Simulated Transverse Momentum [GeV];Fraction of Reconstructed GsfElectrons', colorTwo)
    ptALLStackEta2.SetMaximum(1)

    ptALLLegendEta4 = ptALLSeedHistEtaCategories[4].MakeLegend(0.55, 0.15, 0.85, 0.35, colorTwo)
    ptALLStackEta4  = ptALLSeedHistEtaCategories[4].MakeStack(';Simulated Transverse Momentum [GeV];Fraction of Reconstructed GsfElectrons', colorTwo)
    ptALLStackEta2.SetMaximum(1)

    return (ptECALdrivenStack, ptECALdrivenLegend, etaECALdrivenStack, etaECALdrivenLegend, ptALLStack, ptALLLegend, etaALLStack, etaALLLegend, ptALLStackEta0, ptALLLegendEta0, ptALLStackEta1, ptALLLegendEta1, ptALLStackEta2, ptALLLegendEta2, ptALLStackEta3, ptALLLegendEta3, ptALLStackEta4, ptALLLegendEta4)

def analyzeECALdrivenEfficiency(chain):

    # Define the binning to be logarithmic uniform
    ptBinning = [0.5]
    etaBinning = [-2.5]
    for i in xrange(1,51):
        ptBinning.append(TMath.Power(10,TMath.Log10(5)-1+i*(TMath.Log10(2)+2.-TMath.Log10(5)+1.)/50))
        etaBinning.append(-2.5+(i*5.)/50.)
        
    # Add histograms
    ptECALdrivenSeedEff = EffHistograms('ptECALdrivenSeedEff', ptBinning, seedNames)
    etaECALdrivenSeedEff = EffHistograms('etaECALdrivenSeedEff', etaBinning, seedNames)

    nentries = chain.GetEntries()
    for ientry, events in enumerate(chain):
        if ientry % 1000 == 0:
            print '%d / %d' % (ientry, nentries)
        for trkIdxs, shareFracs, px, py, pz, pdgId, vtxIdx in zip(events.sim_trkIdx, events.sim_shareFrac, events.sim_px, events.sim_py, events.sim_pz, events.sim_pdgId, events.sim_parentVtxIdx):
            if TMath.Abs(pdgId) == 11 and len(events.simvtx_sourceSimIdx[vtxIdx]) == 0:             # an electron that comes from the pp interaction
                simVec = Math.PxPyPzMVector(px, py, pz, 0.000511)
                if TMath.Abs(simVec.eta()) > 2.5: continue
                wgt = 0.
                algo = 0
                if len(trkIdxs) > 0:
                    seedIdx = events.trk_seedIdx[trkIdxs[0]]
                    if events.see_ecalDriven[seedIdx]:                 
                        wgt = 1.
                        algo = events.see_algoOriginal[seedIdx]
                ptECALdrivenSeedEff.Fill(simVec.pt(), simVec.eta(), algo, wgt)
                if simVec.pt() > 10: etaECALdrivenSeedEff.Fill(simVec.eta(), simVec.eta(), algo, wgt)

    ptECALdrivenEffLegend = ptECALdrivenSeedEff.MakeLegend(0.15, 0.49, 0.45, 0.89, None)
    etaECALdrivenEffLegend = etaECALdrivenSeedEff.MakeLegend(0.55, 0.15, 0.85, 0.55, colorSeven)

    ptECALdrivenEffStack  = ptECALdrivenSeedEff.MakeStack(';Simulated Transverse Momentum [GeV];Fraction of Simulated GsfElectrons', colorSeven)
    etaECALdrivenEffStack = etaECALdrivenSeedEff.MakeStack(';Simulated Pseudorapidity;Fraction of Simulated GsfElectrons', colorSeven)
    
    ptECALdrivenEffStack.SetMaximum(1)
    etaECALdrivenEffStack.SetMaximum(1)

    return (ptECALdrivenEffStack, ptECALdrivenEffLegend, etaECALdrivenEffStack, etaECALdrivenEffLegend)

    

if __name__ == '__main__':

    inputFiles = ['/eos/uscms/store/user/rclsa/GsfTrackingNtuple/ValTrkGSF_1_0/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/trackingNtuple_%d.root' % x for x in xrange(1,786)]
    inputChain = TChain('trackingNtuple/tree')
    for inputFile in inputFiles:
        inputChain.Add(inputFile)

    setupGraphics()

    ptECALdrivenEffStack, ptECALdrivenEffLegend, etaECALdrivenEffStack, etaECALdrivenEffLegend = analyzeECALdrivenEfficiency(inputChain)
    ptECALdrivenEffCanvas = printHist(ptECALdrivenEffStack, ptECALdrivenEffLegend)
    ptECALdrivenEffCanvas.SetLogx(1)
    ptECALdrivenEffCanvas.Update()
    ptECALdrivenEffCanvas.Print('ptECALdrivenEffHist.png')
    ptECALdrivenEffCanvas.Print('ptECALdrivenEffHist.pdf')
    etaECALdrivenEffCanvas = printHist(etaECALdrivenEffStack, etaECALdrivenEffLegend)
    etaECALdrivenEffCanvas.Print('etaECALdrivenEffHist.png')
    etaECALdrivenEffCanvas.Print('etaECALdrivenEffHist.pdf')

    ptECALdrivenStack, ptECALdrivenLegend, etaECALdrivenStack, etaECALdrivenLegend, ptALLStack, ptALLLegend, etaALLStack, etaALLLegend, ptALLStackEta0, ptALLLegendEta0, ptALLStackEta1, ptALLLegendEta1, ptALLStackEta2, ptALLLegendEta2, ptALLStackEta3, ptALLLegendEta3, ptALLStackEta4, ptALLLegendEta4 = analyzeSeeds(inputChain)
    ptECALdrivenCanvas = printHist(ptECALdrivenStack, ptECALdrivenLegend)
    ptECALdrivenCanvas.SetLogx(1)
    ptECALdrivenCanvas.Update()
    ptECALdrivenCanvas.Print('ptECALdrivenHist.png')
    ptECALdrivenCanvas.Print('ptECALdrivenHist.pdf')
    etaECALdrivenCanvas = printHist(etaECALdrivenStack, etaECALdrivenLegend)
    etaECALdrivenCanvas.Print('etaECALdrivenHist.png')
    etaECALdrivenCanvas.Print('etaECALdrivenHist.pdf')
    ptALLCanvas = printHist(ptALLStack, ptALLLegend)
    ptALLCanvas.SetLogx(1)
    ptALLCanvas.Update()
    ptALLCanvas.Print('ptALLHist.png')
    ptALLCanvas.Print('ptALLHist.pdf')
    etaALLCanvas = printHist(etaALLStack, etaALLLegend)
    etaALLCanvas.Print('etaALLHist.png')
    etaALLCanvas.Print('etaALLHist.pdf')

    ptALLCanvasEta0 = printHist(ptALLStackEta0, ptALLLegendEta0)
    ptALLCanvasEta0.SetLogx(1)
    ptALLCanvasEta0.Update()
    ptALLCanvasEta0.Print('ptALLHistEta0.png')
    ptALLCanvasEta0.Print('ptALLHistEta0.pdf')
    ptALLCanvasEta1 = printHist(ptALLStackEta1, ptALLLegendEta1)
    ptALLCanvasEta1.SetLogx(1)
    ptALLCanvasEta1.Update()
    ptALLCanvasEta1.Print('ptALLHistEta1.png')
    ptALLCanvasEta1.Print('ptALLHistEta1.pdf')
    ptALLCanvasEta2 = printHist(ptALLStackEta2, ptALLLegendEta2)
    ptALLCanvasEta2.SetLogx(1)
    ptALLCanvasEta2.Update()
    ptALLCanvasEta2.Print('ptALLHistEta2.png')
    ptALLCanvasEta2.Print('ptALLHistEta2.pdf')
    ptALLCanvasEta3 = printHist(ptALLStackEta3, ptALLLegendEta3)
    ptALLCanvasEta3.SetLogx(1)
    ptALLCanvasEta3.Update()
    ptALLCanvasEta3.Print('ptALLHistEta3.png')
    ptALLCanvasEta3.Print('ptALLHistEta3.pdf')
    ptALLCanvasEta4 = printHist(ptALLStackEta4, ptALLLegendEta4)
    ptALLCanvasEta4.SetLogx(1)
    ptALLCanvasEta4.Update()
    ptALLCanvasEta4.Print('ptALLHistEta4.png')
    ptALLCanvasEta4.Print('ptALLHistEta4.pdf')

    
