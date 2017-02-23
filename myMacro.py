from ROOT import *
from array import *

import CMS_lumi, tdrstyle
iPos    = 11
iPeriod = 0

def setupGraphics():

    tdrstyle.setTDRStyle()
    CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"
    CMS_lumi.lumi_8TeV = "18.3 fb^{-1}"
    CMS_lumi.lumi_13TeV = "36.3 fb^{-1}"
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "(Simulation)"
    CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)    
    if( iPos==0 ): CMS_lumi.relPosX = 0.12
    
def printHist(h,legend = None):

    H_ref = 600; 
    W_ref = 800; 
    W = W_ref
    H  = H_ref
    T = 0.08*H_ref
    B = 0.12*H_ref 
    L = 0.12*W_ref
    R = 0.04*W_ref
    canvas = TCanvas("c2","c2",50,50,W,H)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin( L/W )
    canvas.SetRightMargin( R/W )
    canvas.SetTopMargin( T/H )
    canvas.SetBottomMargin( B/H )
    canvas.SetTickx(0)
    canvas.SetTicky(0)

    h.Draw()
    CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)
    canvas.cd()
    canvas.Update()
    canvas.RedrawAxis()
    if legend is not None:
        legend.Draw()
    return canvas

def analyzeSeeds(chain):

    seedName = ['initialStepSeeds',
                'highPtTripletStepSeeds',
                'mixedTripletStepSeeds',
                'pixelLessStepSeeds',
                'tripletElectronSeeds',
                'pixelPairElectronSeeds',
                'stripPairElectronSeeds']
    binning = [0.5]
    for i in xrange(1,50):
        binning.append(TMath.Power(10,TMath.Log10(5)-1+i*(TMath.Log10(2)+2.-TMath.Log10(5)+1.)/50))
    seedHists = []
    for x in seedName:
        seedHists.append(TH1D(x,x,len(binning)-1,array('d',binning)))
    seedHists.append(TH1D('total','total',len(binning)-1,array('d',binning)))
    nentries = chain.GetEntries()
    for ientry, events in enumerate(chain):
        if ientry % 1000 == 0:
            print '%d / %d' % (ientry, nentries)
        usedSimIdx = [] # Apparently we have an irritating duplication of tracks in electronGsfTracks.
                        # I should fix that in the next version 
        for seedIdx, simIdxs in zip(events.trk_seedIdx, events.trk_simTrkIdx):
            if len(simIdxs) > 0:                                                     # matched 
                if usedSimIdx.count(simIdxs[0]) > 0:
                    continue
                if events.see_ecalDriven[seedIdx]:                                       # ECAL-driven
                    if TMath.Abs(events.sim_pdgId[simIdxs[0]]) == 11:                 # matched to electrons
                        usedSimIdx.append(simIdxs[0])
                        algo = events.see_algoOriginal[seedIdx]
                        simVec = Math.PxPyPzMVector(events.sim_px[simIdxs[0]],
                                                       events.sim_py[simIdxs[0]],
                                                       events.sim_pz[simIdxs[0]],
                                                       0.000511)
                        seedHists[algo].Fill(simVec.pt())
                        seedHists[-1].Fill(simVec.pt())
    stack = THStack('seeds', 'Electron seeds')
    colors = [kViolet, kBlue, kCyan, kGreen, kYellow, kOrange, kRed]    

    legend = TLegend(0.55, 0.65, 0.8, 0.90)
    legend.SetBorderSize(0)
    for ih, hist in enumerate(seedHists[:-1]):
        legend.AddEntry(hist, '%s [%.2f%%]' % (seedName[ih], hist.Integral()/seedHists[-1].Integral()), 'f')

    for color, hist in zip(colors, seedHists[:-1]):
        hist.SetLineColor(kBlack)
        hist.SetFillColor(color)
        hist.Divide(seedHists[-1])
        stack.Add(hist)
    stack.SetMaximum(1.5)
    stack.Draw()
    stack.GetXaxis().SetTitle('True Transverse Momentum [GeV]')
    stack.GetYaxis().SetTitle('Fraction of Reconstructed GsfElectrons')
    return stack, legend
                    
if __name__ == '__main__':

    inputFiles = ['/eos/uscms/store/user/rclsa/GsfTrackingNtuple/ValTrkGSF_1_0/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/trackingNtuple_%d.root' % x for x in xrange(1,786)]
    inputChain = TChain('trackingNtuple/tree')
    for inputFile in inputFiles:
        inputChain.Add(inputFile)

    setupGraphics()

    stackSeeds, legendSeeds = analyzeSeeds(inputChain)
    canvasSeeds = printHist(stack, legend)
    canvasSeeds.SetLogx(1)
    canvasSeeds.Print('ECALdrivenSeeds.png')
