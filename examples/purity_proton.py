import numpy as np
from ROOT import TH2F, RooRealVar, TFile, TGraph
from torchic import HistLoadInfo
from torchic import histogram
from torchic.core.roofitter import Roofitter

def purity_proton():

    h2_nsigmaTPC_info = HistLoadInfo('data/AnalysisResults_LiPbPb.root', 'lithium4analysis/QA/h2NsigmaProtonTPC_preselection')
    h2_nsigmaTPC = histogram.load_hist(h2_nsigmaTPC_info)
    xmin = 0.
    xmax = 5.
    ymin = -5.
    ymax = 5.

    output_file = TFile('data/purity_proton.root', 'recreate')

    x = RooRealVar('x', 'x', ymin, ymax)
    roofitter = Roofitter(x, ['gaus', 'exp_offset']) # fit with a gaussian and an exponential bkg

    pts = []
    signal = []
    bkg = []    

    for pt_bin in range(h2_nsigmaTPC.GetXaxis().FindBin(xmin), h2_nsigmaTPC.GetXaxis().FindBin(xmax)):
        h_nsigmaTOF = h2_nsigmaTPC.ProjectionY(f'h_nsigmaTPC_{pt_bin}', pt_bin, pt_bin, 'e')
        pt = h2_nsigmaTPC.GetXaxis().GetBinCenter(pt_bin)

        fractions = roofitter.fit(h_nsigmaTOF, -5., 5., funcs_to_fit=['gaus_0', 'exp_1'])
        pts.append(pt)
        signal.append(fractions[0].getValV())
        bkg.append(fractions[1].getValV())
        print('fraction type:', type(fractions[0]))
        print('fraction getVal type:', type(fractions[0].getValV()))
        print('fraction getVal:', fractions[0].getValV())

        roofitter.plot(output_file, canvas_name=f'cNSigmaTPC_pt{pt}')

    graph_signal = TGraph(len(pts), np.array(pts, dtype=np.float32), np.array(signal, dtype=np.float32))
    graph_signal.SetName('g_signal')
    graph_bkg = TGraph(len(pts), np.array(pts, dtype=np.float32), np.array(bkg, dtype=np.float32))
    graph_bkg.SetName('g_bkg')

    graph_signal.Write()
    graph_bkg.Write()
    h2_nsigmaTPC.Write()
    output_file.Close()


if __name__ == '__main__':
    purity_proton()