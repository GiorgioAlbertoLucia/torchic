from ROOT import TH2F, RooRealVar, TFile
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
    roofitter = Roofitter(x, ['gaus', 'exp']) # fit with a gaussian and an exponential bkg

    for pt_bin in range(h2_nsigmaTPC.GetXaxis().FindBin(xmin), h2_nsigmaTPC.GetXaxis().FindBin(xmax)):
        h_nsigmaTOF = h2_nsigmaTPC.ProjectionY(f'h_nsigmaTPC_{pt_bin}', pt_bin, pt_bin, 'e')
        pt = h2_nsigmaTPC.GetXaxis().GetBinCenter(pt_bin)

        fractions = roofitter.fit(h_nsigmaTOF, -5., 5., funcs_to_fit=['gaus_0', 'exp_1'])
        roofitter.plot(output_file, canvas_name=f'cNSigmaTPC_pt{pt}')

    h2_nsigmaTPC.Write()
    output_file.Close()


if __name__ == '__main__':
    purity_proton()