from ROOT import *
from tdrStyle import *
setTDRStyle()

import os,sys,glob
from array import array

r_exp = array('d',[])
nll_exp = array('d',[])
r_obs = array('d',[])
nll_obs = array('d',[])
zeros = array('d',[])

f_exp = TFile("higgsCombineExp.MultiDimFit.mH120.root","READ")
t_exp = f_exp.Get("limit")
for i in xrange(1,t_exp.GetEntries()):
    t_exp.GetEntry(i)
    r_exp.append(t_exp.r)
    nll_exp.append(2.0*t_exp.deltaNLL)
    zeros.append(0.0)
#f_obs = TFile("higgsCombineObs.MultiDimFit.mH125.root","READ")
#t_obs = f_obs.Get("limit")
#for i in xrange(1,t_obs.GetEntries()):
#    t_obs.GetEntry(i)#
#    r_obs.append(t_obs.r)
#    nll_obs.append(2.0*t_obs.deltaNLL)
#    zeros.append(0.0)

v_r_exp = TVectorD(len(r_exp),r_exp)
#v_r_obs = TVectorD(len(r_obs),r_obs)
v_nll_exp = TVectorD(len(nll_exp),nll_exp)
#v_nll_obs = TVectorD(len(nll_obs),nll_obs)
v_zeros = TVectorD(len(zeros),zeros)

c = TCanvas("c","c",800, 800)
c.SetRightMargin(0.06)
c.SetLeftMargin(0.2)

dummy = TH1D("dummy","dummy", 1, 0.0,10.)
dummy.SetBinContent(1,0.0)
dummy.GetXaxis().SetTitle('#mu_{off-shell}')
dummy.GetYaxis().SetTitle('-2 #Delta lnL')
dummy.SetLineColor(0)
dummy.SetLineWidth(0)
dummy.SetFillColor(0)
dummy.SetMinimum(0.0)
dummy.SetMaximum(5.0)
dummy.Draw()

latexf = TLatex()
latexf.SetTextSize(0.4*c.GetTopMargin())
latexf.SetTextColor(2)
f1 = TF1("f1","1.0",0.0,10.0)
f1.SetLineColor(2)
f1.SetLineWidth(2)
f1.Draw("lsame")
latexf.DrawLatex(2.5, 1.1,"68% CL")
f2 = TF1("f1","3.84",0.0,10.0)
f2.SetLineColor(2)
f2.SetLineWidth(2)
f2.Draw("lsame")
latexf.DrawLatex(2.5, 3.94,"95% CL")

gr_exp = TGraphAsymmErrors(v_r_exp,v_nll_exp,v_zeros,v_zeros,v_zeros,v_zeros)
gr_exp.SetLineColor(1)
gr_exp.SetLineWidth(2)
gr_exp.SetLineStyle(2)
gr_exp.Draw("Lsame")

#gr_obs = TGraphAsymmErrors(v_r_obs,v_nll_obs,v_zeros,v_zeros,v_zeros,v_zeros)
#gr_obs.SetLineColor(1)
#gr_obs.SetLineColor(1)
#gr_obs.SetLineWidth(2)
#gr_obs.Draw("Lsame")

latex2 = TLatex()
latex2.SetNDC()
latex2.SetTextSize(0.5*c.GetTopMargin())
latex2.SetTextFont(42)
latex2.SetTextAlign(31) # align right
latex2.DrawLatex(0.87, 0.95,"35.9 fb^{-1} (13 TeV)")
latex2.SetTextSize(0.7*c.GetTopMargin())
latex2.SetTextFont(62)
latex2.SetTextAlign(11) # align right
latex2.DrawLatex(0.20, 0.95, "CMS")
latex2.SetTextSize(0.6*c.GetTopMargin())
latex2.SetTextFont(52)
latex2.SetTextAlign(11)
latex2.DrawLatex(0.32, 0.95, "Preliminary")

legend = TLegend(.60,.14,.90,.26)
#legend.AddEntry(gr_obs , "Observed", "l")
legend.AddEntry(gr_exp , "Expected", "l")
legend.SetShadowColor(0)
legend.SetFillColor(0)
legend.SetLineColor(0)
legend.Draw("same")

gPad.RedrawAxis()

c.SaveAs("muscan.pdf")
