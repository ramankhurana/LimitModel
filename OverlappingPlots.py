# In this at the end of filevector I am putting the dirname
# so loop over n-1 files and n will give the name of the output dir.
## First version by Raman Khurana: 20 June 2017
# In legend also the n element will give the name for the ratio plot y axis label.
#edited by Monika Mittal 1 March 2018 : improved the style and added possiblity for the ratio
#Script for ratio plot is not stable
import sys

import ROOT 
ROOT.gROOT.SetBatch(True)
sys.argv.append( '-b-' )


from ROOT import TFile, TH1F, gDirectory, TCanvas, TPad, TProfile,TGraph, TGraphAsymmErrors
from ROOT import TH1D, TH1, TH1I
from ROOT import gStyle
from ROOT import gROOT
from ROOT import TStyle
from ROOT import TLegend
from ROOT import TMath
from ROOT import TPaveText
from ROOT import TLatex

import os
colors=[4,1,2,3,5,9,41,46,30,12,28,20,32]
markerStyle=[23,21,22,20,24,25,26,27,28,29,20,21,22,23]            
linestyle=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]


def DrawOverlap(fileVec, histVec, titleVec,legendtext,pngname,logstatus=[0,0],xRange=[-99999,99999,1],legendheader=""):

    gStyle.SetOptTitle(0)
    gStyle.SetOptStat(0)
    gStyle.SetTitleOffset(1.1,"Y");
    #gStyle.SetTitleOffset(1.9,"X");
    gStyle.SetLineWidth(3)
    gStyle.SetFrameLineWidth(3); 

    i=0

    histList_=[]
    histList=[]
    histList1=[]
    maximum=[]
    
    ## Legend    
    
    leg = TLegend(0.1, 0.70, 0.89, 0.89)#,NULL,"brNDC");
    leg.SetHeader(legendheader)
    leg.SetBorderSize(0)
    leg.SetNColumns(2)
    leg.SetLineColor(1)
    leg.SetLineStyle(1)
    leg.SetLineWidth(1)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
     
    c = TCanvas("c1", "c1",0,0,500,500)
    c.SetBottomMargin(0.15)
    c.SetLogy(logstatus[1])
    c.SetLogx(logstatus[0])
    print ("you have provided "+str(len(fileVec))+" files and "+str(len(histVec))+" histograms to make a overlapping plot" )
    print "opening rootfiles"
    c.cd()
    
    ii=0    
    inputfile={}
    print str(fileVec[(len(fileVec)-1)])

    for ifile_ in range(len(fileVec)):
        print ("opening file  "+fileVec[ifile_])
        inputfile[ifile_] = TFile( fileVec[ifile_] )
        print "fetching histograms"
        for ihisto_ in range(len(histVec)):
            print ("printing histo "+str(histVec[ihisto_]))
            histo = inputfile[ifile_].Get(histVec[ihisto_])
            print ("integral: ", histo.Integral())
            #status_ = type(histo) is TGraphAsymmErrors
            histList.append(histo)
            # for ratio plot as they should nt be normalize 
            histList1.append(histo)
            #print histList[ii].Integral()
            #histList[ii].Rebin(xRange[2])
            if ("TAToTTQ_rtc04" in str(histVec[ihisto_])):
                histList[ii].Scale(100.0/histList[ii].Integral())
            
            maximum.append(histList[ii].GetMaximum())
            maximum.sort()
            ii=ii+1

    print histList
    print ("max:", maximum)
    print ()
    for ih in range(len(histList)):
        tt = type(histList[ih])
        if logstatus[1] is 1 :
            #histList[ih].SetMaximum(100) #1.4 for log
            #histList[ih].SetMinimum(0.1) #1.4 for log
            histList[ih].SetMaximum(maximum[0]*10)
            
            histList[ih].SetMinimum(0.005)
        if logstatus[1] is 0 :
            histList[ih].SetMaximum(maximum[-1]*1.4) #1.4 for log
            histList[ih].SetMinimum(0.005) #1.4 for log
#        print "graph_status =" ,(tt is TGraphAsymmErrors)
#        print "hist status =", (tt is TH1D) or (tt is TH1F)
         #histList[ih].Smooth()
        if ih == 0 :      
            if tt is TGraphAsymmErrors : 
                histList[ih].Draw("APL")
            if (tt is TH1D) or (tt is TH1F) or (tt is TH1) or (tt is TH1I) :
                histList[ih].Draw("HIST ")## removed hist   
        if ih > 0 :
            #histList[ih].SetLineWidth(2)
            if tt is TGraphAsymmErrors : 
                histList[ih].Draw("P L same")
            if (tt is TH1D) or (tt is TH1F) or (tt is TH1) or (tt is TH1I) :
                histList[ih].Draw("HIST   same")   ## removed hist 

        if tt is TGraphAsymmErrors :
            histList[ih].SetMaximum(100) 
            histList[ih].SetMarkerColor(colors[ih])
            histList[ih].SetLineColor(colors[ih])
            histList[ih].SetLineWidth(2)
            histList[ih].SetMarkerStyle(markerStyle[ih])
            histList[ih].SetMarkerSize(1)
            leg.AddEntry(histList[ih],legendtext[ih],"PL")
        if (tt is TH1D) or (tt is TH1F) or (tt is TH1) or (tt is TH1I) :
            histList[ih].SetLineStyle(linestyle[ih])
            histList[ih].SetLineColor(colors[ih])
            histList[ih].SetLineWidth(3)
            leg.AddEntry(histList[ih],legendtext[ih],"L")
        histList[ih].GetYaxis().SetTitle(titleVec[1])
        histList[ih].GetYaxis().SetTitleSize(0.052)
        histList[ih].GetYaxis().SetTitleOffset(0.98)
        histList[ih].GetYaxis().SetTitleFont(42)
        histList[ih].GetYaxis().SetLabelFont(42)
        histList[ih].GetYaxis().SetLabelSize(.052)
        histList[ih].GetXaxis().SetRangeUser(xRange[0],xRange[1])
        #     histList[ih].GetXaxis().SetLabelSize(0.0000);
        
        histList[ih].GetXaxis().SetTitle(titleVec[0])
        histList[ih].GetXaxis().SetLabelSize(0.052)
        histList[ih].GetXaxis().SetTitleSize(0.052)
        #histList[ih].GetXaxis().SetTitleOffset(1.14)
        histList[ih].GetXaxis().SetTitleFont(42)

        histList[ih].GetXaxis().SetLabelFont(42)
        histList[ih].GetYaxis().SetLabelFont(42) 
        histList[ih].GetXaxis().SetNdivisions(507)
        #histList[ih].GetXaxis().SetMoreLogLabels(); 
        #histList[ih].GetXaxis().SetNoExponent();
        #histList[ih].GetTGaxis().SetMaxDigits(3);

        i=i+1
    pt = TPaveText(0.01,0.92,0.95,0.96,"brNDC")
    pt.SetBorderSize(0)
    pt.SetTextAlign(12)
    pt.SetFillStyle(0)
    pt.SetTextFont(42)
    pt.SetTextSize(0.046)
    #text = pt.AddText(0.12,0.35,"CMS Internal                     36 fb^{-1} (2016) ")
    text = pt.AddText(0.12,0.35,"CMS Internal         ")#            41.5 fb^{-1} (2017) ")
    #text = pt.AddText(0.12,0.35,"CMS Internal                     59.6 fb^{-1} (2018) ")
    #text = pt.AddText(0.6,0.5,"41.5 fb^{-1} (2017) ")
    pt.Draw()
   
    

    leg.Draw()
    outputdirname = 'plots_NuisnacesShapeComparison/ee/'
    histname=outputdirname+pngname 
    c.SaveAs(histname+'.png')
    c.SaveAs(histname+'.pdf')
    outputname = 'cp  -r '+ outputdirname+'/*' +' /afs/cern.ch/work/k/khurana/public/AnalysisStuff/monoH/LimitModelPlots/plots_limit/limitcomp/'
#    os.system(outputname) 
    


print "calling the plotter"

'''

dirname='../LimitModel/inputs/'

files=[dirname+'rtc01/TMVApp_All_em.root']

legend=['200','300','350','400','500','600','700']
histoname1=['ttc2017_TAToTTQ_rtc01_MA200', 'ttc2017_TAToTTQ_rtc01_MA300','ttc2017_TAToTTQ_rtc01_MA350','ttc2017_TAToTTQ_rtc01_MA400','ttc2017_TAToTTQ_rtc01_MA500', 'ttc2017_TAToTTQ_rtc01_MA600','ttc2017_TAToTTQ_rtc01_MA700' ]

xtitle='BDT discriminant'
ytitle='# of events (normalised to 1)'
axistitle = [xtitle, ytitle]
DrawOverlap(files,histoname1,axistitle,legend,'ttc2017_TAToTTQ_rtc01_MAScan',[0,1],[-1,1],legendheader="#rho_{tc}=0.1")



files=[dirname+'rtc04/TMVApp_All_em.root']
legend=['200','300','350','400','500','600','700']
histoname1=['ttc2017_TAToTTQ_rtc01_MA200', 'ttc2017_TAToTTQ_rtc01_MA300','ttc2017_TAToTTQ_rtc01_MA350','ttc2017_TAToTTQ_rtc01_MA400','ttc2017_TAToTTQ_rtc01_MA500', 'ttc2017_TAToTTQ_rtc01_MA600','ttc2017_TAToTTQ_rtc01_MA700' ]
histoname1  =  [iname.replace("rtc01","rtc04") for iname in histoname1]
xtitle='BDT discriminant'
ytitle='# of events (normalised to 1)'
axistitle = [xtitle, ytitle]
DrawOverlap(files,histoname1,axistitle,legend,'ttc2017_TAToTTQ_rtc04_MAScan',[0,1],[-1,1], legendheader="#rho_{tc}=0.4")





files=[dirname+'rtc08/TMVApp_All_em.root']
legend=['200','300','350','400','500','600','700']
histoname1=['ttc2017_TAToTTQ_rtc01_MA200', 'ttc2017_TAToTTQ_rtc01_MA300','ttc2017_TAToTTQ_rtc01_MA350','ttc2017_TAToTTQ_rtc01_MA400','ttc2017_TAToTTQ_rtc01_MA500', 'ttc2017_TAToTTQ_rtc01_MA600','ttc2017_TAToTTQ_rtc01_MA700' ]
histoname1  =  [iname.replace("rtc01","rtc08") for iname in histoname1]
xtitle='BDT discriminant'
ytitle='# of events (normalised to 1)'
axistitle = [xtitle, ytitle]
DrawOverlap(files,histoname1,axistitle,legend,'ttc2017_TAToTTQ_rtc08_MAScan',[0,1],[-1,1], legendheader="#rho_{tc}=0.8")



files=[dirname+'rtc10/TMVApp_All_em.root']
legend=['200','300','350','400','500','600','700']
histoname1=['ttc2017_TAToTTQ_rtc01_MA200', 'ttc2017_TAToTTQ_rtc01_MA300','ttc2017_TAToTTQ_rtc01_MA350','ttc2017_TAToTTQ_rtc01_MA400','ttc2017_TAToTTQ_rtc01_MA500', 'ttc2017_TAToTTQ_rtc01_MA600','ttc2017_TAToTTQ_rtc01_MA700' ]
histoname1  =  [iname.replace("rtc01","rtc10") for iname in histoname1]
xtitle='BDT discriminant'
ytitle='# of events (normalised to 1)'
axistitle = [xtitle, ytitle]
DrawOverlap(files,histoname1,axistitle,legend,'ttc2017_TAToTTQ_rtc10_MAScan',[0,1],[-1,1], legendheader="#rho_{tc}=1.0")



files=[dirname+'rtcScan_TMVApp_All_em.root']
legend=["0.1", "0.4", "0.8", "1.0"]
histoname1=['ttc2017_TAToTTQ_rtc01_MA200', 'ttc2017_TAToTTQ_rtc04_MA200','ttc2017_TAToTTQ_rtc08_MA200','ttc2017_TAToTTQ_rtc10_MA200']
xtitle='BDT discriminant'
ytitle='# of events (normalised to 1)'
axistitle = [xtitle, ytitle]
DrawOverlap(files,histoname1,axistitle,legend,'ttc2017_TAToTTQ_rtcScan_MA200',[0,1],[-1,1], legendheader="m_{A}=200 GeV")


files=[dirname+'rtcScan_TMVApp_All_em.root']
legend=["0.1", "0.4", "0.8", "1.0"]
histoname1=['ttc2017_TAToTTQ_rtc01_MA500', 'ttc2017_TAToTTQ_rtc04_MA500','ttc2017_TAToTTQ_rtc08_MA500','ttc2017_TAToTTQ_rtc10_MA500']
xtitle='BDT discriminant'
ytitle='# of events (normalised to 1)'
axistitle = [xtitle, ytitle]
DrawOverlap(files,histoname1,axistitle,legend,'ttc2017_TAToTTQ_rtcScan_MA500',[0,1],[-1,1], legendheader="m_{A}=500 GeV")


'''
