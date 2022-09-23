from OverlappingPlots import *



dirname='/afs/cern.ch/user/m/melu/public/For_Raman/2018/rtc04/'
histoname1= ['ttc2018_TTTo1L', 'ttc2018_ttWtoLNu','ttc2018_ttZ','ttc2018_TAToTTQ_rtc04_MA350']
files=[dirname+'/TMVApp_350_mm.root']
legend=[ 'tt1L','ttWLnu','ttZ','SigmA350']
xtitle='BDT discriminant'
ytitle='# of events '
axistitle = [xtitle, ytitle]
DrawOverlap(files,histoname1,axistitle,legend,"sig_bkg_somparison_mm_2018",[0,0],[-1,1], legendheader="Yield Comparison mm") 

