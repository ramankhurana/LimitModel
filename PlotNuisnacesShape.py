from OverlappingPlots import *



dirname='/afs/cern.ch/user/m/melu/public/For_Raman/2016/rtc04/'
for inuis in ['jes2016','jer2016']:
    #['pileup','prefire','eleID','dieleTrigger','jes','jer']:
    for ivar in ['ttc2016_ttZ']: #['ttc2016_ttWZ','ttc2016_ttWW','ttc2016_ttWtoLNu','ttc2016_ttZ','ttc2016_TTTo1L']:
        files=[dirname+'/TMVApp_350_ee.root']
        legend=["up", "central", "down"]
        histoname1=[ivar+'_'+inuis+'Up',ivar, ivar+'_'+inuis+'Down']
        
        xtitle='BDT discriminant'
        ytitle='# of events '
        axistitle = [xtitle, ytitle]
        DrawOverlap(files,histoname1,axistitle,legend,histoname1[1]+'_'+inuis,[0,0],[-1,1], legendheader="m_{A}=200 GeV")
