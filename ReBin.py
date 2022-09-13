'''
Algorithm: 
- get the file name
- get the histograms 
- create a new histogram, bkg_sum 
- use the fine binned histogram for this purpose, however during development, it is not availale, 
- loop over all the bins (from top to bottom). 
- if the event content is <7, sum its content with the neighbour (on left of the histogram), and increase the bin width by 1 unit
- check if the new integral is <7 or not, if <7 add the neighbouring bin on the left and increase the bin width by 1 unit
- keep repeating until the bin content is >=7 and save this bin width. 
- Use this new content and bin width to make the new shape and write it to a new file. 
- use them to get the limits, 
- repeat the whole process with 6, 7, 8, 9, 10 events and see what givess the best limits. 
'''

from ROOT import TH1F, TFile 
import copy
import os
import sys 


#processes=["TTTo1L","ttWW", "ttWZ", "ttWtoLNu", "ttZ", "ttZtoQQ", "tttW", "tttt", "tzq", "WWW", "DY", "WWZ", "WWdps", "WZ", "WZZ", "ZZZ", "osWW", "tW", "tbarW", "ttH", "ttWH", "ttWtoQQ", 
#           "ttZH", "tttJ", "zz2l", "TAToTTQ_rtcCOUPLIING_MAMASS"]

nuisnaces=["_lumiYEAR", "_pileup", "_muIDYEARsys", "_muIDYEARstat",
           "_eleIDYEARsys", "_eleIDYEARstat", "_elemuTriggerYEAR", 
           "_ctagYEARstat", "_ctagYEAREleID", "_ctagYEARLHEmuF", 
           "_ctagYEARLHEmuR", "_ctagYEARmuID", "_ctagYEARPSFSR", 
           "_ctagYEARPU", "_ctagDYXS", "_ctagSTXS", "_ctagVVXS",
           "_ctagWJetXS", "_ctagTTXS", "_ctagJER", "_ctagJES",
           "_chargeflipYEAR",
           "_sigYEARpdf",  ## only signal 
           "_metYEARunclusterE", 
           "_prefireYEAR", "_jesYEAR",
           "_jerYEAR", "_elemuTriggerYEAR", 
           "_dieleTriggerYEAR", "_dimuTriggerYEAR"]

variations=["Up", "Down"]

allvariations= [inuis+iv for iv in variations for inuis in nuisnaces]
allvariations.append("")

print ("allvariations: ", allvariations)
regions=["ee","mm","em"]
couplings=["0p1","0p4","0p8","1p0"]
masses=[200, 300, 350, 400, 500, 600, 700] 
years=["2016postapv"]#"2016", "2016apv", "2017", "2018"]
#inputdir="/eos/cms/store/group/phys_top/ExtraYukawa/BDT/BDT_output/{}/ttc_a_rtc{}_MA{}" ## YEAR, coupling, mass needs to be provided
#inputdir="/afs/cern.ch/work/k/khurana/NTU/ttc/CMSSW_10_6_29/src/ttcbar/LimitModel/BDT_output/{}/ttc_a_rtc{}_MA{}" ## YEAR, coupling, mass needs to be provided
inputdir = "BDT_output/{}/ttc_a_rtc{}_MA{}"

#outputdir="/eos/cms/store/group/phys_top/ExtraYukawa/FinalInputs"

#outputdir="/afs/cern.ch/user/z/zhenggan/work_space/CMSSW_10_2_13/src/limitmodel/FinalInputs"
outputdir="FinalInputs"

filename="TMVApp_{}_{}.root"

##filename_ = filename.format("400","ee")
##print (filename_)


signal_="TAToTTQ_rtcCOUPLIING_MAMASS"

for imass in masses: 
    for ir in regions:
        for iyear in years:
            for ic in couplings:
                
                filename_ = filename.format(str(imass), ir)
                ic_ = ic.replace("p","")
                inputdir_ = inputdir.format(iyear, ic_, str(imass)) 
                print (" fiilename: ", inputdir_+"/"+filename_)
                # print (inputdir+iyear+"/rtc"+ic.replace("p","")+"/"+filename_)
                print (inputdir_+"/"+filename_)
                rootfiilename=inputdir_+"/"+filename_
                f_in = TFile(rootfiilename,"R")
                
                f_in.cd()
                #f_in.ls()

                prefix="ttc"+iyear+"_"
                
                rebin_=5
                ## create the output file in EOS 
                print (rootfiilename.split("/")[-2], rootfiilename.split("/")[-1])
                os.system("mkdir "+outputdir+"/"+iyear+"/"+rootfiilename.split("/")[-2])
                outputfilename=outputdir+"/"+iyear+"/"+rootfiilename.split("/")[-2]+"/"+rootfiilename.split("/")[-1]
                print (outputfilename)
                fout = TFile(outputfilename,"RECREATE")

                
                ## This list needs to be altered for each year, 
                allvariations = [iv.replace("YEAR",iyear) for iv in allvariations]
                print ("systematic variations: ", allvariations)
                for inuis in allvariations:
                    
                    print (prefix+"TTTo1L"+inuis)
                    if (type( f_in.Get(prefix+"TTTo1L"+inuis) ) is TH1F) :
                        h_TTTo1L      = copy.deepcopy ( f_in.Get(prefix+"TTTo1L"+inuis) ) 
                        h_TTTo1L.Rebin(rebin_); h_TTTo1L.SetNameTitle("ttc"+iyear+"_TTTo1L"+inuis,"ttc"+iyear+"_TTTo1L"+inuis)
                        fout.cd()
                        h_TTTo1L.Write()
                        print (prefix+"TTTo1L"+inuis,"hisstogram written")
                            
                    
                    f_in.cd()
                    print (prefix+"ttW"+inuis)
                    temp_fix = ""
                    if iyear=='2018':temp_fix="toLNu"
                    elif iyear=='2016postapv':pass
                    else:raise ValueError('Please fix me.')
                    if type(f_in.Get(prefix+"ttW"+temp_fix+inuis)) is TH1F:
                        h_ttWtoLNu    = copy.deepcopy ( f_in.Get(prefix+"ttW"+temp_fix+inuis) )
                        h_ttWtoLNu.Rebin(rebin_);  h_ttWtoLNu.SetNameTitle("ttc"+iyear+"_ttWtoLNu"+inuis,"ttc"+iyear+"_ttWtoLNu"+inuis)
                        fout.cd()
                        h_ttWtoLNu.Write()
                        print (prefix+"ttWtoLNu"+inuis, "hisstogram written")
                    
                    f_in.cd()
                    if type(f_in.Get(prefix+"ttWW"+inuis)) is TH1F:
                        h_ttVV        = copy.deepcopy (f_in.Get(prefix+"ttWW"+inuis) )
                        h_ttVV.Add (f_in.Get(prefix+"ttWZ"+inuis) )
                        h_ttVV.Rebin(rebin_);  h_ttVV.SetNameTitle("ttc"+iyear+"_ttVV"+inuis,"ttc"+iyear+"_ttVV"+inuis)
                        fout.cd()
                        h_ttVV.Write()
                        
                    f_in.cd()
                    if type(f_in.Get(prefix+"ttZ"+inuis)) is TH1F:
                        h_ttV      = copy.deepcopy(f_in.Get(prefix+"ttZ"+inuis))
                        if iyear=='2018': 
                            h_ttV.Add( f_in.Get(prefix+"ttZtoQQ"+inuis) )
                            h_ttV.Add( f_in.Get(prefix+"ttWtoQQ"+inuis) )
                        elif iyear=='2016postapv': 
                            h_ttV.Add( f_in.Get(prefix+"ttZToQQ"+inuis) )
                            h_ttV.Add( f_in.Get(prefix+"ttWToQQ"+inuis) )
                        else:raise ValueError('Please fix me.')
                        h_ttV.Rebin(rebin_);  h_ttV.SetNameTitle("ttc"+iyear+"_ttV"+inuis,"ttc"+iyear+"_ttV"+inuis)
                        fout.cd()
                        h_ttV.Write()
                    
                    f_in.cd()
                    if iyear=='2018':
                        if (type(f_in.Get(prefix+"tttW"+inuis))) is TH1F:
                            h_tttX     = copy.deepcopy ( f_in.Get(prefix+"tttW"+inuis) )
                            h_tttX.Add (  f_in.Get(prefix+"tttt"+inuis) )
                            h_tttX.Rebin(rebin_);  h_tttX.SetNameTitle("ttc"+iyear+"_tttX"+inuis,"ttc"+iyear+"_tttX"+inuis)
                            fout.cd()
                            h_tttX.Write()   
                    
                    elif iyear=='2016postapv':
                        if (type(f_in.Get(prefix+"tttJ"+inuis))) is TH1F:
                            h_tttX     = copy.deepcopy ( f_in.Get(prefix+"tttJ"+inuis) )
                            h_tttX.Add (  f_in.Get(prefix+"tttt"+inuis) )
                            h_tttX.Rebin(rebin_);  h_tttX.SetNameTitle("ttc"+iyear+"_tttX"+inuis,"ttc"+iyear+"_tttX"+inuis)
                            fout.cd()
                            h_tttX.Write()   
                    else:raise ValueError('Please fix me.')


                    f_in.cd()
                    if iyear=='2018':process_name_y = 'WWW'
                    elif iyear=='2016postapv':process_name_y = 'www1'
                    else:raise ValueError('Please fix me.')

                    if (type(f_in.Get(prefix+process_name_y+inuis))) is TH1F:
                        h_VVV      = copy.deepcopy ( f_in.Get(prefix+process_name_y+inuis) )
                        
                        if iyear=='2018':
                            h_VVV.Add( f_in.Get(prefix+"WWZ"+inuis))
                            h_VVV.Add( f_in.Get(prefix+"WZZ"+inuis))
                        elif iyear=='2016postapv':
                            h_VVV.Add( f_in.Get(prefix+"wwz1"+inuis))
                            h_VVV.Add( f_in.Get(prefix+"wzz1"+inuis))

                        else:raise ValueError('Please specifiy this as well!')
                        h_VVV.Rebin(rebin_);  h_VVV.SetNameTitle("ttc"+iyear+"_VVV"+inuis,"ttc"+iyear+"_VVV"+inuis)
                        fout.cd()   
                        h_VVV.Write()
                        
                    
                    f_in.cd()
                    if (type(f_in.Get(prefix+"tW"+inuis))) is TH1F:
                        h_t_tbar_W = copy.deepcopy ( f_in.Get(prefix+"tW"+inuis) )
                        h_t_tbar_W.Add(  f_in.Get(prefix+"tbarW"+inuis) ) 
                        h_t_tbar_W.Rebin(rebin_);  h_t_tbar_W.SetNameTitle("ttc"+iyear+"_t_tbar_W"+inuis,"ttc"+iyear+"_t_tbar_W"+inuis)
                        fout.cd()
                        h_t_tbar_W.Write()
                        
                    f_in.cd()
                    if (type(f_in.Get(prefix+"ttH"+inuis))) is TH1F:
                        h_tt_V_H   = copy.deepcopy ( f_in.Get(prefix+"ttH"+inuis) )
                        h_tt_V_H.Add ( f_in.Get(prefix+"ttWH"+inuis))
                        h_tt_V_H.Add ( f_in.Get(prefix+"ttZH"+inuis))
                        h_tt_V_H.Rebin(rebin_);  h_tt_V_H.SetNameTitle("ttc"+iyear+"_tt_V_H"+inuis,"ttc"+iyear+"_tt_V_H"+inuis)
                        fout.cd()
                        h_tt_V_H.Write()
                    
                    f_in.cd()
                    if iyear=='2018':tmp_fix=''
                    elif iyear=='2016postapv':tmp_fix = 'nlo'
                    else:raise ValueError('Please fix me.')
                    
                    if (type(f_in.Get(prefix+"DY"+tmp_fix+inuis))) is TH1F:
                        h_DY       = copy.deepcopy (f_in.Get(prefix+"DY"+tmp_fix+inuis))
                        h_DY.Rebin(rebin_);  h_DY.SetNameTitle("ttc"+iyear+"_DY"+inuis,"ttc"+iyear+"_DY"+inuis)
                        fout.cd()
                        h_DY.Write()
                    
                    f_in.cd()
                    if iyear=='2018':
                        H=f_in.Get(prefix+"tzq"+inuis)
                    elif iyear=='2016postapv':
                        H=f_in.Get(prefix+"tZq"+inuis)
                    else:raise ValueError('Please fix me.')
                    if (type(H)) is TH1F:
                        if iyear=='2018':h_tzq      = copy.deepcopy (H)
                        elif iyear=='2016postapv':h_tzq      = copy.deepcopy (H)
                        else:raise ValueError('Please fix me.')

                        h_tzq      = copy.deepcopy (H) 
                        h_tzq.Rebin(rebin_);  h_tzq.SetNameTitle("ttc"+iyear+"_tzq"+inuis,"ttc"+iyear+"_tzq"+inuis)
                        fout.cd()
                        h_tzq.Write()
                        
                    f_in.cd()
                    if (type(f_in.Get(prefix+"WWdps"+inuis))) is TH1F:
                        h_VV       = copy.deepcopy (f_in.Get(prefix+"WWdps"+inuis))
                        h_VV.Add( f_in.Get(prefix+"WZ"+inuis))
                        h_VV.Add( f_in.Get(prefix+"osWW"+inuis))
                        h_VV.Rebin(rebin_);  h_VV.SetNameTitle("ttc"+iyear+"_VV"+inuis,"ttc"+iyear+"_VV"+inuis)
                        fout.cd()
                        h_VV.Write
                        
                    f_in.cd()
                    if (type(f_in.Get(prefix+"data_obs"+inuis))) is TH1F:
                        h_data_obs = copy.deepcopy(f_in.Get(prefix+"data_obs"+inuis))
                        h_data_obs.Rebin(rebin_);  h_data_obs.SetNameTitle("ttc"+iyear+"_data_obs"+inuis,"ttc"+iyear+"_data_obs"+inuis)
                        fout.cd()
                        h_data_obs.Write()
                        
                    
                    sig_name_ = prefix+(signal_.replace("MASS",str(imass))).replace("COUPLIING",ic_)
                    f_in.cd()
                    if (type(f_in.Get(sig_name_+inuis))) is TH1F:
                        h_signal_ = copy.deepcopy( f_in.Get(sig_name_+inuis)); h_signal_.Rebin(rebin_); h_signal_.SetNameTitle(sig_name_+inuis, sig_name_+inuis)
                        fout.cd()
                        h_signal_.Write()


                    
                                        

