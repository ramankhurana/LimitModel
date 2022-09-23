dirname=$1 #SignalExtractionChecks2017
year=$2 #2017, 2016, 2018, 20161718
catg=$3 #C, em, mm, ee 
datacard=$4 #datacards_ttc_2017/ttc_datacard_2017_SR_C_C_MA200_rtc04.txt


## source runallchecks.sh LimitChecks_19Sept2022 2017 C datacards_ttc_2017/ttc_datacard_2017_SR_C_C_MA200_rtc04.txt
 
mode=asimov_t_0

## create workspace 
text2workspace.py $datacard --channel-masks
datacardws=`echo $datacard | sed  's|.txt|.root|g'`
echo $datacardws


## run pulls and impact asimov b-only 
#### pulls 
combine -M FitDiagnostics  $datacardws --saveShapes --saveWithUncertainties -t -1 --expectSignal 0 -n _${catg}_${year}_${mode}_${dirname} --cminDefaultMinimizerStrategy 0 # --freezeParameters CMS2017_EleID #--X-rtd MINIMIZER_MaxCalls=999999999 #--setParameters mask_ZEE=1 ## --cminDefaultMinimizerPrecision 1E-12
python diffNuisances.py fitDiagnostics_${catg}_${year}_${mode}_${dirname}.root --abs --all -g pulls_${catg}_${year}_${mode}_${dirname}.root
root -l -b -q PlotPulls.C\(\"pulls_${catg}_${year}_${mode}_${dirname}.root\",\"${dirname}/\",\"_${catg}_${year}_${mode}_${dirname}\"\)




#### impacts
#--freezeParameters ratett --setParameters ratett=1.2
#text2workspace.py $datacard --channel-masks
#combineTool.py -M Impacts -d $datacardws --doInitialFit --robustFit 1 -m 125 -t -1 --expectSignal 0 --rMin -10 
#combineTool.py -M Impacts -d $datacardws --doFits  --robustFit 1 -m 125 --parallel 32 -t -1 --expectSignal 0 --rMin -10 
#combineTool.py -M Impacts -d  $datacardws -m 125 -o impacts_t0.json
#plotImpacts.py -i  impacts_t0.json -o   ${dirname}/impacts_t0_${dirname}_${catg}

## run pulls and impact asimov signal injected 
### pulls 

#mode=asimov_t_m1

#combine -M FitDiagnostics  $datacardws --saveShapes --saveWithUncertainties -t -1 --expectSignal 1 -n _${catg}_${year}_${mode}_${dirname}
#python diffNuisances.py fitDiagnostics_${catg}_${year}_${mode}_${dirname}.root --abs --all -g pulls_${catg}_${year}_${mode}_${dirname}.root
#root -l -b -q PlotPulls.C\(\"pulls_${catg}_${year}_${mode}_${dirname}.root\",\"${dirname}/\",\"_${catg}_${year}_${mode}_${dirname}\"\)


#### impacts
#combineTool.py -M Impacts -d $datacardws --doInitialFit --robustFit 1 -m 125 -t -1 --expectSignal 1 --rMin -10
#combineTool.py -M Impacts -d $datacardws --doFits  --robustFit 1 -m 125 --parallel 32 -t -1 --expectSignal 1 --rMin -10
#combineTool.py -M Impacts -d  $datacardws -m 125 -o impacts_t1.json
#plotImpacts.py -i  impacts_t1.json -o  ${dirname}/impacts_t1_${dirname}



mode=fit_CRonly_result

## CR only fit pulls 
#combine -M FitDiagnostics -d $datacardws -n _${catg}_${year}_${mode}_${dirname}  --saveShapes --saveWithUncertainties --setParameters mask_SR=1,mask_cat_1b_SR=1,mask_cat_2b_SR=1,mask_d2016_SR=1,mask_d2017_SR=1,mask_d2018_SR=1,mask_d2016_cat_1b_SR=1,mask_d2016_cat_2b_SR=1,mask_d2017_cat_1b_SR=1,mask_d2017_cat_2b_SR=1,mask_d2018_cat_1b_SR=1,mask_d2018_cat_2b_SR=1 --X-rtd MINIMIZER_analytic --cminFallbackAlgo Minuit2,0:1.0  --cminDefaultMinimizerStrategy 0
#root -l -b -q plotPostNuisance_combine.C\(\"fitDiagnostics_${catg}_${year}_${mode}_${dirname}.root\",\"${dirname}/\",\"${catg}_${year}_${mode}_${dirname}\"\)


#combine -M FitDiagnostics -d $datacardws -n _${catg}_${year}_${mode}_${dirname}  --saveShapes --saveWithUncertainties --setParameters mask_SR=1,mask_d2016_cat_1b_SR=1,mask_d2016_cat_2b_SR=1,mask_d2017_cat_1b_SR=1,mask_d2017_cat_2b_SR=1,mask_d2018_cat_1b_SR=1,mask_d2018_cat_2b_SR=1, --X-rtd MINIMIZER_analytic --cminFallbackAlgo Minuit2,0:1.0


## CR only postfit summary plot  ## need to work on this one. 
#python stackhist.py  fitDiagnostics_$catg_$year_$mode.root $catg $mode $year

#python crSummary_postfit.py -i fitDiagnostics_${catg}_${year}_${mode}_${dirname}.root -d b -c 2b -t ${catg}_${year}_${mode}_${dirname} -y 2017



#combine -M FitDiagnostics -d $datacardws -n _${catg}_${year}_${mode}_${dirname}  --saveShapes --saveWithUncertainties --freezeParameters ratett --setParameters mask_SR=1,ratett=1.2 --X-rtd MINIMIZER_analytic --cminFallbackAlgo Minuit2,0:1.0

## https://hypernews.cern.ch/HyperNews/CMS/get/higgs-combination/1621.html
## https://indico.cern.ch/event/976099/contributions/4138476/attachments/2163625/3651175/CombineTutorial-2020-debugging.pdf

##combineTool.py -M FastScan -w datacards_bbDM_2017/datacard_bbDM2017_2b_Merged_sp_0p7_tb_35_mXd_1_mA_600_ma_100.root:w

##--freezeParameters CMS2017_EleID
# --X-rescale-nuisance 'CMS2017_EleID' 0.3 
