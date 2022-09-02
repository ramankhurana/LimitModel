#!/bin/bash -e 
echo "TEST FIRST" 
echo "copy input root file"
eos cp /eos/cms/store/group/phys_top/ExtraYukawa/BDT/BDT_input_application/2018/ttc_s0_rtu01_MS0500.root .
eos cp /eos/cms/store/group/phys_top/ExtraYukawa/BDT/BDT_input_application/2018/aa.tar.gz .
tar zxf aa.tar.gz
rm aa.tar.gz
PWD=`pwd`
HOME=$PWD

# copy BDT weights, one central weight and six systematic weights
mkdir BDT_weights_0
cd BDT_weights_0
cp /eos/cms/store/group/phys_top/ExtraYukawa/BDT/BDT_weights/2018/ttc_s0_rtu01_MS0500/ttc_s0_rtu01_MS0500_0/dataset_ttc_s0_rtu01_MS0500_central/weights/TMVAClassification_BDTG.weights.xml .
cd $HOME

mkdir BDT_weights_1
cd BDT_weights_1
cp /eos/cms/store/group/phys_top/ExtraYukawa/BDT/BDT_weights/2018/ttc_s0_rtu01_MS0500/ttc_s0_rtu01_MS0500_1/dataset_ttc_s0_rtu01_MS0500_jesup/weights/TMVAClassification_BDTG.weights.xml .
cd $HOME

mkdir BDT_weights_2
cd BDT_weights_2
cp /eos/cms/store/group/phys_top/ExtraYukawa/BDT/BDT_weights/2018/ttc_s0_rtu01_MS0500/ttc_s0_rtu01_MS0500_2/dataset_ttc_s0_rtu01_MS0500_jesdo/weights/TMVAClassification_BDTG.weights.xml .
cd $HOME

mkdir BDT_weights_3
cd BDT_weights_3
cp /eos/cms/store/group/phys_top/ExtraYukawa/BDT/BDT_weights/2018/ttc_s0_rtu01_MS0500/ttc_s0_rtu01_MS0500_3/dataset_ttc_s0_rtu01_MS0500_jerup/weights/TMVAClassification_BDTG.weights.xml .
cd $HOME

mkdir BDT_weights_4
cd BDT_weights_4
cp /eos/cms/store/group/phys_top/ExtraYukawa/BDT/BDT_weights/2018/ttc_s0_rtu01_MS0500/ttc_s0_rtu01_MS0500_4/dataset_ttc_s0_rtu01_MS0500_jerdo/weights/TMVAClassification_BDTG.weights.xml .
cd $HOME

mkdir BDT_weights_5
cd BDT_weights_5
cp /eos/cms/store/group/phys_top/ExtraYukawa/BDT/BDT_weights/2018/ttc_s0_rtu01_MS0500/ttc_s0_rtu01_MS0500_5/dataset_ttc_s0_rtu01_MS0500_metUnslusEnup/weights/TMVAClassification_BDTG.weights.xml .
cd $HOME

mkdir BDT_weights_6
cd BDT_weights_6
cp /eos/cms/store/group/phys_top/ExtraYukawa/BDT/BDT_weights/2018/ttc_s0_rtu01_MS0500/ttc_s0_rtu01_MS0500_6/dataset_ttc_s0_rtu01_MS0500_metUnslusEndo/weights/TMVAClassification_BDTG.weights.xml .
cd $HOME

echo $HOME 
export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 project CMSSW CMSSW_10_6_29`
cd $PWD/CMSSW_10_6_29
ls -lrth
eval `scramv1 runtime -sh`

cd #PWD
echo "TEST DIR"

root -b -l TMVAClassificationApplication.C
echo "Run complete!!"
ls -lrth
tar zcf output.tar.gz TMVApp_*.root
ls *.root | grep -v "TMVApp_*.root" | xargs rm
rm -r BDT_weights_*
echo "end!!!"
ls -lrth
rm -rf CMSSW_10_6_29
ls -lrth
