/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides a simple example on how to use the trained classifiers
/// within an analysis module
/// - Project   : TMVA - a Root-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Executable: TMVAClassificationApplication
///
/// \macro_output
/// \macro_code
/// \author Andreas Hoecker

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;
using namespace std;

TH1F* Getoutput( TString myMethodList = "", std::string input_name="",float xs=1.0, float eff_N=1.0, std::string weight_name="", string mass_scan="", string channel="", string type_="", string cp="")
{
   cout<<"start Getoutput!!"<<endl;
   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   TH1F *histnull(0);
   float lumi=59800;
   float norm_scale=xs/eff_N;
   cout<<"norm_scale:"<<norm_scale<<", input_name:"<<input_name<<endl;
   cout<<"xs:"<<xs<<", eff_N:"<<eff_N<<endl;

   // Cut optimisation
   Use["Cuts"]            = 0;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   //
   // 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 0;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 0;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 0; // k-nearest neighbour method
   //
   // Linear Discriminant Analysis
   Use["LD"]              = 0; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // Function Discriminant analysis
   Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   Use["DNN_CPU"] = 0;         // CUDA-accelerated DNN training.
   Use["DNN_GPU"] = 0;         // Multi-core accelerated DNN.
   //
   // Support Vector Machine
   Use["SVM"]             = 0;
   //
   // Boosted Decision Trees
   Use["BDT"]             = 0; // uses Adaptive Boost
   Use["BDTG"]            = 1; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
   //
   // Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 0;
   // ---------------------------------------------------------------
   Use["Plugin"]          = 0;
   Use["Category"]        = 0;
   Use["SVM_Gauss"]       = 0;
   Use["SVM_Poly"]        = 0;
   Use["SVM_Lin"]         = 0;

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;
   std::cout << "++> Check myMethodList:"<< myMethodList << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod
                      << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
               std::cout << it->first << " ";
            }
            std::cout << std::endl;
            return histnull;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   Float_t HT, ttc_l1_pt, ttc_l2_pt, ttc_met, ttc_met_phi, ttc_mll, ttc_mllj1, ttc_mllj2, ttc_mllj3;
   Float_t dr_j1j2, dr_j1j3, dr_j2j3;
   Float_t ttc_l1_eta,ttc_l2_eta;
   Float_t j1_FlavCvB, j1_FlavCvL;
   Float_t j2_FlavCvB, j2_FlavCvL;
   Float_t j3_FlavCvB, j3_FlavCvL;
   Int_t ttc_region;
   reader->AddVariable( "HT", &HT );
   reader->AddVariable( "j1_FlavCvB", &j1_FlavCvB );
   reader->AddVariable( "j1_FlavCvL", &j1_FlavCvL );
   reader->AddVariable( "j2_FlavCvB", &j2_FlavCvB );
   reader->AddVariable( "j2_FlavCvL", &j2_FlavCvL );
   reader->AddVariable( "j3_FlavCvB", &j3_FlavCvB );
   reader->AddVariable( "j3_FlavCvL", &j3_FlavCvL );
   reader->AddVariable( "dr_j1j2", &dr_j1j2);
   reader->AddVariable( "dr_j1j3", &dr_j1j3);
   reader->AddVariable( "dr_j2j3", &dr_j2j3);
   reader->AddVariable( "ttc_l1_pt", &ttc_l1_pt );
   reader->AddVariable( "ttc_l2_pt", &ttc_l2_pt);
   reader->AddVariable( "ttc_met", &ttc_met);
   reader->AddVariable( "ttc_met_phi", &ttc_met_phi);
   reader->AddVariable( "ttc_mll", &ttc_mll);
   reader->AddVariable( "ttc_mllj1", &ttc_mllj1);
   reader->AddVariable( "ttc_mllj2", &ttc_mllj2);
   reader->AddVariable( "ttc_mllj3", &ttc_mllj3);


   // Book the MVA methods


   TString dir    = "./BDT_weights_0/";
   TString prefix = "TMVAClassification";

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = TString(it->first) + TString(" method");
         TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
         reader->BookMVA( methodName, weightfile );
      }
   }

   // Book output histograms
   UInt_t nbin = 200;
   TH1F *histLk(0);
   TH1F *histLkD(0);
   TH1F *histLkPCA(0);
   TH1F *histLkKDE(0);
   TH1F *histLkMIX(0);
   TH1F *histPD(0);
   TH1F *histPDD(0);
   TH1F *histPDPCA(0);
   TH1F *histPDEFoam(0);
   TH1F *histPDEFoamErr(0);
   TH1F *histPDEFoamSig(0);
   TH1F *histKNN(0);
   TH1F *histHm(0);
   TH1F *histFi(0);
   TH1F *histFiG(0);
   TH1F *histFiB(0);
   TH1F *histLD(0);
   TH1F *histNn(0);
   TH1F *histNnbfgs(0);
   TH1F *histNnbnn(0);
   TH1F *histNnC(0);
   TH1F *histNnT(0);
   TH1F *histBdt(0);
   TH1F *histBdtG(0);
   TH1F *histBdtB(0);
   TH1F *histBdtD(0);
   TH1F *histBdtF(0);
   TH1F *histRf(0);
   TH1F *histSVMG(0);
   TH1F *histSVMP(0);
   TH1F *histSVML(0);
   TH1F *histFDAMT(0);
   TH1F *histFDAGA(0);
   TH1F *histCat(0);
   TH1F *histPBdt(0);
   TH1F *histDnnGpu(0);
   TH1F *histDnnCpu(0);

   if (Use["Likelihood"])    histLk      = new TH1F( "MVA_Likelihood",    "MVA_Likelihood",    nbin, -1, 1 );
   if (Use["LikelihoodD"])   histLkD     = new TH1F( "MVA_LikelihoodD",   "MVA_LikelihoodD",   nbin, -1, 0.9999 );
   if (Use["LikelihoodPCA"]) histLkPCA   = new TH1F( "MVA_LikelihoodPCA", "MVA_LikelihoodPCA", nbin, -1, 1 );
   if (Use["LikelihoodKDE"]) histLkKDE   = new TH1F( "MVA_LikelihoodKDE", "MVA_LikelihoodKDE", nbin,  -0.00001, 0.99999 );
   if (Use["LikelihoodMIX"]) histLkMIX   = new TH1F( "MVA_LikelihoodMIX", "MVA_LikelihoodMIX", nbin,  0, 1 );
   if (Use["PDERS"])         histPD      = new TH1F( "MVA_PDERS",         "MVA_PDERS",         nbin,  0, 1 );
   if (Use["PDERSD"])        histPDD     = new TH1F( "MVA_PDERSD",        "MVA_PDERSD",        nbin,  0, 1 );
   if (Use["PDERSPCA"])      histPDPCA   = new TH1F( "MVA_PDERSPCA",      "MVA_PDERSPCA",      nbin,  0, 1 );
   if (Use["KNN"])           histKNN     = new TH1F( "MVA_KNN",           "MVA_KNN",           nbin,  0, 1 );
   if (Use["HMatrix"])       histHm      = new TH1F( "MVA_HMatrix",       "MVA_HMatrix",       nbin, -0.95, 1.55 );
   if (Use["Fisher"])        histFi      = new TH1F( "MVA_Fisher",        "MVA_Fisher",        nbin, -4, 4 );
   if (Use["FisherG"])       histFiG     = new TH1F( "MVA_FisherG",       "MVA_FisherG",       nbin, -1, 1 );
   if (Use["BoostedFisher"]) histFiB     = new TH1F( "MVA_BoostedFisher", "MVA_BoostedFisher", nbin, -2, 2 );
   if (Use["LD"])            histLD      = new TH1F( "MVA_LD",            "MVA_LD",            nbin, -2, 2 );
   if (Use["MLP"])           histNn      = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
   if (Use["MLPBFGS"])       histNnbfgs  = new TH1F( "MVA_MLPBFGS",       "MVA_MLPBFGS",       nbin, -1.25, 1.5 );
   if (Use["MLPBNN"])        histNnbnn   = new TH1F( "MVA_MLPBNN",        "MVA_MLPBNN",        nbin, -1.25, 1.5 );
   if (Use["CFMlpANN"])      histNnC     = new TH1F( "MVA_CFMlpANN",      "MVA_CFMlpANN",      nbin,  0, 1 );
   if (Use["TMlpANN"])       histNnT     = new TH1F( "MVA_TMlpANN",       "MVA_TMlpANN",       nbin, -1.3, 1.3 );
   if (Use["DNN_GPU"]) histDnnGpu = new TH1F("MVA_DNN_GPU", "MVA_DNN_GPU", nbin, -0.1, 1.1);
   if (Use["DNN_CPU"]) histDnnCpu = new TH1F("MVA_DNN_CPU", "MVA_DNN_CPU", nbin, -0.1, 1.1);
   if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -0.8, 0.8 );
   if (Use["BDTG"])          histBdtG    = new TH1F( input_name.c_str(),  input_name.c_str(),          nbin, -1.0, 1.0 );
   if (Use["BDTB"])          histBdtB    = new TH1F( "MVA_BDTB",          "MVA_BDTB",          nbin, -1.0, 1.0 );
   if (Use["BDTD"])          histBdtD    = new TH1F( "MVA_BDTD",          "MVA_BDTD",          nbin, -0.8, 0.8 );
   if (Use["BDTF"])          histBdtF    = new TH1F( "MVA_BDTF",          "MVA_BDTF",          nbin, -1.0, 1.0 );
   if (Use["RuleFit"])       histRf      = new TH1F( "MVA_RuleFit",       "MVA_RuleFit",       nbin, -2.0, 2.0 );
   if (Use["SVM_Gauss"])     histSVMG    = new TH1F( "MVA_SVM_Gauss",     "MVA_SVM_Gauss",     nbin,  0.0, 1.0 );
   if (Use["SVM_Poly"])      histSVMP    = new TH1F( "MVA_SVM_Poly",      "MVA_SVM_Poly",      nbin,  0.0, 1.0 );
   if (Use["SVM_Lin"])       histSVML    = new TH1F( "MVA_SVM_Lin",       "MVA_SVM_Lin",       nbin,  0.0, 1.0 );
   if (Use["FDA_MT"])        histFDAMT   = new TH1F( "MVA_FDA_MT",        "MVA_FDA_MT",        nbin, -2.0, 3.0 );
   if (Use["FDA_GA"])        histFDAGA   = new TH1F( "MVA_FDA_GA",        "MVA_FDA_GA",        nbin, -2.0, 3.0 );
   if (Use["Category"])      histCat     = new TH1F( "MVA_Category",      "MVA_Category",      nbin, -2., 2. );
   if (Use["Plugin"])        histPBdt    = new TH1F( "MVA_PBDT",          "MVA_BDT",           nbin, -0.8, 0.8 );

   // PDEFoam also returns per-event error, fill in histogram, and also fill significance
   if (Use["PDEFoam"]) {
      histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
      histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
      histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
   }

   // Book example histogram for probability (the other methods are done similarly)
   TH1F *probHistFi(0), *rarityHistFi(0);
   if (Use["Fisher"]) {
      probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
      rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
   }

   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //
   TFile *input_f(0);
   cout<<"input file:"<<input_name<<endl;
   std::string filename="./"+input_name+".root";
   input_f=TFile::Open(filename.c_str());
   
   // Event loop

   // Prepare the event tree
   // - Here the variable names have to corresponds to your tree
   // - You can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   std::cout << "--- Select signal sample" << std::endl;
   TTree* theTree = (TTree*)input_f->Get("SlimTree");

   Float_t genweight, puWeight, puWeightUp, puWeightDown, trig_SF, trig_SFup, trig_SFdo, mu_id, mu_id_sysup, mu_id_sysdo, mu_id_statup, mu_id_statdo, ele_id, ele_id_sysup, ele_id_sysdo, ele_id_statup, ele_id_statdo;
   Float_t ctag_SF, ctag_SF_statup, ctag_SF_statdo, ctag_SF_EleIDup, ctag_SF_EleIDdo, ctag_SF_LHEScaleWeightmuFup, ctag_SF_LHEScaleWeightmuFdo, ctag_SF_LHEScaleWeightmuRup, ctag_SF_LHEScaleWeightmuRdo, ctag_SF_MuIDup, ctag_SF_MuIDdo, ctag_SF_PSWeightFSRup, ctag_SF_PSWeightFSRdo, ctag_SF_PUWeightup, ctag_SF_PUWeightdo, ctag_SF_XSec_DYJetsup, ctag_SF_XSec_DYJetsdo, ctag_SF_XSec_STup, ctag_SF_XSec_STdo, ctag_SF_XSec_VVup, ctag_SF_XSec_VVdo, ctag_SF_XSec_WJetsup, ctag_SF_XSec_WJetsdo, ctag_SF_XSec_ttbarup, ctag_SF_XSec_ttbardo, ctag_SF_jerup, ctag_SF_jerdo, ctag_SF_jesTotalup, ctag_SF_jesTotaldo, charFlip_SF, charFlip_SFup, charFlip_SFdo, sig_pdfup, sig_pdfdo, sig_scaleup, sig_scaledo, sig_psup, sig_psdo;

   theTree->SetBranchAddress( "HT", &HT );
   theTree->SetBranchAddress( "ttc_region", &ttc_region );
   theTree->SetBranchAddress( "j1_FlavCvB", &j1_FlavCvB );
   theTree->SetBranchAddress( "j1_FlavCvL", &j1_FlavCvL );
   theTree->SetBranchAddress( "j2_FlavCvB", &j2_FlavCvB );
   theTree->SetBranchAddress( "j2_FlavCvL", &j2_FlavCvL );
   theTree->SetBranchAddress( "j3_FlavCvB", &j3_FlavCvB );
   theTree->SetBranchAddress( "j3_FlavCvL", &j3_FlavCvL );
   theTree->SetBranchAddress( "dr_j1j2", &dr_j1j2);
   theTree->SetBranchAddress( "dr_j1j3", &dr_j1j3);
   theTree->SetBranchAddress( "dr_j2j3", &dr_j2j3);
   theTree->SetBranchAddress( "ttc_l1_pt", &ttc_l1_pt );
   theTree->SetBranchAddress( "ttc_l2_pt", &ttc_l2_pt);
   theTree->SetBranchAddress( "ttc_l1_eta", &ttc_l1_eta );
   theTree->SetBranchAddress( "ttc_l2_eta", &ttc_l2_eta);
   theTree->SetBranchAddress( "ttc_met", &ttc_met);
   theTree->SetBranchAddress( "ttc_met_phi", &ttc_met_phi);
   theTree->SetBranchAddress( "ttc_mll", &ttc_mll);
   theTree->SetBranchAddress( "ttc_mllj1", &ttc_mllj1);
   theTree->SetBranchAddress( "ttc_mllj2", &ttc_mllj2);
   theTree->SetBranchAddress( "ttc_mllj3", &ttc_mllj3);
   theTree->SetBranchAddress( "genweight", &genweight);
   theTree->SetBranchAddress( "puWeight", &puWeight);
   theTree->SetBranchAddress( "puWeightUp", &puWeightUp);
   theTree->SetBranchAddress( "puWeightDown", &puWeightDown);
   theTree->SetBranchAddress( "trig_SF", &trig_SF);
   theTree->SetBranchAddress( "trig_SFup", &trig_SFup);
   theTree->SetBranchAddress( "trig_SFdo", &trig_SFdo);
   theTree->SetBranchAddress( "mu_id", &mu_id);
   theTree->SetBranchAddress( "mu_id_sysup", &mu_id_sysup);
   theTree->SetBranchAddress( "mu_id_sysdo", &mu_id_sysdo);
   theTree->SetBranchAddress( "mu_id_statup", &mu_id_statup);
   theTree->SetBranchAddress( "mu_id_statdo", &mu_id_statdo);
   theTree->SetBranchAddress( "ele_id", &ele_id);
   theTree->SetBranchAddress( "ele_id_sysup", &ele_id_sysup);
   theTree->SetBranchAddress( "ele_id_sysdo", &ele_id_sysdo);
   theTree->SetBranchAddress( "ele_id_statup", &ele_id_statup);
   theTree->SetBranchAddress( "ele_id_statdo", &ele_id_statdo);
   theTree->SetBranchAddress( "ctag_SF", &ctag_SF);
   theTree->SetBranchAddress( "ctag_SF_statup", &ctag_SF_statup);
   theTree->SetBranchAddress( "ctag_SF_statdo", &ctag_SF_statdo);
   theTree->SetBranchAddress( "ctag_SF_EleIDup", &ctag_SF_EleIDup);
   theTree->SetBranchAddress( "ctag_SF_EleIDdo", &ctag_SF_EleIDdo);
   theTree->SetBranchAddress( "ctag_SF_LHEScaleWeightmuFup", &ctag_SF_LHEScaleWeightmuFup);
   theTree->SetBranchAddress( "ctag_SF_LHEScaleWeightmuFdo", &ctag_SF_LHEScaleWeightmuFdo);
   theTree->SetBranchAddress( "ctag_SF_LHEScaleWeightmuRup", &ctag_SF_LHEScaleWeightmuRup);
   theTree->SetBranchAddress( "ctag_SF_LHEScaleWeightmuRdo", &ctag_SF_LHEScaleWeightmuRdo);
   theTree->SetBranchAddress( "ctag_SF_MuIDup", &ctag_SF_MuIDup);
   theTree->SetBranchAddress( "ctag_SF_MuIDdo", &ctag_SF_MuIDdo);
   theTree->SetBranchAddress( "ctag_SF_PSWeightFSRup", &ctag_SF_PSWeightFSRup);
   theTree->SetBranchAddress( "ctag_SF_PSWeightFSRdo", &ctag_SF_PSWeightFSRdo);
   theTree->SetBranchAddress( "ctag_SF_PUWeightup", &ctag_SF_PUWeightup);
   theTree->SetBranchAddress( "ctag_SF_PUWeightdo", &ctag_SF_PUWeightdo);
   theTree->SetBranchAddress( "ctag_SF_XSec_DYJetsup", &ctag_SF_XSec_DYJetsup);
   theTree->SetBranchAddress( "ctag_SF_XSec_DYJetsdo", &ctag_SF_XSec_DYJetsdo);
   theTree->SetBranchAddress( "ctag_SF_XSec_STup", &ctag_SF_XSec_STup);
   theTree->SetBranchAddress( "ctag_SF_XSec_STdo", &ctag_SF_XSec_STdo);
   theTree->SetBranchAddress( "ctag_SF_XSec_VVup", &ctag_SF_XSec_VVup);
   theTree->SetBranchAddress( "ctag_SF_XSec_VVdo", &ctag_SF_XSec_VVdo);
   theTree->SetBranchAddress( "ctag_SF_XSec_WJetsup", &ctag_SF_XSec_WJetsup);
   theTree->SetBranchAddress( "ctag_SF_XSec_WJetsdo", &ctag_SF_XSec_WJetsdo);
   theTree->SetBranchAddress( "ctag_SF_XSec_ttbarup", &ctag_SF_XSec_ttbarup);
   theTree->SetBranchAddress( "ctag_SF_XSec_ttbardo", &ctag_SF_XSec_ttbardo);
   theTree->SetBranchAddress( "ctag_SF_jerup", &ctag_SF_jerup);
   theTree->SetBranchAddress( "ctag_SF_jerdo", &ctag_SF_jerdo);
   theTree->SetBranchAddress( "ctag_SF_jesTotalup", &ctag_SF_jesTotalup);
   theTree->SetBranchAddress( "ctag_SF_jesTotaldo", &ctag_SF_jesTotaldo);
   theTree->SetBranchAddress( "charFlip_SF", &charFlip_SF);
   theTree->SetBranchAddress( "charFlip_SFup", &charFlip_SFup);
   theTree->SetBranchAddress( "charFlip_SFdo", &charFlip_SFdo);
   theTree->SetBranchAddress( "sig_pdfup", &sig_pdfup);
   theTree->SetBranchAddress( "sig_pdfdo", &sig_pdfdo);
   theTree->SetBranchAddress( "sig_scaleup", &sig_scaleup);
   theTree->SetBranchAddress( "sig_scaledo", &sig_scaledo);
   theTree->SetBranchAddress( "sig_psup", &sig_psup);
   theTree->SetBranchAddress( "sig_psdo", &sig_psdo);

   // Efficiency calculator for cut method
   Int_t    nSelCutsGA = 0;
   Double_t effS       = 0.7;

   std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();

   // store normalization for no ctag
   float ctag_norm=0.;

   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

      if (ievt%40000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      theTree->GetEntry(ievt);
      if(channel=="mm" && ttc_region!=1) continue;
      if(channel=="em" && ttc_region!=2) continue;
      if(channel=="ee" && ttc_region!=3) continue;
      if(channel=="ee" && ttc_mll>60 && ttc_mll<120) continue;
      // Return the MVA outputs and fill into histograms

      if (input_name.find("DY")!= string::npos)cout<<"MENG:"<<genweight<<" "<<norm_scale<<" "<<mu_id<<" "<<ele_id<<" "<<trig_SF<<" "<<charFlip_SF<<" "<<lumi<<endl;

      if (Use["CutsGA"]) {
         // Cuts is a special case: give the desired signal efficiency
         Bool_t passed = reader->EvaluateMVA( "CutsGA method", effS );
         if (passed) nSelCutsGA++;
      }
      //set chargeflip SF to be 1.0 for signal 
      if (input_name.find("ttc_a")!= string::npos || input_name.find("ttc_s0")!= string::npos){
	charFlip_SF=1.0;
	charFlip_SFup=1.0;
	charFlip_SFdo=1.0;
      }
      //set chargeflip SF to be 1.0 for only ee channel 
      if(channel=="em" || channel=="mm"){
	charFlip_SF=1.0;
	charFlip_SFup=1.0;
	charFlip_SFdo=1.0;
      }

      if (Use["Likelihood"   ])   histLk     ->Fill( reader->EvaluateMVA( "Likelihood method"    ) );
      if (Use["LikelihoodD"  ])   histLkD    ->Fill( reader->EvaluateMVA( "LikelihoodD method"   ) );
      if (Use["LikelihoodPCA"])   histLkPCA  ->Fill( reader->EvaluateMVA( "LikelihoodPCA method" ) );
      if (Use["LikelihoodKDE"])   histLkKDE  ->Fill( reader->EvaluateMVA( "LikelihoodKDE method" ) );
      if (Use["LikelihoodMIX"])   histLkMIX  ->Fill( reader->EvaluateMVA( "LikelihoodMIX method" ) );
      if (Use["PDERS"        ])   histPD     ->Fill( reader->EvaluateMVA( "PDERS method"         ) );
      if (Use["PDERSD"       ])   histPDD    ->Fill( reader->EvaluateMVA( "PDERSD method"        ) );
      if (Use["PDERSPCA"     ])   histPDPCA  ->Fill( reader->EvaluateMVA( "PDERSPCA method"      ) );
      if (Use["KNN"          ])   histKNN    ->Fill( reader->EvaluateMVA( "KNN method"           ) );
      if (Use["HMatrix"      ])   histHm     ->Fill( reader->EvaluateMVA( "HMatrix method"       ) );
      if (Use["Fisher"       ])   histFi     ->Fill( reader->EvaluateMVA( "Fisher method"        ) );
      if (Use["FisherG"      ])   histFiG    ->Fill( reader->EvaluateMVA( "FisherG method"       ) );
      if (Use["BoostedFisher"])   histFiB    ->Fill( reader->EvaluateMVA( "BoostedFisher method" ) );
      if (Use["LD"           ])   histLD     ->Fill( reader->EvaluateMVA( "LD method"            ) );
      if (Use["MLP"          ])   histNn     ->Fill( reader->EvaluateMVA( "MLP method"           ) );
      if (Use["MLPBFGS"      ])   histNnbfgs ->Fill( reader->EvaluateMVA( "MLPBFGS method"       ) );
      if (Use["MLPBNN"       ])   histNnbnn  ->Fill( reader->EvaluateMVA( "MLPBNN method"        ) );
      if (Use["CFMlpANN"     ])   histNnC    ->Fill( reader->EvaluateMVA( "CFMlpANN method"      ) );
      if (Use["TMlpANN"      ])   histNnT    ->Fill( reader->EvaluateMVA( "TMlpANN method"       ) );
      if (Use["DNN_GPU"]) histDnnGpu->Fill(reader->EvaluateMVA("DNN_GPU method"));
      if (Use["DNN_CPU"]) histDnnCpu->Fill(reader->EvaluateMVA("DNN_CPU method"));
      if (Use["BDT"          ])   histBdt    ->Fill( reader->EvaluateMVA( "BDT method"           ) );
      if (Use["BDTG"         ])
       { 
	 if(weight_name=="nominal_noctag"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF);
	 }
	 else if(weight_name=="central"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF;
	 }
	 else if(weight_name=="pileup_up"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*(puWeightUp/puWeight)*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*(puWeightUp/puWeight)*mu_id*ele_id*trig_SF*charFlip_SF;
	 }
	 else if(weight_name=="pileup_down"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*(puWeightDown/puWeight)*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*(puWeightDown/puWeight)*mu_id*ele_id*trig_SF*charFlip_SF;
	 }
	 else if(weight_name=="muID_sysup"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id_sysup*ele_id*trig_SF*charFlip_SF*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id_sysup*ele_id*trig_SF*charFlip_SF;
	 }
	 else if(weight_name=="muID_sysdown"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id_sysdo*ele_id*trig_SF*charFlip_SF*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id_sysdo*ele_id*trig_SF*charFlip_SF;
	 }
	 else if(weight_name=="muID_statup"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id_statup*ele_id*trig_SF*charFlip_SF*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id_statup*ele_id*trig_SF*charFlip_SF;
	 }
	 else if(weight_name=="muID_statdown"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id_statdo*ele_id*trig_SF*charFlip_SF*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id_statdo*ele_id*trig_SF*charFlip_SF;
	 }
	 else if(weight_name=="eleID_sysup"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id_sysup*trig_SF*charFlip_SF*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id*ele_id_sysup*trig_SF*charFlip_SF;
	 }
	 else if(weight_name=="eleID_sysdown"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id_sysdo*trig_SF*charFlip_SF*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id*ele_id_sysdo*trig_SF*charFlip_SF;
	 }
	 else if(weight_name=="eleID_statup"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id_statup*trig_SF*charFlip_SF*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id*ele_id_statup*trig_SF*charFlip_SF;
	 }
	 else if(weight_name=="eleID_statdown"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id_statdo*trig_SF*charFlip_SF*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id*ele_id_statdo*trig_SF*charFlip_SF;
	 }
	 else if(weight_name=="prefire_up"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF;
	 }
	 else if(weight_name=="prefire_down"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF;
	 }
	 else if(weight_name=="trigger_up"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SFup*charFlip_SF*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id*ele_id*trig_SFup*charFlip_SF;
	 }
	 else if(weight_name=="trigger_down"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SFdo*charFlip_SF*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id*ele_id*trig_SFdo*charFlip_SF;
	 }
	 else if(weight_name=="lumi_up"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*1.025*charFlip_SF*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*1.025*charFlip_SF;
	 }
	 else if(weight_name=="lumi_down"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*0.975*charFlip_SF*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*0.975*charFlip_SF;
	 }
	 else if(weight_name=="charFlip_SFup"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SFup*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SFup;
	 }
	 else if(weight_name=="charFlip_SFdo"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SFdo*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SFdo;
	 }
	 else if(weight_name=="sig_pdfup"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*sig_pdfup*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*sig_pdfup;
	 }
	 else if(weight_name=="sig_pdfdo"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*sig_pdfdo*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*sig_pdfdo;
	 }
	 else if(weight_name=="sig_scaleup"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*sig_scaleup*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*sig_scaleup;
	 }
	 else if(weight_name=="sig_scaledo"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*sig_scaledo*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*sig_scaledo;
	 }
	 else if(weight_name=="sig_psup"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*sig_psup*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*sig_psup;
	 }
	 else if(weight_name=="sig_psdo"){
	   histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*sig_psdo*ctag_SF);
	   ctag_norm+=genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*sig_psdo;
	 }
	 // ctag uncertainty
	 else if(weight_name=="ctag_statup")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_statup);
	 else if(weight_name=="ctag_statdo")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_statdo);
	 else if(weight_name=="ctag_EleIDup")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_EleIDup);
	 else if(weight_name=="ctag_EleIDdo")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_EleIDdo);
	 else if(weight_name=="ctag_LHEScaleWeightmuFup")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_LHEScaleWeightmuFup);
	 else if(weight_name=="ctag_LHEScaleWeightmuFdo")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_LHEScaleWeightmuFdo);
	 else if(weight_name=="ctag_LHEScaleWeightmuRup")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_LHEScaleWeightmuRup);
	 else if(weight_name=="ctag_LHEScaleWeightmuRdo")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_LHEScaleWeightmuRdo);
	 else if(weight_name=="ctag_MuIDup")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_MuIDup);
	 else if(weight_name=="ctag_MuIDdo")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_MuIDdo);
	 else if(weight_name=="ctag_PSWeightFSRup")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_PSWeightFSRup);
	 else if(weight_name=="ctag_PSWeightFSRdo")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_PSWeightFSRdo);
	 else if(weight_name=="ctag_PUWeightup")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_PUWeightup);
	 else if(weight_name=="ctag_PUWeightdo")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_PUWeightdo);
	 else if(weight_name=="ctag_XSec_DYJetsup")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_XSec_DYJetsup);
	 else if(weight_name=="ctag_XSec_DYJetsdo")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_XSec_DYJetsdo);
	 else if(weight_name=="ctag_XSec_STup")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_XSec_STup);
	 else if(weight_name=="ctag_XSec_STdo")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_XSec_STdo);
	 else if(weight_name=="ctag_XSec_VVup")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_XSec_VVup);
	 else if(weight_name=="ctag_XSec_VVdo")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_XSec_VVdo);
	 else if(weight_name=="ctag_XSec_WJetsup")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_XSec_WJetsup);
	 else if(weight_name=="ctag_XSec_WJetsdo")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_XSec_WJetsdo);
	 else if(weight_name=="ctag_XSec_ttbarup")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_XSec_ttbarup);
	 else if(weight_name=="ctag_XSec_ttbardo")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_XSec_ttbardo);
	 else if(weight_name=="ctag_jerup")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_jerup);
	 else if(weight_name=="ctag_jerdo")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_jerdo);
	 else if(weight_name=="ctag_jesTotalup")histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_jesTotalup);
	 else histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*charFlip_SF*ctag_SF_jesTotaldo);
	}
      if (Use["BDTB"         ])   histBdtB   ->Fill( reader->EvaluateMVA( "BDTB method"          ) );
      if (Use["BDTD"         ])   histBdtD   ->Fill( reader->EvaluateMVA( "BDTD method"          ) );
      if (Use["BDTF"         ])   histBdtF   ->Fill( reader->EvaluateMVA( "BDTF method"          ) );
      if (Use["RuleFit"      ])   histRf     ->Fill( reader->EvaluateMVA( "RuleFit method"       ) );
      if (Use["SVM_Gauss"    ])   histSVMG   ->Fill( reader->EvaluateMVA( "SVM_Gauss method"     ) );
      if (Use["SVM_Poly"     ])   histSVMP   ->Fill( reader->EvaluateMVA( "SVM_Poly method"      ) );
      if (Use["SVM_Lin"      ])   histSVML   ->Fill( reader->EvaluateMVA( "SVM_Lin method"       ) );
      if (Use["FDA_MT"       ])   histFDAMT  ->Fill( reader->EvaluateMVA( "FDA_MT method"        ) );
      if (Use["FDA_GA"       ])   histFDAGA  ->Fill( reader->EvaluateMVA( "FDA_GA method"        ) );
      if (Use["Category"     ])   histCat    ->Fill( reader->EvaluateMVA( "Category method"      ) );
      if (Use["Plugin"       ])   histPBdt   ->Fill( reader->EvaluateMVA( "P_BDT method"         ) );

      // Retrieve also per-event error
      if (Use["PDEFoam"]) {
         Double_t val = reader->EvaluateMVA( "PDEFoam method" );
         Double_t err = reader->GetMVAError();
         histPDEFoam   ->Fill( val );
         histPDEFoamErr->Fill( err );
         if (err>1.e-50) histPDEFoamSig->Fill( val/err );
      }

      // Retrieve probability instead of MVA output
      if (Use["Fisher"])   {
         probHistFi  ->Fill( reader->GetProba ( "Fisher method" ) );
         rarityHistFi->Fill( reader->GetRarity( "Fisher method" ) );
      }
   }

   //normalize histo to value without ctag
   if (!(weight_name.find("ctag")!= string::npos))histBdtG->Scale(ctag_norm/histBdtG->Integral());


   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   // Get efficiency for cuts classifier
   if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                                << " (for a required signal efficiency of " << effS << ")" << std::endl;

   if (Use["CutsGA"]) {

      // test: retrieve cuts for particular signal efficiency
      // CINT ignores dynamic_casts so we have to use a cuts-specific Reader function to acces the pointer
      TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;

      if (mcuts) {
         std::vector<Double_t> cutsMin;
         std::vector<Double_t> cutsMax;
         mcuts->GetCuts( 0.7, cutsMin, cutsMax );
         std::cout << "--- -------------------------------------------------------------" << std::endl;
         std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
         for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
            std::cout << "... Cut: "
                      << cutsMin[ivar]
                      << " < \""
                      << mcuts->GetInputVar(ivar)
                      << "\" <= "
                      << cutsMax[ivar] << std::endl;
         }
         std::cout << "--- -------------------------------------------------------------" << std::endl;
      }
   }

   delete reader;

   std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
   if (weight_name=="nominal_noctag")
    {string name_temp="ttc2018_"+input_name+"_nominal_noctag";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="central")
    {string name_temp="ttc2018_"+input_name;
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="pileup_up")
    {string name_temp="ttc2018_"+input_name+"_pileupUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="pileup_down")
    {string name_temp="ttc2018_"+input_name+"_pileupDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="prefire_up")
    {string name_temp="ttc2018_"+input_name+"_prefireUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="prefire_down")
    {string name_temp="ttc2018_"+input_name+"_prefireDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="muID_sysup")
    {string name_temp="ttc2018_"+input_name+"_muID2018sysUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="muID_sysdown")
    {string name_temp="ttc2018_"+input_name+"_muID2018sysDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="muID_statup")
    {string name_temp="ttc2018_"+input_name+"_muID2018statUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="muID_statdown")
    {string name_temp="ttc2018_"+input_name+"_muID2018statDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="eleID_sysup")
    {string name_temp="ttc2018_"+input_name+"_eleID2018sysUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="eleID_sysdown")
    {string name_temp="ttc2018_"+input_name+"_eleID2018sysDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="eleID_statup")
    {string name_temp="ttc2018_"+input_name+"_eleID2018statUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="eleID_statdown")
    {string name_temp="ttc2018_"+input_name+"_eleID2018statDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="lumi_up")
    {string name_temp="ttc2018_"+input_name+"_lumi2018Up";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="lumi_down")
    {string name_temp="ttc2018_"+input_name+"_lumi2018Down";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="trigger_up")
    {
    string name_temp="";
    if(channel=="ee"){name_temp="ttc2018_"+input_name+"_dieleTrigger2018Up";}
    if(channel=="em"){name_temp="ttc2018_"+input_name+"_elemuTrigger2018Up";}
    if(channel=="mm"){name_temp="ttc2018_"+input_name+"_dimuTrigger2018Up";}
    if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
    if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
    histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="trigger_down")
    {
    string name_temp="";
    if(channel=="ee"){name_temp="ttc2018_"+input_name+"_dieleTrigger2018Down";}
    if(channel=="em"){name_temp="ttc2018_"+input_name+"_elemuTrigger2018Down";}
    if(channel=="mm"){name_temp="ttc2018_"+input_name+"_dimuTrigger2018Down";}
    if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
    if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
    histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="charFlip_SFup")
    {string name_temp="ttc2018_"+input_name+"_chargeflip2018Up";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="charFlip_SFdo")
    {string name_temp="ttc2018_"+input_name+"_chargeflip2018Down";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="sig_pdfup")
    {string name_temp="ttc2018_"+input_name+"_sig2018pdfUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="sig_pdfdo")
    {string name_temp="ttc2018_"+input_name+"_sig2018pdfDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="sig_scaleup")
    {string name_temp="ttc2018_"+input_name+"_sig2018scaleUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="sig_scaledo")
    {string name_temp="ttc2018_"+input_name+"_sig2018scaleDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="sig_psup")
    {string name_temp="ttc2018_"+input_name+"_sig2018psUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="sig_psdo")
    {string name_temp="ttc2018_"+input_name+"_sig2018psDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   //ctag
   if (weight_name=="ctag_statup")
    {string name_temp="ttc2018_"+input_name+"_ctag2018statUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_statdo")
    {string name_temp="ttc2018_"+input_name+"_ctag2018statDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_EleIDup")
    {string name_temp="ttc2018_"+input_name+"_ctag2018EleIDUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_EleIDdo")
    {string name_temp="ttc2018_"+input_name+"_ctag2018EleIDDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_LHEScaleWeightmuFup")
    {string name_temp="ttc2018_"+input_name+"_ctag2018LHEmuFUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_LHEScaleWeightmuFdo")
    {string name_temp="ttc2018_"+input_name+"_ctag2018LHEmuFDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_LHEScaleWeightmuRup")
    {string name_temp="ttc2018_"+input_name+"_ctag2018LHEmuRUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_LHEScaleWeightmuRdo")
    {string name_temp="ttc2018_"+input_name+"_ctag2018LHEmuRDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_MuIDup")
    {string name_temp="ttc2018_"+input_name+"_ctag2018muIDUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_MuIDdo")
    {string name_temp="ttc2018_"+input_name+"_ctag2018muIDDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_PSWeightFSRup")
    {string name_temp="ttc2018_"+input_name+"_ctag2018PSFSRUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_PSWeightFSRdo")
    {string name_temp="ttc2018_"+input_name+"_ctag2018PSFSRDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_PUWeightup")
    {string name_temp="ttc2018_"+input_name+"_ctag2018PUUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_PUWeightdo")
    {string name_temp="ttc2018_"+input_name+"_ctag2018PUDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_XSec_DYJetsup")
    {string name_temp="ttc2018_"+input_name+"_ctagDYXSUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_XSec_DYJetsdo")
    {string name_temp="ttc2018_"+input_name+"_ctagDYXSDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_XSec_STup")
    {string name_temp="ttc2018_"+input_name+"_ctagSTXSUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_XSec_STdo")
    {string name_temp="ttc2018_"+input_name+"_ctagSTXSDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_XSec_VVup")
    {string name_temp="ttc2018_"+input_name+"_ctagVVXSUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_XSec_VVdo")
    {string name_temp="ttc2018_"+input_name+"_ctagVVXSDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_XSec_WJetsup")
    {string name_temp="ttc2018_"+input_name+"_ctagWJetXSUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_XSec_WJetsdo")
    {string name_temp="ttc2018_"+input_name+"_ctagWJetXSDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_XSec_ttbarup")
    {string name_temp="ttc2018_"+input_name+"_ctagTTXSUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_XSec_ttbardo")
    {string name_temp="ttc2018_"+input_name+"_ctagTTXSDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_jerup")
    {string name_temp="ttc2018_"+input_name+"_ctagJERUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_jerdo")
    {string name_temp="ttc2018_"+input_name+"_ctagJERDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_jesTotalup")
    {string name_temp="ttc2018_"+input_name+"_ctagJESUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (weight_name=="ctag_jesTotaldo")
    {string name_temp="ttc2018_"+input_name+"_ctagJESDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   input_f->Close();
   return histBdtG;
}

TH1F* Getoutput_sys( TString myMethodList = "", std::string input_name="",float xs=1.0, float eff_N=1.0, std::string system_unc="", string mass_scan="", string channel="", string type_="", string cp="")
{

   cout<<"start Getoutput_sys!!"<<endl;
   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   TH1F *histnull(0);
   float lumi=59800;
   float norm_scale=xs/eff_N;
   cout<<"norm_scale:"<<norm_scale<<", input_name:"<<input_name<<endl;
   cout<<"xs:"<<xs<<", eff_N:"<<eff_N<<endl;

   // Cut optimisation
   Use["Cuts"]            = 0;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   //
   // 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 0;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 0;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 0; // k-nearest neighbour method
   //
   // Linear Discriminant Analysis
   Use["LD"]              = 0; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // Function Discriminant analysis
   Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   Use["DNN_CPU"] = 0;         // CUDA-accelerated DNN training.
   Use["DNN_GPU"] = 0;         // Multi-core accelerated DNN.
   //
   // Support Vector Machine
   Use["SVM"]             = 0;
   //
   // Boosted Decision Trees
   Use["BDT"]             = 0; // uses Adaptive Boost
   Use["BDTG"]            = 1; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
   //
   // Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 0;
   // ---------------------------------------------------------------
   Use["Plugin"]          = 0;
   Use["Category"]        = 0;
   Use["SVM_Gauss"]       = 0;
   Use["SVM_Poly"]        = 0;
   Use["SVM_Lin"]         = 0;

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;
   std::cout << "++> Check myMethodList:"<< myMethodList << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod
                      << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
               std::cout << it->first << " ";
            }
            std::cout << std::endl;
            return histnull;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   Float_t HT, ttc_l1_pt, ttc_l2_pt, ttc_met, ttc_met_phi, ttc_mll, ttc_mllj1, ttc_mllj2, ttc_mllj3;
   Float_t dr_j1j2, dr_j1j3, dr_j2j3;
   Float_t j1_FlavCvB, j1_FlavCvL;
   Float_t j2_FlavCvB, j2_FlavCvL;
   Float_t j3_FlavCvB, j3_FlavCvL;
   Int_t ttc_region;
   

   if(system_unc=="central"){
     reader->AddVariable( "HT", &HT );
     reader->AddVariable( "j1_FlavCvB", &j1_FlavCvB );
     reader->AddVariable( "j1_FlavCvL", &j1_FlavCvL );
     reader->AddVariable( "j2_FlavCvB", &j2_FlavCvB );
     reader->AddVariable( "j2_FlavCvL", &j2_FlavCvL );
     reader->AddVariable( "j3_FlavCvB", &j3_FlavCvB );
     reader->AddVariable( "j3_FlavCvL", &j3_FlavCvL );
     reader->AddVariable( "dr_j1j2", &dr_j1j2);
     reader->AddVariable( "dr_j1j3", &dr_j1j3);
     reader->AddVariable( "dr_j2j3", &dr_j2j3);
     reader->AddVariable( "ttc_l1_pt", &ttc_l1_pt );
     reader->AddVariable( "ttc_l2_pt", &ttc_l2_pt);
     reader->AddVariable( "ttc_met", &ttc_met);
     reader->AddVariable( "ttc_met_phi", &ttc_met_phi);
     reader->AddVariable( "ttc_mll", &ttc_mll);
     reader->AddVariable( "ttc_mllj1", &ttc_mllj1);
     reader->AddVariable( "ttc_mllj2", &ttc_mllj2);
     reader->AddVariable( "ttc_mllj3", &ttc_mllj3);
   }
   if(system_unc=="jesup"){
     reader->AddVariable( "HT_jesup", &HT );
     reader->AddVariable( "j1_FlavCvB", &j1_FlavCvB );
     reader->AddVariable( "j1_FlavCvL", &j1_FlavCvL );
     reader->AddVariable( "j2_FlavCvB", &j2_FlavCvB );
     reader->AddVariable( "j2_FlavCvL", &j2_FlavCvL );
     reader->AddVariable( "j3_FlavCvB", &j3_FlavCvB );
     reader->AddVariable( "j3_FlavCvL", &j3_FlavCvL );
     reader->AddVariable( "dr_j1j2", &dr_j1j2);
     reader->AddVariable( "dr_j1j3", &dr_j1j3);
     reader->AddVariable( "dr_j2j3", &dr_j2j3);
     reader->AddVariable( "ttc_l1_pt", &ttc_l1_pt );
     reader->AddVariable( "ttc_l2_pt", &ttc_l2_pt);
     reader->AddVariable( "ttc_met_jesup", &ttc_met);
     reader->AddVariable( "ttc_met_phi_jesup", &ttc_met_phi);
     reader->AddVariable( "ttc_mll", &ttc_mll);
     reader->AddVariable( "ttc_mllj1_jesup", &ttc_mllj1);
     reader->AddVariable( "ttc_mllj2_jesup", &ttc_mllj2);
     reader->AddVariable( "ttc_mllj3_jesup", &ttc_mllj3);
   }
   if(system_unc=="jesdo"){
     reader->AddVariable( "HT_jesdo", &HT );
     reader->AddVariable( "j1_FlavCvB", &j1_FlavCvB );
     reader->AddVariable( "j1_FlavCvL", &j1_FlavCvL );
     reader->AddVariable( "j2_FlavCvB", &j2_FlavCvB );
     reader->AddVariable( "j2_FlavCvL", &j2_FlavCvL );
     reader->AddVariable( "j3_FlavCvB", &j3_FlavCvB );
     reader->AddVariable( "j3_FlavCvL", &j3_FlavCvL );
     reader->AddVariable( "dr_j1j2", &dr_j1j2);
     reader->AddVariable( "dr_j1j3", &dr_j1j3);
     reader->AddVariable( "dr_j2j3", &dr_j2j3);
     reader->AddVariable( "ttc_l1_pt", &ttc_l1_pt );
     reader->AddVariable( "ttc_l2_pt", &ttc_l2_pt);
     reader->AddVariable( "ttc_met_jesdo", &ttc_met);
     reader->AddVariable( "ttc_met_phi_jesdo", &ttc_met_phi);
     reader->AddVariable( "ttc_mll", &ttc_mll);
     reader->AddVariable( "ttc_mllj1_jesdo", &ttc_mllj1);
     reader->AddVariable( "ttc_mllj2_jesdo", &ttc_mllj2);
     reader->AddVariable( "ttc_mllj3_jesdo", &ttc_mllj3);
   }
   if(system_unc=="jerup"){
     reader->AddVariable( "HT_jerup", &HT );
     reader->AddVariable( "j1_FlavCvB", &j1_FlavCvB );
     reader->AddVariable( "j1_FlavCvL", &j1_FlavCvL );
     reader->AddVariable( "j2_FlavCvB", &j2_FlavCvB );
     reader->AddVariable( "j2_FlavCvL", &j2_FlavCvL );
     reader->AddVariable( "j3_FlavCvB", &j3_FlavCvB );
     reader->AddVariable( "j3_FlavCvL", &j3_FlavCvL );
     reader->AddVariable( "dr_j1j2", &dr_j1j2);
     reader->AddVariable( "dr_j1j3", &dr_j1j3);
     reader->AddVariable( "dr_j2j3", &dr_j2j3);
     reader->AddVariable( "ttc_l1_pt", &ttc_l1_pt );
     reader->AddVariable( "ttc_l2_pt", &ttc_l2_pt);
     reader->AddVariable( "ttc_met_jerup", &ttc_met);
     reader->AddVariable( "ttc_met_phi_jerup", &ttc_met_phi);
     reader->AddVariable( "ttc_mll", &ttc_mll);
     reader->AddVariable( "ttc_mllj1_jerup", &ttc_mllj1);
     reader->AddVariable( "ttc_mllj2_jerup", &ttc_mllj2);
     reader->AddVariable( "ttc_mllj3_jerup", &ttc_mllj3);
   }
   if(system_unc=="jerdo"){
     reader->AddVariable( "HT_jerdo", &HT );
     reader->AddVariable( "j1_FlavCvB", &j1_FlavCvB );
     reader->AddVariable( "j1_FlavCvL", &j1_FlavCvL );
     reader->AddVariable( "j2_FlavCvB", &j2_FlavCvB );
     reader->AddVariable( "j2_FlavCvL", &j2_FlavCvL );
     reader->AddVariable( "j3_FlavCvB", &j3_FlavCvB );
     reader->AddVariable( "j3_FlavCvL", &j3_FlavCvL );
     reader->AddVariable( "dr_j1j2", &dr_j1j2);
     reader->AddVariable( "dr_j1j3", &dr_j1j3);
     reader->AddVariable( "dr_j2j3", &dr_j2j3);
     reader->AddVariable( "ttc_l1_pt", &ttc_l1_pt );
     reader->AddVariable( "ttc_l2_pt", &ttc_l2_pt);
     reader->AddVariable( "ttc_met_jerdo", &ttc_met);
     reader->AddVariable( "ttc_met_phi_jerdo", &ttc_met_phi);
     reader->AddVariable( "ttc_mll", &ttc_mll);
     reader->AddVariable( "ttc_mllj1_jerdo", &ttc_mllj1);
     reader->AddVariable( "ttc_mllj2_jerdo", &ttc_mllj2);
     reader->AddVariable( "ttc_mllj3_jerdo", &ttc_mllj3);
   }
   if(system_unc=="unclusterEup"){
     reader->AddVariable( "HT", &HT );
     reader->AddVariable( "j1_FlavCvB", &j1_FlavCvB );
     reader->AddVariable( "j1_FlavCvL", &j1_FlavCvL );
     reader->AddVariable( "j2_FlavCvB", &j2_FlavCvB );
     reader->AddVariable( "j2_FlavCvL", &j2_FlavCvL );
     reader->AddVariable( "j3_FlavCvB", &j3_FlavCvB );
     reader->AddVariable( "j3_FlavCvL", &j3_FlavCvL );
     reader->AddVariable( "dr_j1j2", &dr_j1j2);
     reader->AddVariable( "dr_j1j3", &dr_j1j3);
     reader->AddVariable( "dr_j2j3", &dr_j2j3);
     reader->AddVariable( "ttc_l1_pt", &ttc_l1_pt );
     reader->AddVariable( "ttc_l2_pt", &ttc_l2_pt);
     reader->AddVariable( "ttc_met_unclusterEup", &ttc_met);
     reader->AddVariable( "ttc_met_phi_unclusterEup", &ttc_met_phi);
     reader->AddVariable( "ttc_mll", &ttc_mll);
     reader->AddVariable( "ttc_mllj1", &ttc_mllj1);
     reader->AddVariable( "ttc_mllj2", &ttc_mllj2);
     reader->AddVariable( "ttc_mllj3", &ttc_mllj3);
   }
   if(system_unc=="unclusterEdo"){
     reader->AddVariable( "HT", &HT );
     reader->AddVariable( "j1_FlavCvB", &j1_FlavCvB );
     reader->AddVariable( "j1_FlavCvL", &j1_FlavCvL );
     reader->AddVariable( "j2_FlavCvB", &j2_FlavCvB );
     reader->AddVariable( "j2_FlavCvL", &j2_FlavCvL );
     reader->AddVariable( "j3_FlavCvB", &j3_FlavCvB );
     reader->AddVariable( "j3_FlavCvL", &j3_FlavCvL );
     reader->AddVariable( "dr_j1j2", &dr_j1j2);
     reader->AddVariable( "dr_j1j3", &dr_j1j3);
     reader->AddVariable( "dr_j2j3", &dr_j2j3);
     reader->AddVariable( "ttc_l1_pt", &ttc_l1_pt );
     reader->AddVariable( "ttc_l2_pt", &ttc_l2_pt);
     reader->AddVariable( "ttc_met_unclusterEdo", &ttc_met);
     reader->AddVariable( "ttc_met_phi_unclusterEdo", &ttc_met_phi);
     reader->AddVariable( "ttc_mll", &ttc_mll);
     reader->AddVariable( "ttc_mllj1", &ttc_mllj1);
     reader->AddVariable( "ttc_mllj2", &ttc_mllj2);
     reader->AddVariable( "ttc_mllj3", &ttc_mllj3);
   }

   // Book the MVA methods

   TString dir    = "";
   if(system_unc=="central"){
     dir    = "./BDT_weights_0/";
   }
   if(system_unc=="jesup"){
     dir    = "./BDT_weights_1/";
   }
   if(system_unc=="jesdo"){
     dir    = "./BDT_weights_2/";
   }
   if(system_unc=="jerup"){
     dir    = "./BDT_weights_3/";
   }
   if(system_unc=="jerdo"){
     dir    = "./BDT_weights_4/";
   }
   if(system_unc=="unclusterEup"){
     dir    = "./BDT_weights_5/";
   }
   if(system_unc=="unclusterEdo"){
     dir    = "./BDT_weights_6/";
   }

   TString prefix = "TMVAClassification";


   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = TString(it->first) + TString(" method");
         TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
         reader->BookMVA( methodName, weightfile );
      }
   }

   // Book output histograms
   UInt_t nbin = 200;
   TH1F *histLk(0);
   TH1F *histLkD(0);
   TH1F *histLkPCA(0);
   TH1F *histLkKDE(0);
   TH1F *histLkMIX(0);
   TH1F *histPD(0);
   TH1F *histPDD(0);
   TH1F *histPDPCA(0);
   TH1F *histPDEFoam(0);
   TH1F *histPDEFoamErr(0);
   TH1F *histPDEFoamSig(0);
   TH1F *histKNN(0);
   TH1F *histHm(0);
   TH1F *histFi(0);
   TH1F *histFiG(0);
   TH1F *histFiB(0);
   TH1F *histLD(0);
   TH1F *histNn(0);
   TH1F *histNnbfgs(0);
   TH1F *histNnbnn(0);
   TH1F *histNnC(0);
   TH1F *histNnT(0);
   TH1F *histBdt(0);
   TH1F *histBdtG(0);
   TH1F *histBdtB(0);
   TH1F *histBdtD(0);
   TH1F *histBdtF(0);
   TH1F *histRf(0);
   TH1F *histSVMG(0);
   TH1F *histSVMP(0);
   TH1F *histSVML(0);
   TH1F *histFDAMT(0);
   TH1F *histFDAGA(0);
   TH1F *histCat(0);
   TH1F *histPBdt(0);
   TH1F *histDnnGpu(0);
   TH1F *histDnnCpu(0);

   if (Use["Likelihood"])    histLk      = new TH1F( "MVA_Likelihood",    "MVA_Likelihood",    nbin, -1, 1 );
   if (Use["LikelihoodD"])   histLkD     = new TH1F( "MVA_LikelihoodD",   "MVA_LikelihoodD",   nbin, -1, 0.9999 );
   if (Use["LikelihoodPCA"]) histLkPCA   = new TH1F( "MVA_LikelihoodPCA", "MVA_LikelihoodPCA", nbin, -1, 1 );
   if (Use["LikelihoodKDE"]) histLkKDE   = new TH1F( "MVA_LikelihoodKDE", "MVA_LikelihoodKDE", nbin,  -0.00001, 0.99999 );
   if (Use["LikelihoodMIX"]) histLkMIX   = new TH1F( "MVA_LikelihoodMIX", "MVA_LikelihoodMIX", nbin,  0, 1 );
   if (Use["PDERS"])         histPD      = new TH1F( "MVA_PDERS",         "MVA_PDERS",         nbin,  0, 1 );
   if (Use["PDERSD"])        histPDD     = new TH1F( "MVA_PDERSD",        "MVA_PDERSD",        nbin,  0, 1 );
   if (Use["PDERSPCA"])      histPDPCA   = new TH1F( "MVA_PDERSPCA",      "MVA_PDERSPCA",      nbin,  0, 1 );
   if (Use["KNN"])           histKNN     = new TH1F( "MVA_KNN",           "MVA_KNN",           nbin,  0, 1 );
   if (Use["HMatrix"])       histHm      = new TH1F( "MVA_HMatrix",       "MVA_HMatrix",       nbin, -0.95, 1.55 );
   if (Use["Fisher"])        histFi      = new TH1F( "MVA_Fisher",        "MVA_Fisher",        nbin, -4, 4 );
   if (Use["FisherG"])       histFiG     = new TH1F( "MVA_FisherG",       "MVA_FisherG",       nbin, -1, 1 );
   if (Use["BoostedFisher"]) histFiB     = new TH1F( "MVA_BoostedFisher", "MVA_BoostedFisher", nbin, -2, 2 );
   if (Use["LD"])            histLD      = new TH1F( "MVA_LD",            "MVA_LD",            nbin, -2, 2 );
   if (Use["MLP"])           histNn      = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
   if (Use["MLPBFGS"])       histNnbfgs  = new TH1F( "MVA_MLPBFGS",       "MVA_MLPBFGS",       nbin, -1.25, 1.5 );
   if (Use["MLPBNN"])        histNnbnn   = new TH1F( "MVA_MLPBNN",        "MVA_MLPBNN",        nbin, -1.25, 1.5 );
   if (Use["CFMlpANN"])      histNnC     = new TH1F( "MVA_CFMlpANN",      "MVA_CFMlpANN",      nbin,  0, 1 );
   if (Use["TMlpANN"])       histNnT     = new TH1F( "MVA_TMlpANN",       "MVA_TMlpANN",       nbin, -1.3, 1.3 );
   if (Use["DNN_GPU"]) histDnnGpu = new TH1F("MVA_DNN_GPU", "MVA_DNN_GPU", nbin, -0.1, 1.1);
   if (Use["DNN_CPU"]) histDnnCpu = new TH1F("MVA_DNN_CPU", "MVA_DNN_CPU", nbin, -0.1, 1.1);
   if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -0.8, 0.8 );
   if (Use["BDTG"])          histBdtG    = new TH1F( input_name.c_str(),  input_name.c_str(),          nbin, -1.0, 1.0 );
   if (Use["BDTB"])          histBdtB    = new TH1F( "MVA_BDTB",          "MVA_BDTB",          nbin, -1.0, 1.0 );
   if (Use["BDTD"])          histBdtD    = new TH1F( "MVA_BDTD",          "MVA_BDTD",          nbin, -0.8, 0.8 );
   if (Use["BDTF"])          histBdtF    = new TH1F( "MVA_BDTF",          "MVA_BDTF",          nbin, -1.0, 1.0 );
   if (Use["RuleFit"])       histRf      = new TH1F( "MVA_RuleFit",       "MVA_RuleFit",       nbin, -2.0, 2.0 );
   if (Use["SVM_Gauss"])     histSVMG    = new TH1F( "MVA_SVM_Gauss",     "MVA_SVM_Gauss",     nbin,  0.0, 1.0 );
   if (Use["SVM_Poly"])      histSVMP    = new TH1F( "MVA_SVM_Poly",      "MVA_SVM_Poly",      nbin,  0.0, 1.0 );
   if (Use["SVM_Lin"])       histSVML    = new TH1F( "MVA_SVM_Lin",       "MVA_SVM_Lin",       nbin,  0.0, 1.0 );
   if (Use["FDA_MT"])        histFDAMT   = new TH1F( "MVA_FDA_MT",        "MVA_FDA_MT",        nbin, -2.0, 3.0 );
   if (Use["FDA_GA"])        histFDAGA   = new TH1F( "MVA_FDA_GA",        "MVA_FDA_GA",        nbin, -2.0, 3.0 );
   if (Use["Category"])      histCat     = new TH1F( "MVA_Category",      "MVA_Category",      nbin, -2., 2. );
   if (Use["Plugin"])        histPBdt    = new TH1F( "MVA_PBDT",          "MVA_BDT",           nbin, -0.8, 0.8 );

   // PDEFoam also returns per-event error, fill in histogram, and also fill significance
   if (Use["PDEFoam"]) {
      histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
      histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
      histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
   }

   // Book example histogram for probability (the other methods are done similarly)
   TH1F *probHistFi(0), *rarityHistFi(0);
   if (Use["Fisher"]) {
      probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
      rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
   }

   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //
   TFile *input_f(0);
   cout<<"input file:"<<input_name<<endl;
   std::string filename="./"+input_name+".root";
   input_f=TFile::Open(filename.c_str());
   
   // Event loop

   // Prepare the event tree
   // - Here the variable names have to corresponds to your tree
   // - You can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   std::cout << "--- Select signal sample" << std::endl;
   TTree* theTree = (TTree*)input_f->Get("SlimTree");

   Float_t genweight, puWeight, trig_SF, mu_id, ele_id, charFlip_SF, ctag_SF;
   Float_t ttc_l1_eta, ttc_l2_eta;

   theTree->SetBranchAddress( "ttc_region", &ttc_region );
   theTree->SetBranchAddress( "j1_FlavCvB", &j1_FlavCvB );
   theTree->SetBranchAddress( "j1_FlavCvL", &j1_FlavCvL );
   theTree->SetBranchAddress( "j2_FlavCvB", &j2_FlavCvB );
   theTree->SetBranchAddress( "j2_FlavCvL", &j2_FlavCvL );
   theTree->SetBranchAddress( "j3_FlavCvB", &j3_FlavCvB );
   theTree->SetBranchAddress( "j3_FlavCvL", &j3_FlavCvL );
   theTree->SetBranchAddress( "dr_j1j2", &dr_j1j2);
   theTree->SetBranchAddress( "dr_j1j3", &dr_j1j3);
   theTree->SetBranchAddress( "dr_j2j3", &dr_j2j3);
   theTree->SetBranchAddress( "ttc_l1_pt", &ttc_l1_pt );
   theTree->SetBranchAddress( "ttc_l2_pt", &ttc_l2_pt);
   theTree->SetBranchAddress( "ttc_l1_eta", &ttc_l1_eta );
   theTree->SetBranchAddress( "ttc_l2_eta", &ttc_l2_eta);
   theTree->SetBranchAddress( "ttc_mll", &ttc_mll);
   theTree->SetBranchAddress( "genweight", &genweight);
   theTree->SetBranchAddress( "puWeight", &puWeight);
   theTree->SetBranchAddress( "trig_SF", &trig_SF);
   theTree->SetBranchAddress( "mu_id", &mu_id);
   theTree->SetBranchAddress( "ele_id", &ele_id);
   theTree->SetBranchAddress( "charFlip_SF", &charFlip_SF);
   theTree->SetBranchAddress( "ctag_SF", &ctag_SF);
   if(system_unc=="central"){
     theTree->SetBranchAddress( "HT", &HT );
     theTree->SetBranchAddress( "ttc_met", &ttc_met);
     theTree->SetBranchAddress( "ttc_met_phi", &ttc_met_phi);
     theTree->SetBranchAddress( "ttc_mllj1", &ttc_mllj1);
     theTree->SetBranchAddress( "ttc_mllj2", &ttc_mllj2);
     theTree->SetBranchAddress( "ttc_mllj3", &ttc_mllj3);
   }
   if(system_unc=="jesup"){
     theTree->SetBranchAddress( "HT_jesup", &HT );
     theTree->SetBranchAddress( "ttc_met_jesup", &ttc_met);
     theTree->SetBranchAddress( "ttc_met_phi_jesup", &ttc_met_phi);
     theTree->SetBranchAddress( "ttc_mllj1_jesup", &ttc_mllj1);
     theTree->SetBranchAddress( "ttc_mllj2_jesup", &ttc_mllj2);
     theTree->SetBranchAddress( "ttc_mllj3_jesup", &ttc_mllj3);
   }
   if(system_unc=="jesdo"){
     theTree->SetBranchAddress( "HT_jesdo", &HT );
     theTree->SetBranchAddress( "ttc_met_jesdo", &ttc_met);
     theTree->SetBranchAddress( "ttc_met_phi_jesdo", &ttc_met_phi);
     theTree->SetBranchAddress( "ttc_mllj1_jesdo", &ttc_mllj1);
     theTree->SetBranchAddress( "ttc_mllj2_jesdo", &ttc_mllj2);
     theTree->SetBranchAddress( "ttc_mllj3_jesdo", &ttc_mllj3);
   }
   if(system_unc=="jerup"){
     theTree->SetBranchAddress( "HT_jerup", &HT );
     theTree->SetBranchAddress( "ttc_met_jerup", &ttc_met);
     theTree->SetBranchAddress( "ttc_met_phi_jerup", &ttc_met_phi);
     theTree->SetBranchAddress( "ttc_mllj1_jerup", &ttc_mllj1);
     theTree->SetBranchAddress( "ttc_mllj2_jerup", &ttc_mllj2);
     theTree->SetBranchAddress( "ttc_mllj3_jerup", &ttc_mllj3);
   }
   if(system_unc=="jerdo"){
     theTree->SetBranchAddress( "HT_jerdo", &HT );
     theTree->SetBranchAddress( "ttc_met_jerdo", &ttc_met);
     theTree->SetBranchAddress( "ttc_met_phi_jerdo", &ttc_met_phi);
     theTree->SetBranchAddress( "ttc_mllj1_jerdo", &ttc_mllj1);
     theTree->SetBranchAddress( "ttc_mllj2_jerdo", &ttc_mllj2);
     theTree->SetBranchAddress( "ttc_mllj3_jerdo", &ttc_mllj3);
   }
   if(system_unc=="unclusterEup"){
     theTree->SetBranchAddress( "HT", &HT );
     theTree->SetBranchAddress( "ttc_met_unclusterEup", &ttc_met);
     theTree->SetBranchAddress( "ttc_met_phi_unclusterEup", &ttc_met_phi);
     theTree->SetBranchAddress( "ttc_mllj1", &ttc_mllj1);
     theTree->SetBranchAddress( "ttc_mllj2", &ttc_mllj2);
     theTree->SetBranchAddress( "ttc_mllj3", &ttc_mllj3);
   }
   if(system_unc=="unclusterEdo"){
     theTree->SetBranchAddress( "HT", &HT );
     theTree->SetBranchAddress( "ttc_met_unclusterEdo", &ttc_met);
     theTree->SetBranchAddress( "ttc_met_phi_unclusterEdo", &ttc_met_phi);
     theTree->SetBranchAddress( "ttc_mllj1", &ttc_mllj1);
     theTree->SetBranchAddress( "ttc_mllj2", &ttc_mllj2);
     theTree->SetBranchAddress( "ttc_mllj3", &ttc_mllj3);
   }
   // Efficiency calculator for cut method
   Int_t    nSelCutsGA = 0;
   Double_t effS       = 0.7;

   std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();

   // store normalization for no ctag
   float ctag_norm=0.;

   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

      if (ievt%40000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      theTree->GetEntry(ievt);
      if(channel=="mm" && ttc_region!=1) continue;
      if(channel=="em" && ttc_region!=2) continue;
      if(channel=="ee" && ttc_region!=3) continue;
      if(channel=="ee" && ttc_mll>60 && ttc_mll<120) continue;

      if(channel=="em" || channel=="mm"){
	charFlip_SF=1.0;
      }
      // Return the MVA outputs and fill into histograms

      if (Use["CutsGA"]) {
         // Cuts is a special case: give the desired signal efficiency
         Bool_t passed = reader->EvaluateMVA( "CutsGA method", effS );
         if (passed) nSelCutsGA++;
      }

      if (Use["Likelihood"   ])   histLk     ->Fill( reader->EvaluateMVA( "Likelihood method"    ) );
      if (Use["LikelihoodD"  ])   histLkD    ->Fill( reader->EvaluateMVA( "LikelihoodD method"   ) );
      if (Use["LikelihoodPCA"])   histLkPCA  ->Fill( reader->EvaluateMVA( "LikelihoodPCA method" ) );
      if (Use["LikelihoodKDE"])   histLkKDE  ->Fill( reader->EvaluateMVA( "LikelihoodKDE method" ) );
      if (Use["LikelihoodMIX"])   histLkMIX  ->Fill( reader->EvaluateMVA( "LikelihoodMIX method" ) );
      if (Use["PDERS"        ])   histPD     ->Fill( reader->EvaluateMVA( "PDERS method"         ) );
      if (Use["PDERSD"       ])   histPDD    ->Fill( reader->EvaluateMVA( "PDERSD method"        ) );
      if (Use["PDERSPCA"     ])   histPDPCA  ->Fill( reader->EvaluateMVA( "PDERSPCA method"      ) );
      if (Use["KNN"          ])   histKNN    ->Fill( reader->EvaluateMVA( "KNN method"           ) );
      if (Use["HMatrix"      ])   histHm     ->Fill( reader->EvaluateMVA( "HMatrix method"       ) );
      if (Use["Fisher"       ])   histFi     ->Fill( reader->EvaluateMVA( "Fisher method"        ) );
      if (Use["FisherG"      ])   histFiG    ->Fill( reader->EvaluateMVA( "FisherG method"       ) );
      if (Use["BoostedFisher"])   histFiB    ->Fill( reader->EvaluateMVA( "BoostedFisher method" ) );
      if (Use["LD"           ])   histLD     ->Fill( reader->EvaluateMVA( "LD method"            ) );
      if (Use["MLP"          ])   histNn     ->Fill( reader->EvaluateMVA( "MLP method"           ) );
      if (Use["MLPBFGS"      ])   histNnbfgs ->Fill( reader->EvaluateMVA( "MLPBFGS method"       ) );
      if (Use["MLPBNN"       ])   histNnbnn  ->Fill( reader->EvaluateMVA( "MLPBNN method"        ) );
      if (Use["CFMlpANN"     ])   histNnC    ->Fill( reader->EvaluateMVA( "CFMlpANN method"      ) );
      if (Use["TMlpANN"      ])   histNnT    ->Fill( reader->EvaluateMVA( "TMlpANN method"       ) );
      if (Use["DNN_GPU"]) histDnnGpu->Fill(reader->EvaluateMVA("DNN_GPU method"));
      if (Use["DNN_CPU"]) histDnnCpu->Fill(reader->EvaluateMVA("DNN_CPU method"));
      if (Use["BDT"          ])   histBdt    ->Fill( reader->EvaluateMVA( "BDT method"           ) );
      if (Use["BDTG"         ])
       { 
	 histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*ele_id*mu_id*charFlip_SF*ctag_SF*trig_SF);
	 ctag_norm+=genweight*norm_scale*lumi*ele_id*mu_id*charFlip_SF*trig_SF;
	}
      if (Use["BDTB"         ])   histBdtB   ->Fill( reader->EvaluateMVA( "BDTB method"          ) );
      if (Use["BDTD"         ])   histBdtD   ->Fill( reader->EvaluateMVA( "BDTD method"          ) );
      if (Use["BDTF"         ])   histBdtF   ->Fill( reader->EvaluateMVA( "BDTF method"          ) );
      if (Use["RuleFit"      ])   histRf     ->Fill( reader->EvaluateMVA( "RuleFit method"       ) );
      if (Use["SVM_Gauss"    ])   histSVMG   ->Fill( reader->EvaluateMVA( "SVM_Gauss method"     ) );
      if (Use["SVM_Poly"     ])   histSVMP   ->Fill( reader->EvaluateMVA( "SVM_Poly method"      ) );
      if (Use["SVM_Lin"      ])   histSVML   ->Fill( reader->EvaluateMVA( "SVM_Lin method"       ) );
      if (Use["FDA_MT"       ])   histFDAMT  ->Fill( reader->EvaluateMVA( "FDA_MT method"        ) );
      if (Use["FDA_GA"       ])   histFDAGA  ->Fill( reader->EvaluateMVA( "FDA_GA method"        ) );
      if (Use["Category"     ])   histCat    ->Fill( reader->EvaluateMVA( "Category method"      ) );
      if (Use["Plugin"       ])   histPBdt   ->Fill( reader->EvaluateMVA( "P_BDT method"         ) );

      // Retrieve also per-event error
      if (Use["PDEFoam"]) {
         Double_t val = reader->EvaluateMVA( "PDEFoam method" );
         Double_t err = reader->GetMVAError();
         histPDEFoam   ->Fill( val );
         histPDEFoamErr->Fill( err );
         if (err>1.e-50) histPDEFoamSig->Fill( val/err );
      }

      // Retrieve probability instead of MVA output
      if (Use["Fisher"])   {
         probHistFi  ->Fill( reader->GetProba ( "Fisher method" ) );
         rarityHistFi->Fill( reader->GetRarity( "Fisher method" ) );
      }
   }

   //normalize histo to value without ctag
   histBdtG->Scale(ctag_norm/histBdtG->Integral());

   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   // Get efficiency for cuts classifier
   if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                                << " (for a required signal efficiency of " << effS << ")" << std::endl;

   if (Use["CutsGA"]) {

      // test: retrieve cuts for particular signal efficiency
      // CINT ignores dynamic_casts so we have to use a cuts-specific Reader function to acces the pointer
      TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;

      if (mcuts) {
         std::vector<Double_t> cutsMin;
         std::vector<Double_t> cutsMax;
         mcuts->GetCuts( 0.7, cutsMin, cutsMax );
         std::cout << "--- -------------------------------------------------------------" << std::endl;
         std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
         for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
            std::cout << "... Cut: "
                      << cutsMin[ivar]
                      << " < \""
                      << mcuts->GetInputVar(ivar)
                      << "\" <= "
                      << cutsMax[ivar] << std::endl;
         }
         std::cout << "--- -------------------------------------------------------------" << std::endl;
      }
   }

   delete reader;

   std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
   if (system_unc=="central")
    {string name_temp="ttc2018_"+input_name;
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (system_unc=="jesup")
    {string name_temp="ttc2018_"+input_name+"_jes2018Up";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (system_unc=="jesdo")
    {string name_temp="ttc2018_"+input_name+"_jes2018Down";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (system_unc=="jerup")
    {string name_temp="ttc2018_"+input_name+"_jer2018Up";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (system_unc=="jerdo")
    {string name_temp="ttc2018_"+input_name+"_jer2018Down";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (system_unc=="unclusterEup")
    {string name_temp="ttc2018_"+input_name+"_met2018unclusterEUp";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   if (system_unc=="unclusterEdo")
    {string name_temp="ttc2018_"+input_name+"_met2018unclusterEDown";
     if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
     if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
     histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
    }
   input_f->Close();
   return histBdtG;
}

TH1F* Getoutput_fakelep( TString myMethodList = "", std::string input_name="",float xs=1.0, float eff_N=1.0, string mass_scan="", string channel="", string type_="", string cp="")
{

   cout<<"start Getoutput_fakelep!!"<<endl;
   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   TH1F *histnull(0);
   float lumi=59800;
   float norm_scale=xs/eff_N;
   cout<<"norm_scale:"<<norm_scale<<endl;

   // Cut optimisation
   Use["Cuts"]            = 0;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   //
   // 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 0;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 0;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 0; // k-nearest neighbour method
   //
   // Linear Discriminant Analysis
   Use["LD"]              = 0; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // Function Discriminant analysis
   Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   Use["DNN_CPU"] = 0;         // CUDA-accelerated DNN training.
   Use["DNN_GPU"] = 0;         // Multi-core accelerated DNN.
   //
   // Support Vector Machine
   Use["SVM"]             = 0;
   //
   // Boosted Decision Trees
   Use["BDT"]             = 0; // uses Adaptive Boost
   Use["BDTG"]            = 1; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
   //
   // Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 0;
   // ---------------------------------------------------------------
   Use["Plugin"]          = 0;
   Use["Category"]        = 0;
   Use["SVM_Gauss"]       = 0;
   Use["SVM_Poly"]        = 0;
   Use["SVM_Lin"]         = 0;

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;
   std::cout << "++> Check myMethodList:"<< myMethodList << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod
                      << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
               std::cout << it->first << " ";
            }
            std::cout << std::endl;
            return histnull;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   Float_t HT, ttc_l1_pt, ttc_l2_pt, ttc_met, ttc_met_phi, ttc_mll, ttc_mllj1, ttc_mllj2, ttc_mllj3;
   Float_t dr_j1j2, dr_j1j3, dr_j2j3;
   Float_t j1_FlavCvB, j1_FlavCvL;
   Float_t j2_FlavCvB, j2_FlavCvL;
   Float_t j3_FlavCvB, j3_FlavCvL;
   Int_t ttc_region;
   
   reader->AddVariable( "HT", &HT );
   reader->AddVariable( "j1_FlavCvB", &j1_FlavCvB );
   reader->AddVariable( "j1_FlavCvL", &j1_FlavCvL );
   reader->AddVariable( "j2_FlavCvB", &j2_FlavCvB );
   reader->AddVariable( "j2_FlavCvL", &j2_FlavCvL );
   reader->AddVariable( "j3_FlavCvB", &j3_FlavCvB );
   reader->AddVariable( "j3_FlavCvL", &j3_FlavCvL );
   reader->AddVariable( "dr_j1j2", &dr_j1j2);
   reader->AddVariable( "dr_j1j3", &dr_j1j3);
   reader->AddVariable( "dr_j2j3", &dr_j2j3);
   reader->AddVariable( "ttc_l1_pt", &ttc_l1_pt );
   reader->AddVariable( "ttc_l2_pt", &ttc_l2_pt);
   reader->AddVariable( "ttc_met", &ttc_met);
   reader->AddVariable( "ttc_met_phi", &ttc_met_phi);
   reader->AddVariable( "ttc_mll", &ttc_mll);
   reader->AddVariable( "ttc_mllj1", &ttc_mllj1);
   reader->AddVariable( "ttc_mllj2", &ttc_mllj2);
   reader->AddVariable( "ttc_mllj3", &ttc_mllj3);

   // Book the MVA methods

   TString dir    = "./BDT_weights_0/";
   TString prefix = "TMVAClassification";


   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = TString(it->first) + TString(" method");
         TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
         reader->BookMVA( methodName, weightfile );
      }
   }

   // Book output histograms
   UInt_t nbin = 200;
   TH1F *histLk(0);
   TH1F *histLkD(0);
   TH1F *histLkPCA(0);
   TH1F *histLkKDE(0);
   TH1F *histLkMIX(0);
   TH1F *histPD(0);
   TH1F *histPDD(0);
   TH1F *histPDPCA(0);
   TH1F *histPDEFoam(0);
   TH1F *histPDEFoamErr(0);
   TH1F *histPDEFoamSig(0);
   TH1F *histKNN(0);
   TH1F *histHm(0);
   TH1F *histFi(0);
   TH1F *histFiG(0);
   TH1F *histFiB(0);
   TH1F *histLD(0);
   TH1F *histNn(0);
   TH1F *histNnbfgs(0);
   TH1F *histNnbnn(0);
   TH1F *histNnC(0);
   TH1F *histNnT(0);
   TH1F *histBdt(0);
   TH1F *histBdtG(0);
   TH1F *histBdtB(0);
   TH1F *histBdtD(0);
   TH1F *histBdtF(0);
   TH1F *histRf(0);
   TH1F *histSVMG(0);
   TH1F *histSVMP(0);
   TH1F *histSVML(0);
   TH1F *histFDAMT(0);
   TH1F *histFDAGA(0);
   TH1F *histCat(0);
   TH1F *histPBdt(0);
   TH1F *histDnnGpu(0);
   TH1F *histDnnCpu(0);

   if (Use["Likelihood"])    histLk      = new TH1F( "MVA_Likelihood",    "MVA_Likelihood",    nbin, -1, 1 );
   if (Use["LikelihoodD"])   histLkD     = new TH1F( "MVA_LikelihoodD",   "MVA_LikelihoodD",   nbin, -1, 0.9999 );
   if (Use["LikelihoodPCA"]) histLkPCA   = new TH1F( "MVA_LikelihoodPCA", "MVA_LikelihoodPCA", nbin, -1, 1 );
   if (Use["LikelihoodKDE"]) histLkKDE   = new TH1F( "MVA_LikelihoodKDE", "MVA_LikelihoodKDE", nbin,  -0.00001, 0.99999 );
   if (Use["LikelihoodMIX"]) histLkMIX   = new TH1F( "MVA_LikelihoodMIX", "MVA_LikelihoodMIX", nbin,  0, 1 );
   if (Use["PDERS"])         histPD      = new TH1F( "MVA_PDERS",         "MVA_PDERS",         nbin,  0, 1 );
   if (Use["PDERSD"])        histPDD     = new TH1F( "MVA_PDERSD",        "MVA_PDERSD",        nbin,  0, 1 );
   if (Use["PDERSPCA"])      histPDPCA   = new TH1F( "MVA_PDERSPCA",      "MVA_PDERSPCA",      nbin,  0, 1 );
   if (Use["KNN"])           histKNN     = new TH1F( "MVA_KNN",           "MVA_KNN",           nbin,  0, 1 );
   if (Use["HMatrix"])       histHm      = new TH1F( "MVA_HMatrix",       "MVA_HMatrix",       nbin, -0.95, 1.55 );
   if (Use["Fisher"])        histFi      = new TH1F( "MVA_Fisher",        "MVA_Fisher",        nbin, -4, 4 );
   if (Use["FisherG"])       histFiG     = new TH1F( "MVA_FisherG",       "MVA_FisherG",       nbin, -1, 1 );
   if (Use["BoostedFisher"]) histFiB     = new TH1F( "MVA_BoostedFisher", "MVA_BoostedFisher", nbin, -2, 2 );
   if (Use["LD"])            histLD      = new TH1F( "MVA_LD",            "MVA_LD",            nbin, -2, 2 );
   if (Use["MLP"])           histNn      = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
   if (Use["MLPBFGS"])       histNnbfgs  = new TH1F( "MVA_MLPBFGS",       "MVA_MLPBFGS",       nbin, -1.25, 1.5 );
   if (Use["MLPBNN"])        histNnbnn   = new TH1F( "MVA_MLPBNN",        "MVA_MLPBNN",        nbin, -1.25, 1.5 );
   if (Use["CFMlpANN"])      histNnC     = new TH1F( "MVA_CFMlpANN",      "MVA_CFMlpANN",      nbin,  0, 1 );
   if (Use["TMlpANN"])       histNnT     = new TH1F( "MVA_TMlpANN",       "MVA_TMlpANN",       nbin, -1.3, 1.3 );
   if (Use["DNN_GPU"]) histDnnGpu = new TH1F("MVA_DNN_GPU", "MVA_DNN_GPU", nbin, -0.1, 1.1);
   if (Use["DNN_CPU"]) histDnnCpu = new TH1F("MVA_DNN_CPU", "MVA_DNN_CPU", nbin, -0.1, 1.1);
   if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -0.8, 0.8 );
   if (Use["BDTG"])          histBdtG    = new TH1F( input_name.c_str(),  input_name.c_str(),          nbin, -1.0, 1.0 );
   if (Use["BDTB"])          histBdtB    = new TH1F( "MVA_BDTB",          "MVA_BDTB",          nbin, -1.0, 1.0 );
   if (Use["BDTD"])          histBdtD    = new TH1F( "MVA_BDTD",          "MVA_BDTD",          nbin, -0.8, 0.8 );
   if (Use["BDTF"])          histBdtF    = new TH1F( "MVA_BDTF",          "MVA_BDTF",          nbin, -1.0, 1.0 );
   if (Use["RuleFit"])       histRf      = new TH1F( "MVA_RuleFit",       "MVA_RuleFit",       nbin, -2.0, 2.0 );
   if (Use["SVM_Gauss"])     histSVMG    = new TH1F( "MVA_SVM_Gauss",     "MVA_SVM_Gauss",     nbin,  0.0, 1.0 );
   if (Use["SVM_Poly"])      histSVMP    = new TH1F( "MVA_SVM_Poly",      "MVA_SVM_Poly",      nbin,  0.0, 1.0 );
   if (Use["SVM_Lin"])       histSVML    = new TH1F( "MVA_SVM_Lin",       "MVA_SVM_Lin",       nbin,  0.0, 1.0 );
   if (Use["FDA_MT"])        histFDAMT   = new TH1F( "MVA_FDA_MT",        "MVA_FDA_MT",        nbin, -2.0, 3.0 );
   if (Use["FDA_GA"])        histFDAGA   = new TH1F( "MVA_FDA_GA",        "MVA_FDA_GA",        nbin, -2.0, 3.0 );
   if (Use["Category"])      histCat     = new TH1F( "MVA_Category",      "MVA_Category",      nbin, -2., 2. );
   if (Use["Plugin"])        histPBdt    = new TH1F( "MVA_PBDT",          "MVA_BDT",           nbin, -0.8, 0.8 );

   // PDEFoam also returns per-event error, fill in histogram, and also fill significance
   if (Use["PDEFoam"]) {
      histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
      histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
      histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
   }

   // Book example histogram for probability (the other methods are done similarly)
   TH1F *probHistFi(0), *rarityHistFi(0);
   if (Use["Fisher"]) {
      probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
      rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
   }

   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //
   TFile *input_f(0);
   cout<<"input file:"<<input_name<<endl;
   std::string filename="./"+input_name+".root";
   input_f=TFile::Open(filename.c_str());
   
   // Event loop

   // Prepare the event tree
   // - Here the variable names have to corresponds to your tree
   // - You can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   std::cout << "--- Select signal sample" << std::endl;
   TTree* theTree = (TTree*)input_f->Get("SlimTree");

   Float_t ttc_l1_eta, ttc_l2_eta, fakeweight;

   theTree->SetBranchAddress( "ttc_region", &ttc_region );
   theTree->SetBranchAddress( "j1_FlavCvB", &j1_FlavCvB );
   theTree->SetBranchAddress( "j1_FlavCvL", &j1_FlavCvL );
   theTree->SetBranchAddress( "j2_FlavCvB", &j2_FlavCvB );
   theTree->SetBranchAddress( "j2_FlavCvL", &j2_FlavCvL );
   theTree->SetBranchAddress( "j3_FlavCvB", &j3_FlavCvB );
   theTree->SetBranchAddress( "j3_FlavCvL", &j3_FlavCvL );
   theTree->SetBranchAddress( "dr_j1j2", &dr_j1j2);
   theTree->SetBranchAddress( "dr_j1j3", &dr_j1j3);
   theTree->SetBranchAddress( "dr_j2j3", &dr_j2j3);
   theTree->SetBranchAddress( "ttc_l1_pt", &ttc_l1_pt );
   theTree->SetBranchAddress( "ttc_l2_pt", &ttc_l2_pt);
   theTree->SetBranchAddress( "ttc_l1_eta", &ttc_l1_eta );
   theTree->SetBranchAddress( "ttc_l2_eta", &ttc_l2_eta);
   theTree->SetBranchAddress( "ttc_mll", &ttc_mll);
   theTree->SetBranchAddress( "HT", &HT );
   theTree->SetBranchAddress( "ttc_met", &ttc_met);
   theTree->SetBranchAddress( "ttc_met_phi", &ttc_met_phi);
   theTree->SetBranchAddress( "ttc_mllj1", &ttc_mllj1);
   theTree->SetBranchAddress( "ttc_mllj2", &ttc_mllj2);
   theTree->SetBranchAddress( "ttc_mllj3", &ttc_mllj3);
   theTree->SetBranchAddress( "fakeweight", &fakeweight);

   // Efficiency calculator for cut method
   Int_t    nSelCutsGA = 0;
   Double_t effS       = 0.7;

   std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

      if (ievt%40000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      theTree->GetEntry(ievt);
      if(channel=="mm" && ttc_region!=1) continue;
      if(channel=="em" && ttc_region!=2) continue;
      if(channel=="ee" && ttc_region!=3) continue;
      if(channel=="ee" && ttc_mll>60 && ttc_mll<120) continue;
      // Return the MVA outputs and fill into histograms

      if (Use["CutsGA"]) {
         // Cuts is a special case: give the desired signal efficiency
         Bool_t passed = reader->EvaluateMVA( "CutsGA method", effS );
         if (passed) nSelCutsGA++;
      }

      if (Use["Likelihood"   ])   histLk     ->Fill( reader->EvaluateMVA( "Likelihood method"    ) );
      if (Use["LikelihoodD"  ])   histLkD    ->Fill( reader->EvaluateMVA( "LikelihoodD method"   ) );
      if (Use["LikelihoodPCA"])   histLkPCA  ->Fill( reader->EvaluateMVA( "LikelihoodPCA method" ) );
      if (Use["LikelihoodKDE"])   histLkKDE  ->Fill( reader->EvaluateMVA( "LikelihoodKDE method" ) );
      if (Use["LikelihoodMIX"])   histLkMIX  ->Fill( reader->EvaluateMVA( "LikelihoodMIX method" ) );
      if (Use["PDERS"        ])   histPD     ->Fill( reader->EvaluateMVA( "PDERS method"         ) );
      if (Use["PDERSD"       ])   histPDD    ->Fill( reader->EvaluateMVA( "PDERSD method"        ) );
      if (Use["PDERSPCA"     ])   histPDPCA  ->Fill( reader->EvaluateMVA( "PDERSPCA method"      ) );
      if (Use["KNN"          ])   histKNN    ->Fill( reader->EvaluateMVA( "KNN method"           ) );
      if (Use["HMatrix"      ])   histHm     ->Fill( reader->EvaluateMVA( "HMatrix method"       ) );
      if (Use["Fisher"       ])   histFi     ->Fill( reader->EvaluateMVA( "Fisher method"        ) );
      if (Use["FisherG"      ])   histFiG    ->Fill( reader->EvaluateMVA( "FisherG method"       ) );
      if (Use["BoostedFisher"])   histFiB    ->Fill( reader->EvaluateMVA( "BoostedFisher method" ) );
      if (Use["LD"           ])   histLD     ->Fill( reader->EvaluateMVA( "LD method"            ) );
      if (Use["MLP"          ])   histNn     ->Fill( reader->EvaluateMVA( "MLP method"           ) );
      if (Use["MLPBFGS"      ])   histNnbfgs ->Fill( reader->EvaluateMVA( "MLPBFGS method"       ) );
      if (Use["MLPBNN"       ])   histNnbnn  ->Fill( reader->EvaluateMVA( "MLPBNN method"        ) );
      if (Use["CFMlpANN"     ])   histNnC    ->Fill( reader->EvaluateMVA( "CFMlpANN method"      ) );
      if (Use["TMlpANN"      ])   histNnT    ->Fill( reader->EvaluateMVA( "TMlpANN method"       ) );
      if (Use["DNN_GPU"]) histDnnGpu->Fill(reader->EvaluateMVA("DNN_GPU method"));
      if (Use["DNN_CPU"]) histDnnCpu->Fill(reader->EvaluateMVA("DNN_CPU method"));
      if (Use["BDT"          ])   histBdt    ->Fill( reader->EvaluateMVA( "BDT method"           ) );
      if (Use["BDTG"         ])
       { 
	 histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"), fakeweight);
	}
      if (Use["BDTB"         ])   histBdtB   ->Fill( reader->EvaluateMVA( "BDTB method"          ) );
      if (Use["BDTD"         ])   histBdtD   ->Fill( reader->EvaluateMVA( "BDTD method"          ) );
      if (Use["BDTF"         ])   histBdtF   ->Fill( reader->EvaluateMVA( "BDTF method"          ) );
      if (Use["RuleFit"      ])   histRf     ->Fill( reader->EvaluateMVA( "RuleFit method"       ) );
      if (Use["SVM_Gauss"    ])   histSVMG   ->Fill( reader->EvaluateMVA( "SVM_Gauss method"     ) );
      if (Use["SVM_Poly"     ])   histSVMP   ->Fill( reader->EvaluateMVA( "SVM_Poly method"      ) );
      if (Use["SVM_Lin"      ])   histSVML   ->Fill( reader->EvaluateMVA( "SVM_Lin method"       ) );
      if (Use["FDA_MT"       ])   histFDAMT  ->Fill( reader->EvaluateMVA( "FDA_MT method"        ) );
      if (Use["FDA_GA"       ])   histFDAGA  ->Fill( reader->EvaluateMVA( "FDA_GA method"        ) );
      if (Use["Category"     ])   histCat    ->Fill( reader->EvaluateMVA( "Category method"      ) );
      if (Use["Plugin"       ])   histPBdt   ->Fill( reader->EvaluateMVA( "P_BDT method"         ) );

      // Retrieve also per-event error
      if (Use["PDEFoam"]) {
         Double_t val = reader->EvaluateMVA( "PDEFoam method" );
         Double_t err = reader->GetMVAError();
         histPDEFoam   ->Fill( val );
         histPDEFoamErr->Fill( err );
         if (err>1.e-50) histPDEFoamSig->Fill( val/err );
      }

      // Retrieve probability instead of MVA output
      if (Use["Fisher"])   {
         probHistFi  ->Fill( reader->GetProba ( "Fisher method" ) );
         rarityHistFi->Fill( reader->GetRarity( "Fisher method" ) );
      }
   }

   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   // Get efficiency for cuts classifier
   if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                                << " (for a required signal efficiency of " << effS << ")" << std::endl;

   if (Use["CutsGA"]) {

      // test: retrieve cuts for particular signal efficiency
      // CINT ignores dynamic_casts so we have to use a cuts-specific Reader function to acces the pointer
      TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;

      if (mcuts) {
         std::vector<Double_t> cutsMin;
         std::vector<Double_t> cutsMax;
         mcuts->GetCuts( 0.7, cutsMin, cutsMax );
         std::cout << "--- -------------------------------------------------------------" << std::endl;
         std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
         for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
            std::cout << "... Cut: "
                      << cutsMin[ivar]
                      << " < \""
                      << mcuts->GetInputVar(ivar)
                      << "\" <= "
                      << cutsMax[ivar] << std::endl;
         }
         std::cout << "--- -------------------------------------------------------------" << std::endl;
      }
   }

   // Write histograms

   delete reader;

   std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
   string name_temp="ttc2018_"+input_name;
   if (input_name.find("ttc_a")!= string::npos) name_temp.replace(8,5,"TAToTTQ");
   if (input_name.find("ttc_s0")!= string::npos) name_temp.replace(8,6,"TS0ToTTQ");
   histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());
   
   input_f->Close();
   return histBdtG;
}

TH1F* Getoutput_fake_mc( TString myMethodList = "", std::string input_name="",float xs=1.0, float eff_N=1.0, std::string weight_name="", string mass_scan="", string channel="", string type_="", string cp="")
{
   cout<<"start Getoutput_fake_mc!!"<<endl;
   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   TH1F *histnull(0);
   float lumi=59800;
   float norm_scale=xs/eff_N;
   cout<<"norm_scale:"<<norm_scale<<", input_name:"<<input_name<<endl;
   cout<<"xs:"<<xs<<", eff_N:"<<eff_N<<endl;

   // Cut optimisation
   Use["Cuts"]            = 0;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   //
   // 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 0;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 0;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 0; // k-nearest neighbour method
   //
   // Linear Discriminant Analysis
   Use["LD"]              = 0; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // Function Discriminant analysis
   Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   Use["DNN_CPU"] = 0;         // CUDA-accelerated DNN training.
   Use["DNN_GPU"] = 0;         // Multi-core accelerated DNN.
   //
   // Support Vector Machine
   Use["SVM"]             = 0;
   //
   // Boosted Decision Trees
   Use["BDT"]             = 0; // uses Adaptive Boost
   Use["BDTG"]            = 1; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
   //
   // Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 0;
   // ---------------------------------------------------------------
   Use["Plugin"]          = 0;
   Use["Category"]        = 0;
   Use["SVM_Gauss"]       = 0;
   Use["SVM_Poly"]        = 0;
   Use["SVM_Lin"]         = 0;

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;
   std::cout << "++> Check myMethodList:"<< myMethodList << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod
                      << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
               std::cout << it->first << " ";
            }
            std::cout << std::endl;
            return histnull;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   Float_t HT, ttc_l1_pt, ttc_l2_pt, ttc_met, ttc_met_phi, ttc_mll, ttc_mllj1, ttc_mllj2, ttc_mllj3;
   Float_t dr_j1j2, dr_j1j3, dr_j2j3;
   Float_t ttc_l1_eta,ttc_l2_eta;
   Float_t j1_FlavCvB, j1_FlavCvL;
   Float_t j2_FlavCvB, j2_FlavCvL;
   Float_t j3_FlavCvB, j3_FlavCvL;
   Int_t ttc_region;
   reader->AddVariable( "HT", &HT );
   reader->AddVariable( "j1_FlavCvB", &j1_FlavCvB );
   reader->AddVariable( "j1_FlavCvL", &j1_FlavCvL );
   reader->AddVariable( "j2_FlavCvB", &j2_FlavCvB );
   reader->AddVariable( "j2_FlavCvL", &j2_FlavCvL );
   reader->AddVariable( "j3_FlavCvB", &j3_FlavCvB );
   reader->AddVariable( "j3_FlavCvL", &j3_FlavCvL );
   reader->AddVariable( "dr_j1j2", &dr_j1j2);
   reader->AddVariable( "dr_j1j3", &dr_j1j3);
   reader->AddVariable( "dr_j2j3", &dr_j2j3);
   reader->AddVariable( "ttc_l1_pt", &ttc_l1_pt );
   reader->AddVariable( "ttc_l2_pt", &ttc_l2_pt);
   reader->AddVariable( "ttc_met", &ttc_met);
   reader->AddVariable( "ttc_met_phi", &ttc_met_phi);
   reader->AddVariable( "ttc_mll", &ttc_mll);
   reader->AddVariable( "ttc_mllj1", &ttc_mllj1);
   reader->AddVariable( "ttc_mllj2", &ttc_mllj2);
   reader->AddVariable( "ttc_mllj3", &ttc_mllj3);

   // Book the MVA methods


   TString dir    = "./BDT_weights_0/";
   TString prefix = "TMVAClassification";

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = TString(it->first) + TString(" method");
         TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
         reader->BookMVA( methodName, weightfile );
      }
   }

   // Book output histograms
   UInt_t nbin = 200;
   TH1F *histLk(0);
   TH1F *histLkD(0);
   TH1F *histLkPCA(0);
   TH1F *histLkKDE(0);
   TH1F *histLkMIX(0);
   TH1F *histPD(0);
   TH1F *histPDD(0);
   TH1F *histPDPCA(0);
   TH1F *histPDEFoam(0);
   TH1F *histPDEFoamErr(0);
   TH1F *histPDEFoamSig(0);
   TH1F *histKNN(0);
   TH1F *histHm(0);
   TH1F *histFi(0);
   TH1F *histFiG(0);
   TH1F *histFiB(0);
   TH1F *histLD(0);
   TH1F *histNn(0);
   TH1F *histNnbfgs(0);
   TH1F *histNnbnn(0);
   TH1F *histNnC(0);
   TH1F *histNnT(0);
   TH1F *histBdt(0);
   TH1F *histBdtG(0);
   TH1F *histBdtB(0);
   TH1F *histBdtD(0);
   TH1F *histBdtF(0);
   TH1F *histRf(0);
   TH1F *histSVMG(0);
   TH1F *histSVMP(0);
   TH1F *histSVML(0);
   TH1F *histFDAMT(0);
   TH1F *histFDAGA(0);
   TH1F *histCat(0);
   TH1F *histPBdt(0);
   TH1F *histDnnGpu(0);
   TH1F *histDnnCpu(0);

   if (Use["Likelihood"])    histLk      = new TH1F( "MVA_Likelihood",    "MVA_Likelihood",    nbin, -1, 1 );
   if (Use["LikelihoodD"])   histLkD     = new TH1F( "MVA_LikelihoodD",   "MVA_LikelihoodD",   nbin, -1, 0.9999 );
   if (Use["LikelihoodPCA"]) histLkPCA   = new TH1F( "MVA_LikelihoodPCA", "MVA_LikelihoodPCA", nbin, -1, 1 );
   if (Use["LikelihoodKDE"]) histLkKDE   = new TH1F( "MVA_LikelihoodKDE", "MVA_LikelihoodKDE", nbin,  -0.00001, 0.99999 );
   if (Use["LikelihoodMIX"]) histLkMIX   = new TH1F( "MVA_LikelihoodMIX", "MVA_LikelihoodMIX", nbin,  0, 1 );
   if (Use["PDERS"])         histPD      = new TH1F( "MVA_PDERS",         "MVA_PDERS",         nbin,  0, 1 );
   if (Use["PDERSD"])        histPDD     = new TH1F( "MVA_PDERSD",        "MVA_PDERSD",        nbin,  0, 1 );
   if (Use["PDERSPCA"])      histPDPCA   = new TH1F( "MVA_PDERSPCA",      "MVA_PDERSPCA",      nbin,  0, 1 );
   if (Use["KNN"])           histKNN     = new TH1F( "MVA_KNN",           "MVA_KNN",           nbin,  0, 1 );
   if (Use["HMatrix"])       histHm      = new TH1F( "MVA_HMatrix",       "MVA_HMatrix",       nbin, -0.95, 1.55 );
   if (Use["Fisher"])        histFi      = new TH1F( "MVA_Fisher",        "MVA_Fisher",        nbin, -4, 4 );
   if (Use["FisherG"])       histFiG     = new TH1F( "MVA_FisherG",       "MVA_FisherG",       nbin, -1, 1 );
   if (Use["BoostedFisher"]) histFiB     = new TH1F( "MVA_BoostedFisher", "MVA_BoostedFisher", nbin, -2, 2 );
   if (Use["LD"])            histLD      = new TH1F( "MVA_LD",            "MVA_LD",            nbin, -2, 2 );
   if (Use["MLP"])           histNn      = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
   if (Use["MLPBFGS"])       histNnbfgs  = new TH1F( "MVA_MLPBFGS",       "MVA_MLPBFGS",       nbin, -1.25, 1.5 );
   if (Use["MLPBNN"])        histNnbnn   = new TH1F( "MVA_MLPBNN",        "MVA_MLPBNN",        nbin, -1.25, 1.5 );
   if (Use["CFMlpANN"])      histNnC     = new TH1F( "MVA_CFMlpANN",      "MVA_CFMlpANN",      nbin,  0, 1 );
   if (Use["TMlpANN"])       histNnT     = new TH1F( "MVA_TMlpANN",       "MVA_TMlpANN",       nbin, -1.3, 1.3 );
   if (Use["DNN_GPU"]) histDnnGpu = new TH1F("MVA_DNN_GPU", "MVA_DNN_GPU", nbin, -0.1, 1.1);
   if (Use["DNN_CPU"]) histDnnCpu = new TH1F("MVA_DNN_CPU", "MVA_DNN_CPU", nbin, -0.1, 1.1);
   if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -0.8, 0.8 );
   if (Use["BDTG"])          histBdtG    = new TH1F( input_name.c_str(),  input_name.c_str(),          nbin, -1.0, 1.0 );
   if (Use["BDTB"])          histBdtB    = new TH1F( "MVA_BDTB",          "MVA_BDTB",          nbin, -1.0, 1.0 );
   if (Use["BDTD"])          histBdtD    = new TH1F( "MVA_BDTD",          "MVA_BDTD",          nbin, -0.8, 0.8 );
   if (Use["BDTF"])          histBdtF    = new TH1F( "MVA_BDTF",          "MVA_BDTF",          nbin, -1.0, 1.0 );
   if (Use["RuleFit"])       histRf      = new TH1F( "MVA_RuleFit",       "MVA_RuleFit",       nbin, -2.0, 2.0 );
   if (Use["SVM_Gauss"])     histSVMG    = new TH1F( "MVA_SVM_Gauss",     "MVA_SVM_Gauss",     nbin,  0.0, 1.0 );
   if (Use["SVM_Poly"])      histSVMP    = new TH1F( "MVA_SVM_Poly",      "MVA_SVM_Poly",      nbin,  0.0, 1.0 );
   if (Use["SVM_Lin"])       histSVML    = new TH1F( "MVA_SVM_Lin",       "MVA_SVM_Lin",       nbin,  0.0, 1.0 );
   if (Use["FDA_MT"])        histFDAMT   = new TH1F( "MVA_FDA_MT",        "MVA_FDA_MT",        nbin, -2.0, 3.0 );
   if (Use["FDA_GA"])        histFDAGA   = new TH1F( "MVA_FDA_GA",        "MVA_FDA_GA",        nbin, -2.0, 3.0 );
   if (Use["Category"])      histCat     = new TH1F( "MVA_Category",      "MVA_Category",      nbin, -2., 2. );
   if (Use["Plugin"])        histPBdt    = new TH1F( "MVA_PBDT",          "MVA_BDT",           nbin, -0.8, 0.8 );

   // PDEFoam also returns per-event error, fill in histogram, and also fill significance
   if (Use["PDEFoam"]) {
      histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
      histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
      histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
   }

   // Book example histogram for probability (the other methods are done similarly)
   TH1F *probHistFi(0), *rarityHistFi(0);
   if (Use["Fisher"]) {
      probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
      rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
   }

   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //
   TFile *input_f(0);
   cout<<"input file:"<<input_name<<endl;
   std::string filename="./"+input_name+".root";
   input_f=TFile::Open(filename.c_str());
   
   // Event loop

   // Prepare the event tree
   // - Here the variable names have to corresponds to your tree
   // - You can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   std::cout << "--- Select signal sample" << std::endl;
   TTree* theTree = (TTree*)input_f->Get("SlimTree");
   if(theTree->GetEntries()==0) return histBdtG;

   Float_t genweight, puWeight, trig_SF, mu_id, ele_id;
   Float_t ctag_SF, fakeweight;

   theTree->SetBranchAddress( "HT", &HT );
   theTree->SetBranchAddress( "ttc_region", &ttc_region );
   theTree->SetBranchAddress( "j1_FlavCvB", &j1_FlavCvB );
   theTree->SetBranchAddress( "j1_FlavCvL", &j1_FlavCvL );
   theTree->SetBranchAddress( "j2_FlavCvB", &j2_FlavCvB );
   theTree->SetBranchAddress( "j2_FlavCvL", &j2_FlavCvL );
   theTree->SetBranchAddress( "j3_FlavCvB", &j3_FlavCvB );
   theTree->SetBranchAddress( "j3_FlavCvL", &j3_FlavCvL );
   theTree->SetBranchAddress( "dr_j1j2", &dr_j1j2);
   theTree->SetBranchAddress( "dr_j1j3", &dr_j1j3);
   theTree->SetBranchAddress( "dr_j2j3", &dr_j2j3);
   theTree->SetBranchAddress( "ttc_l1_pt", &ttc_l1_pt );
   theTree->SetBranchAddress( "ttc_l2_pt", &ttc_l2_pt);
   theTree->SetBranchAddress( "ttc_l1_eta", &ttc_l1_eta );
   theTree->SetBranchAddress( "ttc_l2_eta", &ttc_l2_eta);
   theTree->SetBranchAddress( "ttc_met", &ttc_met);
   theTree->SetBranchAddress( "ttc_met_phi", &ttc_met_phi);
   theTree->SetBranchAddress( "ttc_mll", &ttc_mll);
   theTree->SetBranchAddress( "ttc_mllj1", &ttc_mllj1);
   theTree->SetBranchAddress( "ttc_mllj2", &ttc_mllj2);
   theTree->SetBranchAddress( "ttc_mllj3", &ttc_mllj3);
   theTree->SetBranchAddress( "genweight", &genweight);
   theTree->SetBranchAddress( "puWeight", &puWeight);
   theTree->SetBranchAddress( "trig_SF", &trig_SF);
   theTree->SetBranchAddress( "mu_id", &mu_id);
   theTree->SetBranchAddress( "ele_id", &ele_id);
   theTree->SetBranchAddress( "ctag_SF", &ctag_SF);
   theTree->SetBranchAddress( "fakeweight", &fakeweight);

   // Efficiency calculator for cut method
   Int_t    nSelCutsGA = 0;
   Double_t effS       = 0.7;

   std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();

   // store normalization for no ctag
   float ctag_norm=0.;

   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

      if (ievt%40000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      theTree->GetEntry(ievt);
      if(channel=="mm" && ttc_region!=1) continue;
      if(channel=="em" && ttc_region!=2) continue;
      if(channel=="ee" && ttc_region!=3) continue;
      if(channel=="ee" && ttc_mll>60 && ttc_mll<120) continue;
      // Return the MVA outputs and fill into histograms

      if (Use["CutsGA"]) {
         // Cuts is a special case: give the desired signal efficiency
         Bool_t passed = reader->EvaluateMVA( "CutsGA method", effS );
         if (passed) nSelCutsGA++;
      }

      if (Use["Likelihood"   ])   histLk     ->Fill( reader->EvaluateMVA( "Likelihood method"    ) );
      if (Use["LikelihoodD"  ])   histLkD    ->Fill( reader->EvaluateMVA( "LikelihoodD method"   ) );
      if (Use["LikelihoodPCA"])   histLkPCA  ->Fill( reader->EvaluateMVA( "LikelihoodPCA method" ) );
      if (Use["LikelihoodKDE"])   histLkKDE  ->Fill( reader->EvaluateMVA( "LikelihoodKDE method" ) );
      if (Use["LikelihoodMIX"])   histLkMIX  ->Fill( reader->EvaluateMVA( "LikelihoodMIX method" ) );
      if (Use["PDERS"        ])   histPD     ->Fill( reader->EvaluateMVA( "PDERS method"         ) );
      if (Use["PDERSD"       ])   histPDD    ->Fill( reader->EvaluateMVA( "PDERSD method"        ) );
      if (Use["PDERSPCA"     ])   histPDPCA  ->Fill( reader->EvaluateMVA( "PDERSPCA method"      ) );
      if (Use["KNN"          ])   histKNN    ->Fill( reader->EvaluateMVA( "KNN method"           ) );
      if (Use["HMatrix"      ])   histHm     ->Fill( reader->EvaluateMVA( "HMatrix method"       ) );
      if (Use["Fisher"       ])   histFi     ->Fill( reader->EvaluateMVA( "Fisher method"        ) );
      if (Use["FisherG"      ])   histFiG    ->Fill( reader->EvaluateMVA( "FisherG method"       ) );
      if (Use["BoostedFisher"])   histFiB    ->Fill( reader->EvaluateMVA( "BoostedFisher method" ) );
      if (Use["LD"           ])   histLD     ->Fill( reader->EvaluateMVA( "LD method"            ) );
      if (Use["MLP"          ])   histNn     ->Fill( reader->EvaluateMVA( "MLP method"           ) );
      if (Use["MLPBFGS"      ])   histNnbfgs ->Fill( reader->EvaluateMVA( "MLPBFGS method"       ) );
      if (Use["MLPBNN"       ])   histNnbnn  ->Fill( reader->EvaluateMVA( "MLPBNN method"        ) );
      if (Use["CFMlpANN"     ])   histNnC    ->Fill( reader->EvaluateMVA( "CFMlpANN method"      ) );
      if (Use["TMlpANN"      ])   histNnT    ->Fill( reader->EvaluateMVA( "TMlpANN method"       ) );
      if (Use["DNN_GPU"]) histDnnGpu->Fill(reader->EvaluateMVA("DNN_GPU method"));
      if (Use["DNN_CPU"]) histDnnCpu->Fill(reader->EvaluateMVA("DNN_CPU method"));
      if (Use["BDT"          ])   histBdt    ->Fill( reader->EvaluateMVA( "BDT method"           ) );
      if (Use["BDTG"         ])
       { 
	histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*ctag_SF*fakeweight);
	ctag_norm+=genweight*norm_scale*lumi*mu_id*ele_id*trig_SF*fakeweight;
	}
      if (Use["BDTB"         ])   histBdtB   ->Fill( reader->EvaluateMVA( "BDTB method"          ) );
      if (Use["BDTD"         ])   histBdtD   ->Fill( reader->EvaluateMVA( "BDTD method"          ) );
      if (Use["BDTF"         ])   histBdtF   ->Fill( reader->EvaluateMVA( "BDTF method"          ) );
      if (Use["RuleFit"      ])   histRf     ->Fill( reader->EvaluateMVA( "RuleFit method"       ) );
      if (Use["SVM_Gauss"    ])   histSVMG   ->Fill( reader->EvaluateMVA( "SVM_Gauss method"     ) );
      if (Use["SVM_Poly"     ])   histSVMP   ->Fill( reader->EvaluateMVA( "SVM_Poly method"      ) );
      if (Use["SVM_Lin"      ])   histSVML   ->Fill( reader->EvaluateMVA( "SVM_Lin method"       ) );
      if (Use["FDA_MT"       ])   histFDAMT  ->Fill( reader->EvaluateMVA( "FDA_MT method"        ) );
      if (Use["FDA_GA"       ])   histFDAGA  ->Fill( reader->EvaluateMVA( "FDA_GA method"        ) );
      if (Use["Category"     ])   histCat    ->Fill( reader->EvaluateMVA( "Category method"      ) );
      if (Use["Plugin"       ])   histPBdt   ->Fill( reader->EvaluateMVA( "P_BDT method"         ) );

      // Retrieve also per-event error
      if (Use["PDEFoam"]) {
         Double_t val = reader->EvaluateMVA( "PDEFoam method" );
         Double_t err = reader->GetMVAError();
         histPDEFoam   ->Fill( val );
         histPDEFoamErr->Fill( err );
         if (err>1.e-50) histPDEFoamSig->Fill( val/err );
      }

      // Retrieve probability instead of MVA output
      if (Use["Fisher"])   {
         probHistFi  ->Fill( reader->GetProba ( "Fisher method" ) );
         rarityHistFi->Fill( reader->GetRarity( "Fisher method" ) );
      }
   }

   //normalize histo to value without ctag
   histBdtG->Scale(ctag_norm/histBdtG->Integral());


   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   // Get efficiency for cuts classifier
   if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                                << " (for a required signal efficiency of " << effS << ")" << std::endl;

   if (Use["CutsGA"]) {

      // test: retrieve cuts for particular signal efficiency
      // CINT ignores dynamic_casts so we have to use a cuts-specific Reader function to acces the pointer
      TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;

      if (mcuts) {
         std::vector<Double_t> cutsMin;
         std::vector<Double_t> cutsMax;
         mcuts->GetCuts( 0.7, cutsMin, cutsMax );
         std::cout << "--- -------------------------------------------------------------" << std::endl;
         std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
         for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
            std::cout << "... Cut: "
                      << cutsMin[ivar]
                      << " < \""
                      << mcuts->GetInputVar(ivar)
                      << "\" <= "
                      << cutsMax[ivar] << std::endl;
         }
         std::cout << "--- -------------------------------------------------------------" << std::endl;
      }
   }

   delete reader;

   std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
   string name_temp="ttc2018_"+input_name+"_fake";
   histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());

   input_f->Close();
   return histBdtG;
}

TH1F* Getoutput_data( TString myMethodList = "", std::string input_name="",float xs=1.0, float eff_N=1.0, string mass_scan="", string channel="", string type_="", string cp="")
{
   cout<<"start Getoutput_data!!"<<endl;
   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   TH1F *histnull(0);
   float lumi=59800;
   float norm_scale=xs/eff_N;
   cout<<"norm_scale:"<<norm_scale<<", input_name:"<<input_name<<endl;
   cout<<"xs:"<<xs<<", eff_N:"<<eff_N<<endl;

   // Cut optimisation
   Use["Cuts"]            = 0;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   //
   // 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 0;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 0;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 0; // k-nearest neighbour method
   //
   // Linear Discriminant Analysis
   Use["LD"]              = 0; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // Function Discriminant analysis
   Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   Use["DNN_CPU"] = 0;         // CUDA-accelerated DNN training.
   Use["DNN_GPU"] = 0;         // Multi-core accelerated DNN.
   //
   // Support Vector Machine
   Use["SVM"]             = 0;
   //
   // Boosted Decision Trees
   Use["BDT"]             = 0; // uses Adaptive Boost
   Use["BDTG"]            = 1; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
   //
   // Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 0;
   // ---------------------------------------------------------------
   Use["Plugin"]          = 0;
   Use["Category"]        = 0;
   Use["SVM_Gauss"]       = 0;
   Use["SVM_Poly"]        = 0;
   Use["SVM_Lin"]         = 0;

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;
   std::cout << "++> Check myMethodList:"<< myMethodList << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod
                      << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
               std::cout << it->first << " ";
            }
            std::cout << std::endl;
            return histnull;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   Float_t HT, ttc_l1_pt, ttc_l2_pt, ttc_met, ttc_met_phi, ttc_mll, ttc_mllj1, ttc_mllj2, ttc_mllj3;
   Float_t dr_j1j2, dr_j1j3, dr_j2j3;
   Float_t ttc_l1_eta,ttc_l2_eta;
   Float_t j1_FlavCvB, j1_FlavCvL;
   Float_t j2_FlavCvB, j2_FlavCvL;
   Float_t j3_FlavCvB, j3_FlavCvL;
   Int_t ttc_region;
   reader->AddVariable( "HT", &HT );
   reader->AddVariable( "j1_FlavCvB", &j1_FlavCvB );
   reader->AddVariable( "j1_FlavCvL", &j1_FlavCvL );
   reader->AddVariable( "j2_FlavCvB", &j2_FlavCvB );
   reader->AddVariable( "j2_FlavCvL", &j2_FlavCvL );
   reader->AddVariable( "j3_FlavCvB", &j3_FlavCvB );
   reader->AddVariable( "j3_FlavCvL", &j3_FlavCvL );
   reader->AddVariable( "dr_j1j2", &dr_j1j2);
   reader->AddVariable( "dr_j1j3", &dr_j1j3);
   reader->AddVariable( "dr_j2j3", &dr_j2j3);
   reader->AddVariable( "ttc_l1_pt", &ttc_l1_pt );
   reader->AddVariable( "ttc_l2_pt", &ttc_l2_pt);
   reader->AddVariable( "ttc_met", &ttc_met);
   reader->AddVariable( "ttc_met_phi", &ttc_met_phi);
   reader->AddVariable( "ttc_mll", &ttc_mll);
   reader->AddVariable( "ttc_mllj1", &ttc_mllj1);
   reader->AddVariable( "ttc_mllj2", &ttc_mllj2);
   reader->AddVariable( "ttc_mllj3", &ttc_mllj3);

   // Book the MVA methods


   TString dir    = "./BDT_weights_0/";
   TString prefix = "TMVAClassification";

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = TString(it->first) + TString(" method");
         TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
         reader->BookMVA( methodName, weightfile );
      }
   }

   // Book output histograms
   UInt_t nbin = 200;
   TH1F *histLk(0);
   TH1F *histLkD(0);
   TH1F *histLkPCA(0);
   TH1F *histLkKDE(0);
   TH1F *histLkMIX(0);
   TH1F *histPD(0);
   TH1F *histPDD(0);
   TH1F *histPDPCA(0);
   TH1F *histPDEFoam(0);
   TH1F *histPDEFoamErr(0);
   TH1F *histPDEFoamSig(0);
   TH1F *histKNN(0);
   TH1F *histHm(0);
   TH1F *histFi(0);
   TH1F *histFiG(0);
   TH1F *histFiB(0);
   TH1F *histLD(0);
   TH1F *histNn(0);
   TH1F *histNnbfgs(0);
   TH1F *histNnbnn(0);
   TH1F *histNnC(0);
   TH1F *histNnT(0);
   TH1F *histBdt(0);
   TH1F *histBdtG(0);
   TH1F *histBdtB(0);
   TH1F *histBdtD(0);
   TH1F *histBdtF(0);
   TH1F *histRf(0);
   TH1F *histSVMG(0);
   TH1F *histSVMP(0);
   TH1F *histSVML(0);
   TH1F *histFDAMT(0);
   TH1F *histFDAGA(0);
   TH1F *histCat(0);
   TH1F *histPBdt(0);
   TH1F *histDnnGpu(0);
   TH1F *histDnnCpu(0);

   if (Use["Likelihood"])    histLk      = new TH1F( "MVA_Likelihood",    "MVA_Likelihood",    nbin, -1, 1 );
   if (Use["LikelihoodD"])   histLkD     = new TH1F( "MVA_LikelihoodD",   "MVA_LikelihoodD",   nbin, -1, 0.9999 );
   if (Use["LikelihoodPCA"]) histLkPCA   = new TH1F( "MVA_LikelihoodPCA", "MVA_LikelihoodPCA", nbin, -1, 1 );
   if (Use["LikelihoodKDE"]) histLkKDE   = new TH1F( "MVA_LikelihoodKDE", "MVA_LikelihoodKDE", nbin,  -0.00001, 0.99999 );
   if (Use["LikelihoodMIX"]) histLkMIX   = new TH1F( "MVA_LikelihoodMIX", "MVA_LikelihoodMIX", nbin,  0, 1 );
   if (Use["PDERS"])         histPD      = new TH1F( "MVA_PDERS",         "MVA_PDERS",         nbin,  0, 1 );
   if (Use["PDERSD"])        histPDD     = new TH1F( "MVA_PDERSD",        "MVA_PDERSD",        nbin,  0, 1 );
   if (Use["PDERSPCA"])      histPDPCA   = new TH1F( "MVA_PDERSPCA",      "MVA_PDERSPCA",      nbin,  0, 1 );
   if (Use["KNN"])           histKNN     = new TH1F( "MVA_KNN",           "MVA_KNN",           nbin,  0, 1 );
   if (Use["HMatrix"])       histHm      = new TH1F( "MVA_HMatrix",       "MVA_HMatrix",       nbin, -0.95, 1.55 );
   if (Use["Fisher"])        histFi      = new TH1F( "MVA_Fisher",        "MVA_Fisher",        nbin, -4, 4 );
   if (Use["FisherG"])       histFiG     = new TH1F( "MVA_FisherG",       "MVA_FisherG",       nbin, -1, 1 );
   if (Use["BoostedFisher"]) histFiB     = new TH1F( "MVA_BoostedFisher", "MVA_BoostedFisher", nbin, -2, 2 );
   if (Use["LD"])            histLD      = new TH1F( "MVA_LD",            "MVA_LD",            nbin, -2, 2 );
   if (Use["MLP"])           histNn      = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
   if (Use["MLPBFGS"])       histNnbfgs  = new TH1F( "MVA_MLPBFGS",       "MVA_MLPBFGS",       nbin, -1.25, 1.5 );
   if (Use["MLPBNN"])        histNnbnn   = new TH1F( "MVA_MLPBNN",        "MVA_MLPBNN",        nbin, -1.25, 1.5 );
   if (Use["CFMlpANN"])      histNnC     = new TH1F( "MVA_CFMlpANN",      "MVA_CFMlpANN",      nbin,  0, 1 );
   if (Use["TMlpANN"])       histNnT     = new TH1F( "MVA_TMlpANN",       "MVA_TMlpANN",       nbin, -1.3, 1.3 );
   if (Use["DNN_GPU"]) histDnnGpu = new TH1F("MVA_DNN_GPU", "MVA_DNN_GPU", nbin, -0.1, 1.1);
   if (Use["DNN_CPU"]) histDnnCpu = new TH1F("MVA_DNN_CPU", "MVA_DNN_CPU", nbin, -0.1, 1.1);
   if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -0.8, 0.8 );
   if (Use["BDTG"])          histBdtG    = new TH1F( input_name.c_str(),  input_name.c_str(),          nbin, -1.0, 1.0 );
   if (Use["BDTB"])          histBdtB    = new TH1F( "MVA_BDTB",          "MVA_BDTB",          nbin, -1.0, 1.0 );
   if (Use["BDTD"])          histBdtD    = new TH1F( "MVA_BDTD",          "MVA_BDTD",          nbin, -0.8, 0.8 );
   if (Use["BDTF"])          histBdtF    = new TH1F( "MVA_BDTF",          "MVA_BDTF",          nbin, -1.0, 1.0 );
   if (Use["RuleFit"])       histRf      = new TH1F( "MVA_RuleFit",       "MVA_RuleFit",       nbin, -2.0, 2.0 );
   if (Use["SVM_Gauss"])     histSVMG    = new TH1F( "MVA_SVM_Gauss",     "MVA_SVM_Gauss",     nbin,  0.0, 1.0 );
   if (Use["SVM_Poly"])      histSVMP    = new TH1F( "MVA_SVM_Poly",      "MVA_SVM_Poly",      nbin,  0.0, 1.0 );
   if (Use["SVM_Lin"])       histSVML    = new TH1F( "MVA_SVM_Lin",       "MVA_SVM_Lin",       nbin,  0.0, 1.0 );
   if (Use["FDA_MT"])        histFDAMT   = new TH1F( "MVA_FDA_MT",        "MVA_FDA_MT",        nbin, -2.0, 3.0 );
   if (Use["FDA_GA"])        histFDAGA   = new TH1F( "MVA_FDA_GA",        "MVA_FDA_GA",        nbin, -2.0, 3.0 );
   if (Use["Category"])      histCat     = new TH1F( "MVA_Category",      "MVA_Category",      nbin, -2., 2. );
   if (Use["Plugin"])        histPBdt    = new TH1F( "MVA_PBDT",          "MVA_BDT",           nbin, -0.8, 0.8 );

   // PDEFoam also returns per-event error, fill in histogram, and also fill significance
   if (Use["PDEFoam"]) {
      histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
      histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
      histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
   }

   // Book example histogram for probability (the other methods are done similarly)
   TH1F *probHistFi(0), *rarityHistFi(0);
   if (Use["Fisher"]) {
      probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
      rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
   }

   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //
   TFile *input_f(0);
   cout<<"input file:"<<input_name<<endl;
   std::string filename="./"+input_name+".root";
   input_f=TFile::Open(filename.c_str());
   
   // Event loop

   // Prepare the event tree
   // - Here the variable names have to corresponds to your tree
   // - You can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   std::cout << "--- Select signal sample" << std::endl;
   TTree* theTree = (TTree*)input_f->Get("SlimTree");

   theTree->SetBranchAddress( "HT", &HT );
   theTree->SetBranchAddress( "ttc_region", &ttc_region );
   theTree->SetBranchAddress( "j1_FlavCvB", &j1_FlavCvB );
   theTree->SetBranchAddress( "j1_FlavCvL", &j1_FlavCvL );
   theTree->SetBranchAddress( "j2_FlavCvB", &j2_FlavCvB );
   theTree->SetBranchAddress( "j2_FlavCvL", &j2_FlavCvL );
   theTree->SetBranchAddress( "j3_FlavCvB", &j3_FlavCvB );
   theTree->SetBranchAddress( "j3_FlavCvL", &j3_FlavCvL );
   theTree->SetBranchAddress( "dr_j1j2", &dr_j1j2);
   theTree->SetBranchAddress( "dr_j1j3", &dr_j1j3);
   theTree->SetBranchAddress( "dr_j2j3", &dr_j2j3);
   theTree->SetBranchAddress( "ttc_l1_pt", &ttc_l1_pt );
   theTree->SetBranchAddress( "ttc_l2_pt", &ttc_l2_pt);
   theTree->SetBranchAddress( "ttc_l1_eta", &ttc_l1_eta );
   theTree->SetBranchAddress( "ttc_l2_eta", &ttc_l2_eta);
   theTree->SetBranchAddress( "ttc_met", &ttc_met);
   theTree->SetBranchAddress( "ttc_met_phi", &ttc_met_phi);
   theTree->SetBranchAddress( "ttc_mll", &ttc_mll);
   theTree->SetBranchAddress( "ttc_mllj1", &ttc_mllj1);
   theTree->SetBranchAddress( "ttc_mllj2", &ttc_mllj2);
   theTree->SetBranchAddress( "ttc_mllj3", &ttc_mllj3);

   // Efficiency calculator for cut method
   Int_t    nSelCutsGA = 0;
   Double_t effS       = 0.7;

   std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();

   // store normalization for no ctag
   float ctag_norm=0.;

   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

      if (ievt%40000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      theTree->GetEntry(ievt);
      if(channel=="mm" && ttc_region!=1) continue;
      if(channel=="em" && ttc_region!=2) continue;
      if(channel=="ee" && ttc_region!=3) continue;
      if(channel=="ee" && ttc_mll>60 && ttc_mll<120) continue;
      // Return the MVA outputs and fill into histograms

      if (Use["CutsGA"]) {
         // Cuts is a special case: give the desired signal efficiency
         Bool_t passed = reader->EvaluateMVA( "CutsGA method", effS );
         if (passed) nSelCutsGA++;
      }

      if (Use["Likelihood"   ])   histLk     ->Fill( reader->EvaluateMVA( "Likelihood method"    ) );
      if (Use["LikelihoodD"  ])   histLkD    ->Fill( reader->EvaluateMVA( "LikelihoodD method"   ) );
      if (Use["LikelihoodPCA"])   histLkPCA  ->Fill( reader->EvaluateMVA( "LikelihoodPCA method" ) );
      if (Use["LikelihoodKDE"])   histLkKDE  ->Fill( reader->EvaluateMVA( "LikelihoodKDE method" ) );
      if (Use["LikelihoodMIX"])   histLkMIX  ->Fill( reader->EvaluateMVA( "LikelihoodMIX method" ) );
      if (Use["PDERS"        ])   histPD     ->Fill( reader->EvaluateMVA( "PDERS method"         ) );
      if (Use["PDERSD"       ])   histPDD    ->Fill( reader->EvaluateMVA( "PDERSD method"        ) );
      if (Use["PDERSPCA"     ])   histPDPCA  ->Fill( reader->EvaluateMVA( "PDERSPCA method"      ) );
      if (Use["KNN"          ])   histKNN    ->Fill( reader->EvaluateMVA( "KNN method"           ) );
      if (Use["HMatrix"      ])   histHm     ->Fill( reader->EvaluateMVA( "HMatrix method"       ) );
      if (Use["Fisher"       ])   histFi     ->Fill( reader->EvaluateMVA( "Fisher method"        ) );
      if (Use["FisherG"      ])   histFiG    ->Fill( reader->EvaluateMVA( "FisherG method"       ) );
      if (Use["BoostedFisher"])   histFiB    ->Fill( reader->EvaluateMVA( "BoostedFisher method" ) );
      if (Use["LD"           ])   histLD     ->Fill( reader->EvaluateMVA( "LD method"            ) );
      if (Use["MLP"          ])   histNn     ->Fill( reader->EvaluateMVA( "MLP method"           ) );
      if (Use["MLPBFGS"      ])   histNnbfgs ->Fill( reader->EvaluateMVA( "MLPBFGS method"       ) );
      if (Use["MLPBNN"       ])   histNnbnn  ->Fill( reader->EvaluateMVA( "MLPBNN method"        ) );
      if (Use["CFMlpANN"     ])   histNnC    ->Fill( reader->EvaluateMVA( "CFMlpANN method"      ) );
      if (Use["TMlpANN"      ])   histNnT    ->Fill( reader->EvaluateMVA( "TMlpANN method"       ) );
      if (Use["DNN_GPU"]) histDnnGpu->Fill(reader->EvaluateMVA("DNN_GPU method"));
      if (Use["DNN_CPU"]) histDnnCpu->Fill(reader->EvaluateMVA("DNN_CPU method"));
      if (Use["BDT"          ])   histBdt    ->Fill( reader->EvaluateMVA( "BDT method"           ) );
      if (Use["BDTG"         ])
       { 
	histBdtG->Fill( reader->EvaluateMVA( "BDTG method"), 1.0);
	}
      if (Use["BDTB"         ])   histBdtB   ->Fill( reader->EvaluateMVA( "BDTB method"          ) );
      if (Use["BDTD"         ])   histBdtD   ->Fill( reader->EvaluateMVA( "BDTD method"          ) );
      if (Use["BDTF"         ])   histBdtF   ->Fill( reader->EvaluateMVA( "BDTF method"          ) );
      if (Use["RuleFit"      ])   histRf     ->Fill( reader->EvaluateMVA( "RuleFit method"       ) );
      if (Use["SVM_Gauss"    ])   histSVMG   ->Fill( reader->EvaluateMVA( "SVM_Gauss method"     ) );
      if (Use["SVM_Poly"     ])   histSVMP   ->Fill( reader->EvaluateMVA( "SVM_Poly method"      ) );
      if (Use["SVM_Lin"      ])   histSVML   ->Fill( reader->EvaluateMVA( "SVM_Lin method"       ) );
      if (Use["FDA_MT"       ])   histFDAMT  ->Fill( reader->EvaluateMVA( "FDA_MT method"        ) );
      if (Use["FDA_GA"       ])   histFDAGA  ->Fill( reader->EvaluateMVA( "FDA_GA method"        ) );
      if (Use["Category"     ])   histCat    ->Fill( reader->EvaluateMVA( "Category method"      ) );
      if (Use["Plugin"       ])   histPBdt   ->Fill( reader->EvaluateMVA( "P_BDT method"         ) );

      // Retrieve also per-event error
      if (Use["PDEFoam"]) {
         Double_t val = reader->EvaluateMVA( "PDEFoam method" );
         Double_t err = reader->GetMVAError();
         histPDEFoam   ->Fill( val );
         histPDEFoamErr->Fill( err );
         if (err>1.e-50) histPDEFoamSig->Fill( val/err );
      }

      // Retrieve probability instead of MVA output
      if (Use["Fisher"])   {
         probHistFi  ->Fill( reader->GetProba ( "Fisher method" ) );
         rarityHistFi->Fill( reader->GetRarity( "Fisher method" ) );
      }
   }


   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   // Get efficiency for cuts classifier
   if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                                << " (for a required signal efficiency of " << effS << ")" << std::endl;

   if (Use["CutsGA"]) {

      // test: retrieve cuts for particular signal efficiency
      // CINT ignores dynamic_casts so we have to use a cuts-specific Reader function to acces the pointer
      TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;

      if (mcuts) {
         std::vector<Double_t> cutsMin;
         std::vector<Double_t> cutsMax;
         mcuts->GetCuts( 0.7, cutsMin, cutsMax );
         std::cout << "--- -------------------------------------------------------------" << std::endl;
         std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
         for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
            std::cout << "... Cut: "
                      << cutsMin[ivar]
                      << " < \""
                      << mcuts->GetInputVar(ivar)
                      << "\" <= "
                      << cutsMax[ivar] << std::endl;
         }
         std::cout << "--- -------------------------------------------------------------" << std::endl;
      }
   }

   delete reader;

   std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
   string name_temp="ttc2018_"+input_name;
   histBdtG->SetNameTitle(name_temp.c_str(),name_temp.c_str());

   input_f->Close();
   return histBdtG;

}

int TMVAClassificationApplication()
{
   std::cout<<"START!"<<std::endl;
   string sample_path="/eos/cms/store/group/phys_top/ExtraYukawa/2018/";

   string type_="S0"; 
   string cp="rtc10";
   string mass="400";

   float eff_N_signal;
   string flags="GenModel_T"+type_+"ToTTQ_M"+type_+"_400_TuneCP5_13TeV_G2HDM_"+cp+"_madgraphMLM_pythia8";

   string ntemp_signal=sample_path;
   if(type_=="A") ntemp_signal=ntemp_signal+"ttc_a_"+cp+".root";
   if(type_=="S0") ntemp_signal=ntemp_signal+"ttc_s0_"+cp+".root";
   TFile*ftemp_signal=TFile::Open(ntemp_signal.c_str());
   TH1D*htemp_signal=(TH1D*)ftemp_signal->Get("nEventsGenWeighted");
   TTree*ttemp_signal=(TTree*)ftemp_signal->Get("Events");
   int nsignal_total=ttemp_signal->GetEntriesFast();
   eff_N_signal=(0.5*ttemp_signal->GetEntries(flags.c_str())*(htemp_signal->GetBinContent(1))/nsignal_total);
   ftemp_signal->Close();

   vector<std::string> samples;
   samples.push_back("TTTo2L");
   samples.push_back("ttWW");
   samples.push_back("ttWZ");
   samples.push_back("ttZZ");
   samples.push_back("ttWtoLNu");
   samples.push_back("ttZ");
   samples.push_back("ttZtoQQ");
   samples.push_back("tttW");
   samples.push_back("tttt");
   samples.push_back("tzq");
   samples.push_back("WWW");
   samples.push_back("DY");
   samples.push_back("WWZ");
   samples.push_back("WWdps");
   samples.push_back("WZ");
   samples.push_back("WZZ");
   samples.push_back("ZZZ");
   samples.push_back("osWW");
   samples.push_back("tW");
   samples.push_back("tbarW");
   samples.push_back("ttH");
   samples.push_back("ttWH");
   samples.push_back("ttWtoQQ");
   samples.push_back("ttZH");
   samples.push_back("tttJ");
   samples.push_back("zz2l");

   vector<float> xss;
   xss.push_back(88.3419);//TTTo2L
   xss.push_back(0.007003);//ttWW
   xss.push_back(0.002453);//ttWZ
   xss.push_back(0.001386);//ttZZ
   xss.push_back(0.1792);//ttWtoLNu
   xss.push_back(0.2589);//ttZ
   xss.push_back(0.6012);//ttZtoQQ
   xss.push_back(0.00073);//tttW
   xss.push_back(0.0082);//tttt
   xss.push_back(0.07561);//tzq
   xss.push_back(0.2086);//WWW
   xss.push_back(6077.22);//DY
   xss.push_back(0.1707);//WWZ
   xss.push_back(1.62);//WWdps
   xss.push_back(5.213);//WZ
   xss.push_back(0.05709);//WZZ
   xss.push_back(0.01476);//ZZZ
   xss.push_back(11.09);//osWW
   xss.push_back(35.85);//tW
   xss.push_back(35.85);//tbarW
   xss.push_back(0.5269);//ttH
   xss.push_back(0.00114);//ttWH
   xss.push_back(0.3708);//ttWtoQQ
   xss.push_back(0.00113);//ttZH
   xss.push_back(0.0004);//tttJ
   xss.push_back(0.9738);//zz2l
   
   vector<float> eff_N;
   for(int i=0;i<26;i++){
     string ntemp=sample_path+samples[i]+".root";
     if(i==11)ntemp=sample_path+"DYnlo.root";
     TFile*ftemp=TFile::Open(ntemp.c_str());
     TH1D*ttemp=(TH1D*)ftemp->Get("nEventsGenWeighted");
     if (i<11)eff_N.push_back(0.5*ttemp->GetBinContent(1));//half of the events are used for BDT training
     else eff_N.push_back(ttemp->GetBinContent(1));
     ftemp->Close();
   }
//   eff_N.push_back(3.4087701e+08);
//   eff_N.push_back(349000.);
//   eff_N.push_back(175000.);
//   eff_N.push_back(1935527);
//   eff_N.push_back(3455733.0);
   
  string weights[52]={"nominal_noctag","central","pileup_up","pileup_down","muID_sysup","muID_sysdown","muID_statup","muID_statdown","eleID_sysup","eleID_sysdown","eleID_statup","eleID_statdown","trigger_up","trigger_down","lumi_up","lumi_down","ctag_statup","ctag_statdo","ctag_EleIDup","ctag_EleIDdo","ctag_LHEScaleWeightmuFup","ctag_LHEScaleWeightmuFdo","ctag_LHEScaleWeightmuRup","ctag_LHEScaleWeightmuRdo","ctag_MuIDup","ctag_MuIDdo","ctag_PSWeightFSRup","ctag_PSWeightFSRdo","ctag_PUWeightup","ctag_PUWeightdo","ctag_XSec_DYJetsup","ctag_XSec_DYJetsdo","ctag_XSec_STup","ctag_XSec_STdo","ctag_XSec_VVup","ctag_XSec_VVdo","ctag_XSec_WJetsup","ctag_XSec_WJetsdo","ctag_XSec_ttbarup","ctag_XSec_ttbardo","ctag_jerup","ctag_jerdo","ctag_jesTotalup","ctag_jesTotaldo","charFlip_SFup","charFlip_SFdo","sig_pdfup","sig_pdfdo","sig_scaleup","sig_scaledo","sig_psup","sig_psdo"};

  string system_unc[6]={"jesup","jesdo","jerup","jerdo","unclusterEup","unclusterEdo"};
  string channels[3]={"ee","em","mm"};

  for (int ic=0;ic<3;ic++){
    string output_name="TMVApp_"+mass+"_"+channels[ic]+".root";
    TFile *target  = new TFile( output_name.c_str(),"RECREATE" );
    string signal_input="";
    if(type_=="A")signal_input=signal_input+"ttc_a_"+cp+"_M"+type_+mass;
    if(type_=="S0")signal_input=signal_input+"ttc_s0_"+cp+"_M"+type_+mass;
    TH1F*htemp;
    TH1F*hfake_no_mcsubtraction;
    TH1F*hfake;
    TH1F*hfake_up;
    TH1F*hfake_down;

    std::vector<float> ctagnorms;
    std::cout<<"start looping weights"<<std::endl;
    for(int iw=0;iw<52;iw++){
      // signal don't need charge flip SF
      if(!(weights[iw].find("charFlip")!= string::npos))
      {
        htemp=Getoutput("",signal_input,1.0,eff_N_signal,weights[iw],mass,channels[ic],type_,cp);
        if(iw==0)ctagnorms.push_back(htemp->Integral());
        if(iw>0 && weights[iw].find("ctag")!= string::npos)htemp->Scale(ctagnorms[0]/htemp->Integral());
        target->cd();
        htemp->Write();
      }
      //for bkgs, no need of signal theoretic uncertainty
      if(weights[iw].find("sig")!= string::npos)continue;
      for(int is=0;is<26;is++){
        cout<<"start loop process:"<<samples[is]<<endl;
        htemp=Getoutput("",samples[is],xss[is],eff_N[is],weights[iw],mass,channels[ic],type_,cp);
        if(iw==0)ctagnorms.push_back(htemp->Integral());
        if(iw>0 && weights[iw].find("ctag")!= string::npos)htemp->Scale(ctagnorms[is+1]/htemp->Integral());
        target->cd();
        htemp->Write();
     }
   }

   std::cout<<"start looping systematics"<<std::endl;
   for(int isys=0;isys<6;isys++){
      htemp=Getoutput_sys("",signal_input,1.0,eff_N_signal,system_unc[isys],mass,channels[ic],type_,cp);
      target->cd();
      htemp->Write();
      for(int is=0;is<26;is++){
        htemp=Getoutput_sys("",samples[is],xss[is],eff_N[is],system_unc[isys],mass,channels[ic],type_,cp);
        target->cd();
        htemp->Write();
     }
   }
  
   std::cout<<"start looping fake"<<std::endl;
   if(ic==0) {
     hfake_no_mcsubtraction=Getoutput_fakelep("","fakelep_ee",1.,1,mass,channels[ic],type_,cp);
   }
   if(ic==1) {
     hfake_no_mcsubtraction=Getoutput_fakelep("","fakelep_em",1.,1,mass,channels[ic],type_,cp);
   }
   if(ic==2) {
     hfake_no_mcsubtraction=Getoutput_fakelep("","fakelep_mm",1.,1,mass,channels[ic],type_,cp);
   }

   hfake=(TH1F*)hfake_no_mcsubtraction->Clone();

   std::cout<<"start looping fake mc"<<std::endl;
   for(int is=0;is<26;is++){
     htemp=Getoutput_fake_mc("",samples[is]+"_fake_"+channels[ic],xss[is],eff_N[is],"central",mass,channels[ic],type_,cp);
     hfake->Add(htemp);
   }

   hfake_up=(TH1F*)hfake->Clone();
   hfake_up->Scale(1.3);
   hfake_down=(TH1F*)hfake->Clone();
   hfake_down->Scale(0.7);

   hfake_no_mcsubtraction->SetNameTitle("ttc2018_TTTo1L_noMCsub","ttc2018_TTTo1L_noMCsub");
   hfake->SetNameTitle("ttc2018_TTTo1L","ttc2018_TTTo1L");
   hfake_up->SetNameTitle("ttc2018_TTTo1L_fakeUp","ttc2018_TTTo1L_fakeUp");
   hfake_down->SetNameTitle("ttc2018_TTTo1L_fakeDown","ttc2018_TTTo1L_fakeDown");
   target->cd();
   hfake->Write();
   hfake_up->Write();
   hfake_down->Write();

   // check fake subtraction of MC
   std::cout<<"WITHOUT subtraction"<<hfake_no_mcsubtraction->Integral()<<std::endl;
   std::cout<<"WITH subtractiob"<<hfake->Integral()<<std::endl;

   TH1F*hdata;
   std::cout<<"start looping data"<<std::endl;
   if(ic==0) {
     hdata=Getoutput_data("","data_ee",1.,1,mass,channels[ic],type_,cp);
   }
   if(ic==1) {
     hdata=Getoutput_data("","data_em",1.,1,mass,channels[ic],type_,cp);
   }
   if(ic==2) {
     hdata=Getoutput_data("","data_mm",1.,1,mass,channels[ic],type_,cp);
   }
   hdata->SetNameTitle("ttc2018_data_obs","ttc2018_data_obs");
   target->cd();
   hdata->Write();
   target->Close();
  }
  
  return 0;
}
