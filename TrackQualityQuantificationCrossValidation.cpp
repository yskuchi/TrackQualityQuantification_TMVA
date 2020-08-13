#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/CrossValidation.h"
#include "TMVA/HyperParameterOptimisation.h"

using namespace TMVA;
int TrackQualityQuantificationCrossValidation(bool useRandomSplitting = true)
{
   // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
   // if you use your private .rootrc, or run from a different directory, please copy the
   // corresponding lines from .rootrc
   // methods to be processed can be given as an argument; use format:
   //
   //     mylinux~> root  ./TrackQualityQuantificationCrossValidation.cpp
   //
   //---------------------------------------------------------------

   // This loads the library
   TMVA::Tools::Instance();


   std::cout << std::endl;
   std::cout << "==> Start TMVACrossValidation" << std::endl;

   // Here the preparation phase begins

   // Create a new root output file
   TString outfileName( "TrackQualityQuantificationCrossValidation.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //
   //     (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //     (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";
   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
   dataloader->AddVariable( "ngoodhits", "Number of good DCHHits", "", 'I' );
   dataloader->AddVariable( "redChi2", "Chi2/NDF of DCH tracking", "", 'F' );
   dataloader->AddVariable( "EPositronUncert", "Uncertainty E", "GeV", 'F' );
   dataloader->AddVariable( "ThetaPositronUncert", "Uncertainty theta", "degree", 'F' );
   dataloader->AddVariable( "PhiPositronUncert", "Uncertainty phi", "degree", 'F' );
   dataloader->AddVariable( "extrapolationLengthTarget", "Length between target to 1st DCHHit", "cm", 'F' );
   dataloader->AddVariable( "extrapolationLengthSPX", "Length between last DCHHit to SPX", "cm", 'F' );
   dataloader->AddVariable( "nSPXHits", "Number of good SPXHits", "", 'I' );
   dataloader->AddVariable( "matchingChi2", "Chi2 of matching", "", 'F' );
   dataloader->AddVariable( "matchingDT", "Time difference b/w DCHTrack and SPXTrack", "s", 'F' );

   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables
   // dataloader->AddSpectator( "spec1:=var1*2",  "Spectator 1", "units", 'F' );
   // dataloader->AddSpectator( "spec2:=var1*3",  "Spectator 2", "units", 'F' );


   // Read training and test data (see TMVACrossValidation for reading ASCII files)
   // load the signal and background event samples from ROOT trees
   TFile *input(0);
   TString fname = "./TrackQualityQuantificationInput.root";
   if (!gSystem->AccessPathName( fname )) {
      input = TFile::Open( fname ); // check if file in local directory exists
   }
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVACrossValidation           : Using input file: " << input->GetName() << std::endl;

   // --- Register the training and test trees
   TTree* t_sig = (TTree*)input->Get("treeGood");
   TTree* t_bkg = (TTree*)input->Get("treeBad");

   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 6.0;
    
   // You can add an arbitrary number of signal or background trees
   dataloader->AddSignalTree(t_sig, signalWeight);
   dataloader->AddBackgroundTree(t_bkg, backgroundWeight);

   // Set individual event weights (the variables must exist in the original TTree)
   dataloader->SetBackgroundWeightExpression("weight");

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";
 
   // The CV mechanism of TMVA splits up the training set into several folds.
   // The test set is currently left unused. The `nTest_ClassName=1` assigns
   // one event to the the test set for each class and puts the rest in the
   // training set. A value of 0 is a special value and would split the
   // datasets 50 / 50.
   dataloader->PrepareTrainingAndTestTree("", "",
                                          "nTest_Signal=1"
                                          ":nTest_Background=1"
                                          ":SplitMode=Random"
                                          ":NormMode=NumEvents"
                                          ":!V");
   //
   //     dataloader->PrepareTrainingAndTestTree( mycut,
   //            "nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V" );
   // If no numbers of events are given, half of the events in the tree are used
   // for training, and the other half for testing:
   //
   //     dataloader->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );

   // --------------------------------------------------------------------------
 
   //
   // This sets up a CrossValidation class (which wraps a TMVA::Factory
   // internally) for k-fold cross validation.
   //
   // The split type can be "Random", "RandomStratified" or "Deterministic".
   // For the last option, check the comment below. Random splitting randomises
   // the order of events and distributes events as evenly as possible.
   // RandomStratified applies the same logic but distributes events within a
   // class as evenly as possible over the folds.
   //
   UInt_t numFolds = 4;
   TString analysisType = "Classification";
   
   TString splitType = (useRandomSplitting) ? "Random" : "Deterministic";
 
   //
   // One can also use a custom splitting function for producing the folds.
   // The example uses a dataset spectator `eventID`.
   //
   // The idea here is that eventID should be an event number that is integral,
   // random and independent of the data, generated only once. This last
   // property ensures that if a calibration is changed the same event will
   // still be assigned the same fold.
   // 
   // This can be used to use the cross validated classifiers in application,
   // a technique that can simplify statistical analysis.
   // 
   // If you want to run TMVACrossValidationApplication, make sure you have 
   // run this tutorial with Deterministic splitting type, i.e.
   // with the option useRandomSPlitting = false
   // 
 
   TString splitExpr = (!useRandomSplitting) ? "int(fabs([eventID]))%int([NumFolds])" : "";
 
   TString cvOptions = Form("!V"
                            ":!Silent"
                            ":ModelPersistence"
                            ":AnalysisType=%s"
                            ":SplitType=%s"
                            ":NumFolds=%i"
                            ":SplitExpr=%s",
                            analysisType.Data(), splitType.Data(), numFolds,
                            splitExpr.Data());
 
   TMVA::CrossValidation cv{"TrackQualityQuantificationCrossValidation", dataloader, outputFile, cvOptions};
 

   // --------------------------------------------------------------------------
   // Hyper-parameter tuning
   // Currently only available for (adaBoost) BDT and SVM


   TMVA::HyperParameterOptimisation * HPO = new
         TMVA::HyperParameterOptimisation(dataloader);
   // HPO->BookMethod(TMVA::Types::kBDT, "BDT",
   //                 "");
   HPO->BookMethod(TMVA::Types::kSVM, "SVM",
                   "");
   // HPO->SetNumFolds(numFolds);
   // HPO->SetFitter("Minuit");
   // Available FOM types:
   // "Separation"
   // "ROCIntegral"
   // "SigEffAtBkgEff01"
   // "SigEffAtBkgEff001"
   // "SigEffAtBkgEff002"
   // "BkgRejAtSigEff05"
   // "BkgEffAtSigEff05"
   // HPO->SetFOMType("Separation");
   cout<<"Evaluate!!"<<endl;
   HPO->Evaluate();
   cout<<"Print!!"<<endl;
   const TMVA::HyperParameterOptimisationResult & HPOResult = HPO->GetResults();
   HPOResult.Print();


   // --------------------------------------------------------------------------
 
   //
   // Books a method to use for evaluation
   //
   // cv.BookMethod(TMVA::Types::kBDT, "BDTG",
   //               "!H:!V:NTrees=100:MinNodeSize=2.5%:BoostType=Grad"
   //               ":NegWeightTreatment=Pray:Shrinkage=0.10:nCuts=20"
   //               ":MaxDepth=2");
   
   cv.BookMethod(TMVA::Types::kBDT, "BDT",
                 "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:"
                 "BoostType=AdaBoost:AdaBoostBeta=0.5:"
                 "UseBaggedBoost:BaggedSampleFraction=0.5:"
                 "SeparationType=GiniIndex:nCuts=20" );
  
   // cv.BookMethod(TMVA::Types::kFisher, "Fisher",
   //               "!H:!V:Fisher:VarTransform=None");
 
   // --------------------------------------------------------------------------
 
   //
   // Train, test and evaluate the booked methods.
   // Evaluates the booked methods once for each fold and aggregates the result
   // in the specified output file.
   //
   cv.Evaluate();
 
   // --------------------------------------------------------------------------
 
   //
   // Process some output programatically, printing the ROC score for each
   // booked method.
   //
   size_t iMethod = 0;
   for (auto && result : cv.GetResults()) {
      std::cout << "Summary for method " << cv.GetMethods()[iMethod++].GetValue<TString>("MethodName")
                << std::endl;
      for (UInt_t iFold = 0; iFold<cv.GetNumFolds(); ++iFold) {
         std::cout << "\tFold " << iFold << ": "
                   << "ROC int: " << result.GetROCValues()[iFold]
                   << ", "
                   << "BkgEff@SigEff=0.3: " << result.GetEff30Values()[iFold]
                   << std::endl;
      }
   }
 
   // --------------------------------------------------------------------------
 
   //
   // Save the output
   //
   outputFile->Close();
 
   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVACrossValidation is done!" << std::endl;
 
   // --------------------------------------------------------------------------
 
   //
   // Launch the GUI for the root macros
   //
   if (!gROOT->IsBatch()) {
      // Draw cv-specific graphs
      // cv.GetResults()[0].DrawAvgROCCurve(kTRUE, "Avg ROC for BDTG");
      // cv.GetResults()[0].DrawAvgROCCurve(kTRUE, "Avg ROC for Fisher");
      cv.GetResults()[0].DrawAvgROCCurve(kTRUE, "Avg ROC for BDT");
 
      // You can also use the classical gui
      TMVA::TMVAGui(outfileName);
   }
 
   return 0;
}

//
// This is used if the macro is compiled. If run through ROOT with
// `root -l -b -q MACRO.C` or similar it is unused.
//
int main( int argc, char** argv )
{

   TrackQualityQuantificationCrossValidation();
   return 0;
}
