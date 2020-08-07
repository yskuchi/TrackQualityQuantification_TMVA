#include <macros/positron/setting.h>
#include <TF1.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TMath.h>
#include <TChain.h>
#include <TStyle.h>
#include <TFile.h>
#include <TVector3.h>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include <set>
#include "../common/include/units/MEGSystemOfUnits.h"
#include "../common/include/units/MEGPhysicalConstants.h"
#if !defined(__CLING__) || defined(__ROOTCLING__)
#   include <ROMEString.h>
#   include <ROMETreeInfo.h>
#   include "include/generated/MEGAllFolders.h"
#   include "include/glb/MEGPhysicsSelection.h"
#else
class ROMETreeInfo;
class MEGMCKineEvent;
class MEGMCTrackEvent;
class MEGMCSPXEvent;
#endif
using namespace MEG;

void TruthCheck(vector<Bool_t> &signalDCHTrackFlag, Short_t signalSevID, TClonesArray* pDCHTrackArray, 
                TClonesArray* pDCHHitArray, MEGMCCYLDCHEvent* pMCDCH, TClonesArray* pSPXClusterArray, 
                TClonesArray* pSPXHitArray, MEGMCSPXEvent* pMCSPX);
Double_t FitHist(TH1* hist0, TH1* hist1);
vector<Bool_t> AND(vector<Bool_t> &a, vector<Bool_t> &b)
{
   vector<Bool_t> c(a.size());
   for (size_t i = 0; i < a.size(); i++) {
      c[i] = a[i] && b[i];
   }
   return c;
}

//______________________________________________________________________________
void TrackQualificationDataset() 
{
/*
  This macro make a tree file to be used to define the track-qualification variable

  The output ROOT file TrackQualificationInput.root contains a tree 'reg', 'treeGood', and 
  'treeBad' that contain input variables and the target variables for the regression
  and classification.

  Usage:
  - your setup should be written in macros/positron/setting.h
  - cd $MEG2SYS/analyzer
  - ./meganalyzer -I path/to/track_qualification_TMVA/TrackQualificationDataset.cpp
  -  or
  - ./meganalyzer -I '../../track_qualification_TMVA/TrackQualificationDataset.cpp+'
  - rec/sim files should include signal events.

*/

   //-----------------------
   //Variable, branch adrress and histgram setting
   //-----------------------

   // number of events
   Int_t nGammaAcceptance = 0;
   Int_t nEffTCHits = 0;
   Int_t nEffEvents=0;
   Int_t nTrueSignal=0;
   
   // Set Style
   gStyle->SetOptStat(1110);
   // gStyle->SetOptStat(0);
   gStyle->SetStatW(0.3);
   gStyle->SetOptFit(11);
   gStyle->SetOptTitle(1);
   gStyle->SetPadGridX(0);
   gStyle->SetPadGridY(0);
   gStyle->SetPadTickX(0);
   gStyle->SetPadTickY(0);
   gStyle->SetPadTopMargin(0.1);
   gStyle->SetPadBottomMargin(0.15);
   gStyle->SetPadRightMargin(0.15);
   gStyle->SetNdivisions(404, "X");
   gStyle->SetNdivisions(406, "Y");
   gStyle->SetTitleSize(0.05, "X");
   gStyle->SetTitleSize(0.05, "Y");
   gStyle->SetLabelSize(0.06, "X");
   gStyle->SetLabelSize(0.06, "Y");


   // Output tree file
   TFile output("TrackQualificationInput.root", "RECREATE", "TrackQualification outputs");

   // Output tree
   TTree *outputTree = new TTree("reg", "Tree for regression");
   TTree *outputGoodTree = new TTree("treeGood", "Tree for classification good sample");
   TTree *outputBadTree = new TTree("treeBad", "Tree for classification bad sample");
   vector<TTree*> trees {outputTree, outputGoodTree, outputBadTree};

   // Input variables
   Float_t weight_(1.);
   Bool_t  good_;
   Short_t ngoodhits_;
   Float_t redChi2_;
   Float_t EPositronUncert_;
   Float_t PhiPositronUncert_;
   Float_t ThetaPositronUncert_;
   Float_t ZPositronUncert_;
   Float_t YPositronUncert_;
   Float_t extrapolationLengthTarget_;
   Float_t extrapolationLengthSPX_;
   Short_t nSPXHits_;
   Float_t matchingDT_;
   Float_t matchingDW_;
   Float_t matchingDV_;
   Float_t matchingChi2_;
   // Target variable
   Float_t target_;
   for (auto &&tree : trees) {
      tree->Branch("ngoodhits", &ngoodhits_, "ngoodhits/s");
      tree->Branch("redChi2", &redChi2_, "redChi2/F");
      tree->Branch("EPositronUncert", &EPositronUncert_, "EPositronUncert/F");
      tree->Branch("PhiPositronUncert", &PhiPositronUncert_, "PhiPositronUncert/F");
      tree->Branch("ThetaPositronUncert", &ThetaPositronUncert_, "ThetaPositronUncert/F");
      tree->Branch("ZPositronUncert", &ZPositronUncert_, "ZPositronUncert/F");
      tree->Branch("YPositronUncert", &YPositronUncert_, "YPositronUncert/F");
      tree->Branch("extrapolationLengthTarget", &extrapolationLengthTarget_, "extrapolationLengthTarget/F");
      tree->Branch("extrapolationLengthSPX", &extrapolationLengthSPX_, "extrapolationLengthSPX/F");
      tree->Branch("nSPXHits", &nSPXHits_, "nSPXHits/s");
      tree->Branch("matchingDT", &matchingDT_, "matchingDT/F");
      tree->Branch("matchingDW", &matchingDW_, "matchingDW/F");
      tree->Branch("matchingDV", &matchingDV_, "matchingDV/F");
      tree->Branch("matchingChi2", &matchingChi2_, "matchingChi2/F");
      tree->Branch("good", &good_, "good/O");
      tree->Branch("weight", &weight_, "weight/F");
      tree->Branch("target", &target_, "target/F");   
   }

   // Open the files
   TChain *sim = new TChain("sim");
   TFile *file0 = 0;
   TFile *file0_sim = 0;
   TChain *rec = new TChain("rec");
   
   for (Int_t iRun = 0; iRun < nfile; iRun++) {
      Int_t iRunNumber = iRun + sRun;
      // path to your rec/sim file
      TString recfile = Form(inputrecdir + "rec%05d.root", iRunNumber);
      TString simfile = Form(inputsimdir + "sim%05d.root", iRunNumber);

      sim->Add(simfile.Data());
      rec->Add(recfile.Data());
      
      if (!file0_sim) {
         file0_sim = TFile::Open(simfile.Data());
      }
      if (!file0) {
         file0 = TFile::Open(recfile.Data());
      }
      
   }

   TBranch *bMCTrack;
   TClonesArray *pMCTrackArray = new TClonesArray("MEGMCTrackEvent");
   MEGMCTrackEvent *pMCTrack;
   TBranch *bMCKine;
   MEGMCKineEvent *pMCKine = new MEGMCKineEvent();
   TBranch *bMCSPX;
   MEGMCSPXEvent *pMCSPX = new MEGMCSPXEvent();
   MEGSPXMCHit *pSPXMCHit;
   TBranch *bMCMix;
   MEGMCMixtureInfoEvent *pMCMix = new MEGMCMixtureInfoEvent();
   MEGMCCYLDCHEvent *pMCDCH = new MEGMCCYLDCHEvent();
   TBranch *bMCDCH;
   MEGCYLDCHMCHit *pDCHMCHit;
   sim->SetBranchStatus("*", 0); 
   sim->SetBranchAddress("mckine.", &pMCKine, &bMCKine);
   sim->SetBranchStatus("mckine.*", 1);
   sim->SetBranchAddress("mcspx.", &pMCSPX, &bMCSPX);
   sim->SetBranchStatus("mcspx.*", 1);
   sim->SetBranchAddress("mccyldch.", &pMCDCH, &bMCDCH);
   sim->SetBranchStatus("mccyldch.*", 1);
   sim->SetBranchAddress("mcmixevent.", &pMCMix, &bMCMix);
   sim->SetBranchStatus("mcmixevent.*", 1);
   
   TBranch *bInfoS_sim;
   ROMETreeInfo *pInfoS_sim;
   pInfoS_sim = new ROMETreeInfo;
   sim->SetBranchAddress("Info.", &pInfoS_sim, &bInfoS_sim);
   sim->SetBranchStatus("Info.", 1);
   sim->SetBranchStatus("Info.*", 1);

   TBranch *bEventHeader;
   MEGEventHeader *pEventHeader;
   TBranch *bSPXCluster;
   TClonesArray *pSPXClusterArray = new TClonesArray("MEGSPXCluster");
   MEGSPXCluster *pSPXCluster;
   TBranch *bSPXTrack;
   TClonesArray *pSPXTrackArray = new TClonesArray("MEGSPXTrack");
   MEGSPXTrack *pSPXTrack;
   TBranch *bMatchedTrack;
   TClonesArray *pMatchedTrackArray = new TClonesArray("MEGDCHSPXMatchedTrack");
   MEGDCHSPXMatchedTrack *pMatchedTrack;
   TBranch *bSPXHit;
   TClonesArray *pSPXHitArray = new TClonesArray("MEGSPXHit");
   MEGSPXHit *pSPXHit;
   TBranch *bDCHTrackCand;
   TClonesArray *pDCHTrackCandArray = new TClonesArray("MEGDCHTrackCandidate");
   MEGDCHTrackCandidate *pDCHTrackCand;
   TBranch *bDCHTrack;
   TClonesArray *pDCHTrackArray = new TClonesArray("MEGDCHTrack");
   MEGDCHTrack *pDCHTrack;
   TBranch *bDCHHit;
   TClonesArray *pDCHHitArray = new TClonesArray("MEGDCHHit");
   MEGDCHHit *pDCHHit;
   TBranch *bInfoS;
   ROMETreeInfo *pInfoS;
   TClonesArray *pPixelRunHeaders = (TClonesArray*)gROOT->FindObject("SPXPixelRunHeader");
   MEGSPXPixelRunHeader *pPixelRunHeader;
   MEGTargetRunHeader *pTargetRunHeader = (MEGTargetRunHeader*)gROOT->FindObject("TargetRunHeader");
   pInfoS = new ROMETreeInfo;
   pEventHeader = new MEGEventHeader;
   rec->SetBranchStatus("*", 0); 
   rec->SetBranchStatus("Info.", 1); 
   rec->SetBranchStatus("Info.*", 1); 
   rec->SetBranchStatus("eventheader.", 1); 
   rec->SetBranchStatus("eventheader.mask", 1);
   rec->SetBranchStatus("spxclusters", 1); 
   rec->SetBranchStatus("spxclusters.*", 1); 
   rec->SetBranchStatus("spxtracks", 1); 
   rec->SetBranchStatus("spxtracks.*", 1); 
   rec->SetBranchStatus("dchspxmatchedtrack", 1); 
   rec->SetBranchStatus("dchspxmatchedtrack.*", 1); 
   rec->SetBranchStatus("spxhits", 1); 
   rec->SetBranchStatus("spxhits.*", 1);
   rec->SetBranchStatus("dchhits", 1); 
   rec->SetBranchStatus("dchhits.*", 1);
   rec->SetBranchStatus("dchtracks", 1);
   rec->SetBranchStatus("dchtracks.*", 1);
   
   rec->SetBranchAddress("dchtrackcandidates", &pDCHTrackCandArray, &bDCHTrackCand);
   rec->SetBranchAddress("dchtracks", &pDCHTrackArray, &bDCHTrack);
   rec->SetBranchAddress("Info.", &pInfoS, &bInfoS);
   rec->SetBranchAddress("eventheader.", &pEventHeader, &bEventHeader);
   rec->SetBranchAddress("spxclusters", &pSPXClusterArray, &bSPXCluster);
   rec->SetBranchAddress("spxtracks", &pSPXTrackArray, &bSPXTrack);
   rec->SetBranchAddress("spxhits", &pSPXHitArray, &bSPXHit);
   rec->SetBranchAddress("dchhits", &pDCHHitArray, &bDCHHit);
   rec->SetBranchAddress("dchspxmatchedtrack", &pMatchedTrackArray, &bMatchedTrack);

   Double_t targetSlantAngle = pTargetRunHeader->GetEulerAnglesAt(1);
   TVector3 targetOffset(pTargetRunHeader->GetTargetPosition());
   Double_t targetSemiAxis[2] = {pTargetRunHeader->GetLongSemiAxis() / 2. 
                                 * TMath::Cos(targetSlantAngle*TMath::DegToRad()),
                                 pTargetRunHeader->GetShortSemiAxis() / 2. };

   if (nEvent == 0) nEvent = rec->GetEntries(); 
   cout<<"Number of Event "<<nEvent<<endl;


   // Setup MEGPhysicsSelection
   MEGPhysicsSelection selector(kFALSE, 0, kTRUE);
   selector.SetThresholds(EBeamPeriodID::kMCPositronAnalysis);
   // update threshoulds and parameters to custom values
   selector.fTargetZOffset = targetOffset[2];
   selector.fTargetYOffset = targetOffset[1];
   selector.fTargetZ = targetSemiAxis[0];
   selector.fTargetY = targetSemiAxis[1];
   selector.fEPositronUncert = 100 * MeV;
   selector.fPhiPositronUncert = 200 * milliradian;
   selector.fThetaPositronUncert = 200 * milliradian;
   selector.fZPositronUncert = 20 * centimeter;
   selector.fYPositronUncert = 20 * centimeter;
   selector.fTargetFiducialLimit = 2.;

   for (Int_t iEvent = 0; iEvent < nEvent; iEvent++) {
      
      if (iEvent % 500 == 1) {
         cout<<iEvent<<" events finished..."<<endl;
      } 
      rec->GetEntry(iEvent);
      sim->GetEntry(iEvent);
      
      /////////////////////////////
      //                         //
      //   Acceptance cut        //
      //                         // 
      /////////////////////////////
      
      Int_t nsev = pMCMix->Getnsev();
      Int_t signalSevID = -1;//This is the id for the focusing particle in MCMix folder.
      Double_t signalTimeOffset(0);
      for (Int_t isev = 0; isev < nsev; isev++) {
         if (pMCMix->GetPicoSecAt(isev) == 0.) {
            signalSevID = pMCMix->GetSevIDAt(isev);
            signalTimeOffset = pMCMix->GetOffsetAt(isev);
         }
      }
      Int_t signalKineID = -1;//This is the id for the focusing particle in MCKine folder.
      for (Int_t i = 0; i < pMCKine->Getnprimary(); i++) {
         if (pMCKine->GetsevidAt(i) == signalSevID) {
            signalKineID = i;
         }
      }
      Double_t xvtx = pMCKine->GetxvtxAt(signalKineID);
      Double_t yvtx = pMCKine->GetyvtxAt(signalKineID);
      Double_t zvtx = pMCKine->GetzvtxAt(signalKineID);
      Double_t xmom = pMCKine->GetxmomAt(signalKineID);
      Double_t ymom = pMCKine->GetymomAt(signalKineID);
      Double_t zmom = pMCKine->GetzmomAt(signalKineID);        
      Double_t theta = TMath::ATan2(TMath::Sqrt(xmom*xmom + ymom*ymom), zmom);
      Double_t costheta = TMath::Cos(theta);
      Double_t phi = TMath::ATan2(ymom, xmom);
      
      // XEC acceptance
      // use MC truth for the gamma direction
      TVector3 gammaDir(-xmom, -ymom, -zmom);
      gammaDir = gammaDir.Unit();
      TVector3 vertexPos(xvtx, yvtx, zvtx);
      selector.XECAcceptance(gammaDir, vertexPos);

      if (selector.XECAcceptance(gammaDir, vertexPos)) {
         nGammaAcceptance++;     // count acceptance of gamma side
      } else {
         continue; // go next event
      }
      Double_t gammaTime = kSPXDRSTiming + signalTimeOffset;

      // SPX acceptance
      Bool_t bSPXHitFlg = false;
      Short_t nSPXAllHit = pSPXHitArray->GetEntriesFast();
      for (Short_t iSPXAllHit = 0; iSPXAllHit < nSPXAllHit; iSPXAllHit++) {
         pSPXHit = (MEGSPXHit*)pSPXHitArray->At(iSPXAllHit);
         Int_t mchitindex = pSPXHit->Getmchitindex();
         if (mchitindex >= 0) {
            pSPXMCHit = (MEGSPXMCHit*)pMCSPX->GetSPXMCHitAt(mchitindex);
            if (pSPXMCHit->Getsevid() == signalSevID) bSPXHitFlg = true;
         }
      }
      if (bSPXHitFlg == true) {
         nEffTCHits++;
      } else {
         continue; // go to next event
      }


      // Check signal track
      vector<Bool_t> signalDCHTrackFlag;
      TruthCheck(signalDCHTrackFlag,signalSevID, pDCHTrackArray, pDCHHitArray,
                 pMCDCH, pSPXClusterArray, pSPXHitArray, pMCSPX);
      Int_t nTrueTracks = std::count(signalDCHTrackFlag.begin(), signalDCHTrackFlag.end(), kTRUE);
      if (nTrueTracks) {
         nTrueSignal++;
      }



      vector<Bool_t> selected, selectedTruth;

  

      // Track selection
      selector.fSelectAPositron = kFALSE;
      selector.TrackSelection(selected, pDCHTrackArray, pSPXTrackArray, 
                              true, false, gammaTime, gammaDir);
      Int_t nSelectedTracks = std::count(selected.begin(), selected.end(), kTRUE);
      selectedTruth = AND(selected, signalDCHTrackFlag);
      Int_t nSelectedTracksTruth = std::count(selectedTruth.begin(), selectedTruth.end(), kTRUE);
      if (nSelectedTracks <= 0) {
         continue;
      }

      nEffEvents++;

 
      Double_t yvtx_rec(999.), zvtx_rec(999.); 
      Double_t theta_rec(999.), phi_rec(999.), mom_rec(999.);
      Double_t time_rec = -999.;
      for (size_t iDCHTrack = 0; iDCHTrack < selected.size(); iDCHTrack++) {
         if (!selected[iDCHTrack]) {
            continue;
         }
 
         pDCHTrack = (MEGDCHTrack*)pDCHTrackArray->At(iDCHTrack);
         Int_t iSPXTrack = pDCHTrack->GetbestSPXTrackIndex();
         if (iSPXTrack < 0) {
            continue;
         }
         pSPXTrack = (MEGSPXTrack*)pSPXTrackArray->At(iSPXTrack);
         time_rec = pSPXTrack->GettimeAtTarget();
         time_rec -= gammaTime;

         auto SelectedStateTarget = pDCHTrack->GetStateVectorTarget();
         yvtx_rec = SelectedStateTarget->GetY();
         zvtx_rec = SelectedStateTarget->GetZ();
         theta_rec = SelectedStateTarget->GetTheta();
         phi_rec = SelectedStateTarget->GetPhi();
         mom_rec = SelectedStateTarget->GetP();
      
         auto pMatchedTrack = (MEGDCHSPXMatchedTrack*)pMatchedTrackArray->At(iDCHTrack);
         Int_t clusterIndexInMatchedTrack = pSPXTrack->GetclusterIndexInMatchedTrack();
         

         auto stateFirst = pDCHTrack->GetStateVectorFirst();
         auto stateLast = pDCHTrack->GetStateVectorLast();
         //auto stateMatched = pSPXTrack->GetmatchedState();

         ngoodhits_ = pDCHTrack->Getngoodhits();
         nSPXHits_ = pSPXTrack->Getngoodhits();
         redChi2_ = pDCHTrack->Getchi2() / pDCHTrack->Getdof();
         YPositronUncert_ = TMath::Sqrt(SelectedStateTarget->GetCovarianceMEGAt(1,1));
         ZPositronUncert_ = TMath::Sqrt(SelectedStateTarget->GetCovarianceMEGAt(2,2));
         EPositronUncert_ = TMath::Sqrt(SelectedStateTarget->GetCovarianceMEGAt(3,3));
         ThetaPositronUncert_ = TMath::Sqrt(SelectedStateTarget->GetCovarianceMEGAt(4,4)) * radian;
         PhiPositronUncert_ = TMath::Sqrt(SelectedStateTarget->GetCovarianceMEGAt(5,5)) * radian;
         extrapolationLengthTarget_ = (stateFirst->GetT0() - SelectedStateTarget->GetT0()) * c_light;
         extrapolationLengthSPX_ = pMatchedTrack->GetpropagationLengthAt(clusterIndexInMatchedTrack);
         matchingDT_ = pMatchedTrack->GetTimeDifferenceAt(clusterIndexInMatchedTrack);
         matchingDW_ = pMatchedTrack->GetdwAt(clusterIndexInMatchedTrack);
         matchingDV_ = pMatchedTrack->GetdvAt(clusterIndexInMatchedTrack);
         matchingChi2_ = pMatchedTrack->Getchi2At(clusterIndexInMatchedTrack);
         
         Double_t eMomDiff = TMath::Abs(mom_rec - 52.828 * MeV); 
         Double_t ePhiDiff = TMath::Abs((phi_rec - phi) * 1.e3); //mrad
         Double_t eThetaDiff = TMath::Abs((theta_rec - theta) * 1.e3); //mrad
         Double_t eVerZDiff = TMath::Abs(zvtx_rec - zvtx);
         Double_t eVerYDiff = TMath::Abs(yvtx_rec - yvtx);
         Double_t targetTimeDiff = TMath::Abs(time_rec);
         good_ = kFALSE;
         if ( (Bool_t)((eMomDiff < cutRegion * mom_sigma * MeV) + !TailCutFlg[0]) 
              && (Bool_t)((ePhiDiff < cutRegion * phi_sigma) + !TailCutFlg[1])
              && (Bool_t)((eThetaDiff < cutRegion * theta_sigma) + !TailCutFlg[2])
              && (Bool_t)((eVerZDiff < cutRegion * verZ_sigma * millimeter) + !TailCutFlg[3])
              && (Bool_t)((eVerYDiff < cutRegion*verY_sigma * millimeter) + !TailCutFlg[4])
              && (Bool_t)((targetTimeDiff < cutRegion*t_sigma) + !TailCutFlg[5])
              ) {
            good_ = kTRUE;
         }


         // Define the target variable by the 'standardized distance' of the variables at target
         // Correlation among the variables are not considered.
         Double_t targetVar
               = eMomDiff * eMomDiff / (mom_sigma * MeV * mom_sigma * MeV)
               + ePhiDiff * ePhiDiff / (phi_sigma * phi_sigma)
               + eThetaDiff * eThetaDiff / (theta_sigma * theta_sigma)
               + eVerZDiff * eVerZDiff / (verZ_sigma * millimeter * verZ_sigma * millimeter)
               + eVerYDiff * eVerYDiff / (verY_sigma * millimeter * verY_sigma * millimeter)
               + targetTimeDiff * targetTimeDiff / (t_sigma * t_sigma);
         target_ = TMath::Sqrt(targetVar);
         
         outputTree->Fill();
         if (good_) {
            outputGoodTree->Fill();
         } else {
            outputBadTree->Fill();
         }

      }  
   }//End of event loop//

   cout<<nEffEvents<<" events used, "<<outputTree->GetEntries()<<" tracks used."<<endl;
   cout<<"Output: "<<output.GetName()<<endl;

   output.Write();
   output.Close();

}

//______________________________________________________________________________
Double_t FitHist(TH1* hist0, TH1* hist1) 
{

   hist0->Draw();
   hist1->SetFillStyle(1001);
   hist1->SetFillColor(kRed);
   hist1->Draw("same");

   Double_t mean = hist1->GetMean();
   Double_t sigma = hist1->GetStdDev();

   TF1 *gaus = new TF1("gaus", "gausn", mean - 10 * sigma, mean + 10 * sigma);
   gaus->SetLineColor(kCyan);
   gaus->SetLineWidth(2);
   gaus->SetParameters(hist1->Integral() / hist1->GetBinWidth(1),
                       mean, sigma);

   hist1->Fit(gaus, "", "", mean - 5*sigma, mean + 5*sigma);
   sigma = gaus->GetParameter(2);
   mean = gaus->GetParameter(1);
   Double_t norm = gaus->GetParameter(0);

   TF1 *dgaus = new TF1("dgaus", "gausn(0) + gausn(3)", mean - 10 * sigma, mean + 10 * sigma);
   dgaus->SetLineWidth(2);
   dgaus->SetLineColor(kAzure);
   dgaus->SetParameters(norm * 0.9, mean, sigma, norm * 0.1, mean, sigma * 2);
   hist0->Fit(dgaus, "L", "", mean - 10*sigma, mean + 10*sigma);
   hist0->Draw("same");

   return sigma; // 
}

//______________________________________________________________________________
void TruthCheck(vector<Bool_t> &signalDCHTrackFlag, Short_t signalSevID, TClonesArray* pDCHTrackArray, 
                TClonesArray* pDCHHitArray, MEGMCCYLDCHEvent* pMCDCH, TClonesArray* pSPXClusterArray,
                TClonesArray* pSPXHitArray, MEGMCSPXEvent* pMCSPX)
{
   // Determine the signal track as a track which contains hits of signal sub-event with
   // required purity.

    
   Int_t nDCHTracks = pDCHTrackArray->GetEntriesFast();
   signalDCHTrackFlag.assign(nDCHTracks, kFALSE);
   for (Int_t iDCHTrack = 0; iDCHTrack < nDCHTracks; iDCHTrack++) {
      MEGDCHTrack* pDCHTrack = (MEGDCHTrack*)pDCHTrackArray->At(iDCHTrack);
      Short_t nDCHHit = pDCHTrack->Getnhits();
      Short_t nDCHGoodHit = pDCHTrack->Getngoodhits();
      Short_t nDCHSignalHit = 0;
      for (Short_t iDCHHit = 0; iDCHHit < nDCHHit; iDCHHit++) {
         if (pDCHTrack->GetSkippedAt(iDCHHit)) {
            continue;
         }
         Short_t dchHitId = pDCHTrack->GethitindexAt(iDCHHit);
         MEGDCHHit* pDCHHit = (MEGDCHHit*)pDCHHitArray->At(dchHitId);
         Short_t dchMCHitId = pDCHHit->GetMCidx();
         if (dchMCHitId < 0) {
            continue;
         }
         MEGCYLDCHMCHit* pDCHMCHit = (MEGCYLDCHMCHit*)pMCDCH->GetCYLDCHMCHitAt(dchMCHitId);
         if (pDCHMCHit->Getsevid() == signalSevID) {
            nDCHSignalHit++;   
         }
      }      
      // if (Double_t(nDCHSignalHit)/Double_t(nDCHHit) >= FakeThreshold) {
      if (Double_t(nDCHSignalHit)/Double_t(nDCHGoodHit) >= FakeThreshold) {
         signalDCHTrackFlag[iDCHTrack] = kTRUE;
      }
   }

   // Bool_t bDCHAfterContaminationCutFlg=false;
   // Bool_t bSPXAfterContaminationCutFlg=false;
   // Int_t spxTrackIndex = pDCHTrack->GetbestSPXTrackIndex();
   // if (spxTrackIndex >= 0) {
   // }
   // MEGSPXCluster* pSPXCluster = (MEGSPXCluster*)pSPXClusterArray->At(SelectID[1]);
   // Short_t nSPXHitInCluster = pSPXCluster->Getnhits();
   // Short_t nSignalHitInCluster = 0;
   // for (Short_t iSPXHitInCluster = 0; iSPXHitInCluster < nSPXHitInCluster; iSPXHitInCluster++) {
   //    Short_t hitid = pSPXCluster->GethitidAt(iSPXHitInCluster);
   //    MEGSPXHit* pSPXHit = (MEGSPXHit*)pSPXHitArray->At(hitid);
   //    Int_t mchitindex = pSPXHit->Getmchitindex();
   //    if (mchitindex >= 0) {
   //       MEGSPXMCHit* pSPXMCHit = (MEGSPXMCHit*)pMCSPX->GetSPXMCHitAt(pSPXHit->Getmchitindex());            
   //       if (pSPXMCHit->Getsevid() == signalSevID) {
   //          nSignalHitInCluster++;
   //       }
   //    }
   // }
   // if (Double_t(nSignalHitInCluster)/Double_t(nSPXHitInCluster) >= FakeThreshold) bSPXAfterContaminationCutFlg = true;//signal SPX Cluster    
   // if(!bDCHAfterContaminationCutFlg||!bSPXAfterContaminationCutFlg){
   //    SelectID[0]=-1;
   //    SelectID[1]=-1;
   // }

}

