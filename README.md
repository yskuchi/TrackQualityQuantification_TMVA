# TrackQualityQuantification_TMVA
Quantify track quality with TMVA 

Author: Yusuke Uchiyama


### Environment ###
* ROOT 6.20.00
   * TMVA
   * ROOT interpreter (CLING)

## Quantification with classification methods ##
Use classifiers to quantify the track quality.
The classifiers are trained to separate 'good track' and 'bad track'.
Note that we do not intend to distinguish them but to quantify the quality of tracks by using the output of the classifiers.

The 'good' track is defined as a track satisfying the track selection and falling into the 5-sigma window 
(in EPositron, ThetaPositron, PhiPositron, YPositron, ZPositron, and TPositron).

The 'bad' track is defined as a track satisfying a loose track selection (at least successfully propagated to both target and SPX) but not falling into the 5-sigma window.

The output variable is used to compare and rank tracks.
It could be used to apply a cut or to categorize the tracks.

## Procedure ##
1. Make dataset files with `TrackQualityQuantificationDataset.cpp`. 
   Move `TrackQualityQuantificationInput.root` here.
2. Train models, evaluate the performance, tune hyper-parameters, and select a model 
   with `TrackQualityQuantificationClassifier.cpp`
   or `TrackQualityQuantificationCrossValidation.cpp`.
3. Apply the model to data. An example is `TrackQualityQuantificationApplication.cpp`. 
   See also "how to use in ROME analyzer".

```console
$ cd $MEG2SYS/analyzer
$ ./meganalyzer -I path/to/track_qualification_TMVA/TrackQualityQuantificationDataset.cpp
or
$ ./meganalyzer -I '../../track_qualification_TMVA/TrackQualityQuantificationDataset.cpp+'
```

```console
$ root  ./TrackQualityQuantificationClassification.cpp | tee TrackQualityQuantificationClassification.log
or
$ root './TrackQualityQuantificationClassification.cpp("BDT")' | tee TrackQualityQuantificationClassification.log
```

```console
$ root './TrackQualityQuantificationApplication.cpp("BDT")'
```


## Input variables ##

* `ngoodhits`: The number of good hits used in CDCH track fit.
* `redChi2`: The reduced chisquare of the CDCH track fit.
* `EPositronUncert`: The uncertainty on EPositron evaluated from the covariance matrix of the CDCH track fit (GeV).
* `PhiPositronUncert`: The uncertainty on PhiPositron evaluated from the covariance matrix of the CDCH track fit (deg).
* `ThetaPositronUncert`: The uncertainty on ThetaPositron evaluated from the covariance matrix of the CDCH track fit (deg).
* `extrapolationLengthTarget`: The extrapolation length from the first CDCH hit to the target (cm).
* `extrapolationLengthSPX`: The extrapolation length from the last CDCH hit to the SPX (matched point) (cm).
* `nSPXHits`: The number of good hits used in SPX track fit.
* `matchingDT`: The time difference between CDCH track and SPX track (s). 
* `matchingChi2`: The chisquare of CDCH-SPX matching.

![](fig/variables_id_c1.png "Input variable 1")
![](fig/variables_id_c2.png "Input variable 2")

## Analyze the trained data ##
When you train data (with  `TrackQualityQuantificationClassifier.cpp`
    or `TrackQualityQuantificationCrossValidation.cpp`)
you will get an output ROOT file (`TrackQualityQuantificationClassification.root` etc.),
which contains the dataset and results of the training.

You can use the GUI to see the output.
```
$ root -e 'TMVA::TMVAGui("./TrackQualityQuantificationClassification.root")'
or 
$ root
root[] TMVA::TMVAGui("./TrackQualityQuantificationClassification.root")
```

You can also access to the data with commands (or macro) like the following:
```
$ root ./TrackQualityQuantificationClassification.root
root [0]
Attaching file TrackQualityQuantificationClassification.root as _file0...
(TFile *) 0x364aec0
root [1] .ls
TFile**         TrackQualityQuantificationClassification.root
 TFile*         TrackQualityQuantificationClassification.root
  KEY: TDirectoryFile   dataset;1       dataset
root [2] dataset->cd()
(bool) true
root [3] .ls
TDirectoryFile*         dataset dataset
 KEY: TH2F      CorrelationMatrixS;1    Correlation Matrix (signal)
 KEY: TH2F      CorrelationMatrixB;1    Correlation Matrix (background)
 KEY: TDirectoryFile    InputVariables_Id;1     InputVariables_Id
 KEY: TDirectoryFile    Method_BDT;1    Directory for all BDT methods
 KEY: TTree     TestTree;1      TestTree
 KEY: TTree     TrainTree;1     TrainTree
root [4] CorrelationMatrixS->Draw("colz")

root [5] TrainTree->Print()

root [6] TrainTree->Draw("BDT", "classID==0")
Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1
(long long) 17110
root [7] TrainTree->Draw("BDT", "classID==1", "same")
(long long) 4169
root [8] CorrelationMatrixS->Draw("colz")
```

## How to use in ROME analyzer ##

To use TMVA in ROME analyzer, add the library in XML file
```xml
<UnixLibraries>
  <Library>
    <LibraryName>TMVA</LibraryName>
  </Library>
</UnixLibraries>
```

The application of the trained model is implemented in `analyzer/src/glb/MEGPhysicsSelection.cpp`.
See this code for an example of usage. The following is specific to the usage in this class in meg2/analyzer.

Move the weights files to meg2lib under an appropriate subdirectory.
```
#!console
$ rsync -auvv --chmod Dgo+w,Fgo+w to_meg2lib/20200800 meg-l-01.psi.ch:/meg/home/meg/offline/meg2lib/analyzer/cyldch/track_classification/
```
Then, set the directory, name, and method:
```
   MEGPhysicsSelection *selector = new MEGPhysicsSelection(kFALSE, 0, kTRUE);
   selector->SetThresholds(EBeamPeriod::kMCPositronAnalysis);
   selector->fTrackQualityQuantifier.fName  = "TrackQualityQuantificationClassification";
   selector->fTrackQualityQuantifier.fDirectory = "$(MEG2SYS)/meg2lib/analyzer/cyldch/track_classification/20200800/";
   selector->fTrackQualityQuantifier.fMVAMethod = "BDT";

```


## Notes ##

### Aug. 2020 ###
Use MC, signal e+ mixed to 7e7 Michel e+, with a sample of 40k events for training and testing.
17282 'good track' samples and 4211 'bad track' samples.

![](fig/rejBvsS.png "ROC curves")
![](fig/overtrain_BDT.png "BDT output")  
This is trained and tested with 50:50 samples with AdaBoost BDT.


```
LD                       : Results for LD coefficients:
                         : -----------------------------------------------------
                         :                  Variable:               Coefficient:
                         : -----------------------------------------------------
                         :                 ngoodhits:                    +0.004
                         :                   redChi2:                    -0.124
                         :           EPositronUncert:                  -443.178
                         :       ThetaPositronUncert:                    +0.117
                         :         PhiPositronUncert:                    -0.029
                         : extrapolationLengthTarget:                    -0.008
                         :    extrapolationLengthSPX:                    -0.003
                         :                  nSPXHits:                    +0.006
                         :              matchingChi2:                    -0.041
                         :                matchingDT:             -15841643.701
                         :                  (offset):                    +0.112
                         : -----------------------------------------------------
```

```
BDT                      : Ranking result (top variable is best ranked)
                         : -----------------------------------------------------------
                         : Rank : Variable                  : Variable Importance
                         : -----------------------------------------------------------
                         :    1 : ngoodhits                 : 1.853e-01
                         :    2 : nSPXHits                  : 1.350e-01
                         :    3 : matchingChi2              : 1.172e-01
                         :    4 : redChi2                   : 9.971e-02
                         :    5 : extrapolationLengthSPX    : 9.819e-02
                         :    6 : extrapolationLengthTarget : 9.424e-02
                         :    7 : PhiPositronUncert         : 7.512e-02
                         :    8 : matchingDT                : 6.881e-02
                         :    9 : ThetaPositronUncert       : 6.695e-02
                         :   10 : EPositronUncert           : 5.957e-02
                         : -----------------------------------------------------------
```
```
                         : Evaluation results ranked by best signal efficiency and purity (area)
                         : -------------------------------------------------------------------------------------------------------------------
                         : DataSet       MVA
                         : Name:         Method:          ROC-integ
                         : dataset       BDT            : 0.861
                         : dataset       BDTG           : 0.859
                         : dataset       PDEFoamBoost   : 0.851
                         : dataset       BDTD           : 0.847
                         : dataset       MLPBNN         : 0.828
                         : dataset       LD             : 0.827
                         : dataset       SVM            : 0.801
                         : -------------------------------------------------------------------------------------------------------------------
                         :
                         : Testing efficiency compared to training efficiency (overtraining check)
                         : -------------------------------------------------------------------------------------------------------------------
                         : DataSet              MVA              Signal efficiency: from test sample (from training sample)
                         : Name:                Method:          @B=0.01             @B=0.10            @B=0.30
                         : -------------------------------------------------------------------------------------------------------------------
                         : dataset              BDT            : 0.080 (0.220)       0.509 (0.686)      0.882 (0.918)
                         : dataset              BDTG           : 0.063 (0.215)       0.507 (0.666)      0.880 (0.909)
                         : dataset              PDEFoamBoost   : 0.062 (0.098)       0.506 (0.556)      0.867 (0.879)
                         : dataset              BDTD           : 0.072 (0.200)       0.449 (0.616)      0.858 (0.894)
                         : dataset              MLPBNN         : 0.058 (0.052)       0.408 (0.424)      0.831 (0.853)
                         : dataset              LD             : 0.050 (0.038)       0.386 (0.428)      0.839 (0.847)
                         : dataset              SVM            : 0.017 (0.031)       0.301 (0.321)      0.763 (0.794)
                         : -------------------------------------------------------------------------------------------------------------------
```

Cross validation

![](fig/rejBvsS_BDT_crossvalidation.png "Cross validation")
```
Summary for method BDT
        Fold 0: ROC int: 0.867373, BkgEff@SigEff=0.3: 0.879
        Fold 1: ROC int: 0.873474, BkgEff@SigEff=0.3: 0.903
        Fold 2: ROC int: 0.857859, BkgEff@SigEff=0.3: 0.887
        Fold 3: ROC int: 0.866336, BkgEff@SigEff=0.3: 0.887
```
We can check the robustness of model with the cross validation.
e.g if the k ROC-curves from the k-fold cross validation show large variation, then the model
is not stable; you should choose other models.

Tried hyper-parameter tuning with TMVA::HyperParameterOptimisation class and method->OptimizeTuningParameters()
but none of them worked well.
