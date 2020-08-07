# track_qualification_TMVA
Track qualification with TMVA 

Author: Yusuke Uchiyama


### Environment ###
* ROOT 6.20.00
   * TMVA
   * ROOT interpreter (CLING)

## Procedure ##
1. Make dataset files with `TrackQualificationDataset.cpp`. Move `TrackQualificationInput.root` here.
2. Train models, evaluate the performance, and select a model with 'TrackQualificationClassifier.cpp'
    or 'TrackQualificationCrossValidation.cpp.
3. Apply the model to data

```console
$ cd $MEG2SYS/analyzer
$ ./meganalyzer -I path/to/track_qualification_TMVA/TrackQualificationDataset.cpp
or
$ ./meganalyzer -I '../../track_qualification_TMVA/TrackQualificationDataset.cpp+'
```

```console
$ root  ./TrackQualificationClassification.cpp | tee TrackQualificationClassification.log
or
$ root './TrackQualificationClassification.cpp("BDT")' | tee TrackQualificationClassification.log
```

### Input variables ###

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


