:Begin:
:Function: FrequencyModifiedFourierTransform
:Pattern: FrequencyModifiedFourierTransform[n_Integer,data_?(VectorQ[#, (MatchQ[#, {Repeated[_?NumberQ, {3}]}] &)]&) ]
:Arguments: {n,data}
:ArgumentTypes: {Integer,Manual}
:ReturnType: Manual
:End:

:Evaluate: FrequencyModifiedFourierTransform::usage="FrequencyModifiedFourierTransform[n,{{t,x,y},...}] returns n frequencies from the FMFT of the data given"
:Evaluate: FrequencyModifiedFourierTransform::OddNFreq="`` is odd.  The FMFT is currently set to require an even number of frequencies only"
:Evaluate: FrequencyModifiedFourierTransform::FMFTFailure="The FMFT routine returned an error value"
:Evaluate: FrequencyModifiedFourierTransform::BadDimensions="The FMFT routine recieved an input with improper dimensions"
:Evaluate: FrequencyModifiedFourierTransform::mlink="MathLink failed with error message: ``"
