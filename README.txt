This directory contains:


* Matlab implementations of several temporally weighted linear predictive analysis methods listed below

WLP (Weighted Linear Prediction)
Ma, C., Kamp, Y. and Willems, L. F.,
Robust signal selection for linear prediction analysis of voiced speech, 
Speech Communication, 12(2):69-81, 1993
-generalizes conventional LP
-the first temporally weighted LP method where the error of each prediction is weighted to focus on certain parts of the signal frame
-default weighting scheme is a sequence of short-time energies with window length corresponding to the prediction order

SWLP (Stabilized WLP)
Magi, C., Pohjalainen, J., Bäckström, T. and Alku, P.,
Stabilised weighted linear prediction,
Speech Communication, 51(5):401-411, 2009
-guaranteed to produce a stable all-pole filter

XLP (eXtended weighted Linear Prediction)
Pohjalainen, J., Saeidi, R., Kinnunen, T. and Alku, P.,
Extended Weighted Linear Prediction (XLP) Analysis of Speech and its Application to Speaker Verification in Adverse Conditions,
Proc. Interspeech, Makuhari, Japan, September 2010
-generalizes LP and WLP
-the weighting is applied separately to each lag in the prediction of each sample, offering more flexibility in focusing the model on important data
-the default weighting scheme involves a matrix with time-averaged sums of absolute values of the predicted signal and the lag from which it is predicted

SXLP (Stabilized XLP)
Pohjalainen, J., Saeidi, R., Kinnunen, T. and Alku, P.,
Extended Weighted Linear Prediction (XLP) Analysis of Speech and its Application to Speaker Verification in Adverse Conditions,
Proc. Interspeech, Makuhari, Japan, September 2010

XLP-S (eXtended weighted Linear Prediction using autocorrelation Snapshot)
Pohjalainen, J. and Alku, P.,
Extended Weighted Linear Prediction Using the Autocorrelation Snapshot - A Robust Speech Analysis Method and its Application to Recognition of Vocal Emotions,
Proc. Interspeech, Lyon, France, August 2013
-generalizes LP, WLP and the previous constrained version of XLP
-two different weighting schemes from the paper are implemented

GMLP (Gaussian Mixture Linear Prediction)
Pohjalainen, J. and Alku, P.,
Gaussian Mixture Linear Prediction, Proc. ICASSP, Florence, Italy, May 2014
-a new, mixture autoregression based approach to linear predictive modeling


* Matlab implementation of multi-scale autoregressive filter postprocessing of feature vector sequences, as described in

Pohjalainen, J. and Alku, P.,
Multi-scale modulation filtering in automatic detection of emotions in telephone speech, Proc. ICASSP, Florence, Italy, May 2014