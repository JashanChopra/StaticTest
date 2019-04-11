% staticUncertainty is a script to analyze confidence in test stand data
% Jaret Anderson 106331457 for ASEN 2004 Lab 2
% uses classDeltaV.mat

% setup
A = load('classDeltaV.mat'); % load deltaV data from class dataset
deltaV = A.deltaV;
sigmaDeltaV = std(deltaV); % calculate standard deviation in dataset
meanDeltaV = mean(deltaV);
SEMclass = sigmaDeltaV / sqrt(length(deltaV)); % calculate SEM

% confidence intervals 
z = [1.96; 2.24; 2.58]; % scalar multipliers (95, 97.5, 99)
CIclass = [meanDeltaV - z*SEMclass meanDeltaV + z*SEMclass]; % Upper bound Lower bound