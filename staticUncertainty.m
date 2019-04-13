% staticUncertainty is a script to analyze confidence in test stand data
% Jaret Anderson 106331457 for ASEN 2004 Lab 2
% uses classDeltaV.mat, classTime.mat, IspVec.mat, 
% peakThrust.mat, meanThrust.mat
close all; clear all; clc

%% file input
A = load('classDeltaV.mat');    % load deltaV data from class dataset
deltaV = A.deltaV;
B = load('classTime.mat');    % load time data from class dataset
timeVec = B.timeVec;
C = load('IspVec.mat');    % load time data from class dataset
IspVec = C.IspVec;
D = load('peakThrust.mat');    % load time data from class dataset
peakThrust = D.peakThrust;
E = load('meanThrust.mat');    % load time data from class dataset
meanThrust = E.meanThrust;

%% means and uncertainties
sigmaIsp = std(IspVec);      % calculate standard deviation in Isp dataset
meanIsp = mean(IspVec);
SEMisp = sigmaIsp / sqrt(numel(IspVec)); % calculate SEM for Isp

sigmaMeanThrust = std(meanThrust);
sigmaPeakThrust = std(peakThrust);  % uncertainty values
sigmaTime = std(timeVec);
meanMeanThrust = mean(meanThrust);
meanPeakThrust = mean(peakThrust);  % mean values
meanTime = mean(timeVec);

%% incremental plot of SEM vs N
SEMinc = zeros(1000,1);         % vector of SEM values to be plotted vs N
N = 1:1:1000;                   % number of trials vector

for i = N
   SEMinc(i) = sigmaIsp/sqrt(i); % calculate SEM for each number of trials
end

figure(1) 
plot(N(2:end),SEMinc(2:end))
title('SEM vs Number of Trials')
xlabel('Number of Trials')
ylabel('SEM')
hold on
plot(40,SEMinc(40),'ro')      % usable class data
plot(65,SEMinc(65),'bo')      % total class data
legend('Incremental SEM', 'Number of Usable Class Data Files (40)',...
    'Total Number of Class Data Files (65)')

%% plot histograms
figure(2)
hold on
subplot(2,2,1)% isp
histogram(IspVec)
title('Specific Impulse of Trials')
xlabel('Isp (s)')
ylabel('Number of Trials')

subplot(2,2,2)% time of thrust
histogram(timeVec)
title('Time of Thrust of Trials')
xlabel('Time (s)')
ylabel('Number of Trials')

subplot(2,2,3)% peak thrust
histogram(peakThrust)
title('Peak Thrust of Trials')
xlabel('Thrust (N)')
ylabel('Number of Trials')

subplot(2,2,4)% mean thrust
histogram(meanThrust)
title('Mean Thrust of Trials')
xlabel('Thrust (N)')
ylabel('Number of Trials')

%% confidence intervals 
z = [1.96; 2.24; 2.58]; % scalar multipliers (95, 97.5, 99)
CIclass = [meanIsp - z*SEMisp meanIsp + z*SEMisp]; 
                  % Upper bound         %Lower bound

% calculate number of trials needed for specified confidence interval
rangeTolerance = 0.1;                             % range for 95/97.5/99% confidence interval
Ntenth = (sigmaIsp.*z./rangeTolerance).^2;     % calculate number of trials for CI
Ntenth = ceil(Ntenth);                            % cannot do a partial trial

rangeTolerance = 0.01;                            % range for 95/97.5/99% confidence interval
Nhundredth = (sigmaIsp.*z./rangeTolerance).^2; % calculate number of trials for CI
Nhundredth = ceil(Nhundredth);                    % cannot do a partial trial