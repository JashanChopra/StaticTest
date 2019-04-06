%% ASEN 2004 - Rocket Bottle Lab
%{

    Authors: Jashan Chopra (107689146)
    Date Created: April 4th, 2019

Script Purposes and goals:
  1) Read static test stand data
  2) edit data to appropriate ranges
  3) automate the data curtailing process
  4) output data to matrix for statistical usage

%}

% General Housekeeping
clc; clear; close all;

%% Import the test stand data, define initial conditions
filenames = ["LA8am_test3","LA8am_test4"];  % Create a vector of all filenames
numFiles = size(filenames);                 % Number of files

g = 9.80665;        % gravitational constant
mProp = 1;          % [kg] [1 L = 1000g = 1kg]
mBottle=.15;        %[kg] % Mass of Empty Bottle
Patm=83427;         %[Pa] % Atmospheric Pressure
Pgage=345738;       %[Pa] % Intial Gauge Pressure of air in bottle

volBottle=.002;             %[m^3] % Empty bottle Volume
volWater=.001;              %[m^3] % Initial Volume of Water
volAir=volBottle-volWater;  % Volume of air in Bottle Initial

R=287;      % Universal Gas Constant Air
Tair=300;   %[K] % Initital Temp of Air

mAir=((Pgage+Patm)*volAir)/(R*Tair);    % Mass of Air Initital
mAirFinal=((Patm)*volAir)/(R*Tair);     % Mass of Air Final

%% Iteration of deltaV calculation
deltaV = zeros(numFiles,1);     % preallocate deltaV
for f = filenames               % loop iterates over all files

    % automated data trim
    data = fileLoad(f);                                                 % Load data
    indicies = find(data1 <= 0);                                        % negative indices
    data1(indicies) = [];                                               % remove negative indices

    % calculations
    time = (1 / frequency) * linspace(0,length(data),length(data))';    % time vector
    fitobject = fit(time,data,'cubicinterp');                           % curve fit
    area = integrate(fitobject,time(end),time(1));                      % integration
    isp = area / (mProp*g);                                             % isp calc

    mInitial = mBottle + mProp + mAir;                                  % initial mass
    mFinal = mInitial - mProp - mAir + mAirFinal;                       % final mass

    deltaV(f,:) = isp*g*log(mInitial / mFinal);                         % [m/s] % DeltaV

end
