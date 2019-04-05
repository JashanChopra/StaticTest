function [monteMat] = monte(numTrials,vw)
%% ASEN 2004 - Rocket Bottle Lab - Monte Carlo Simulation Function
%{

    Authors: Jashan Chopra (107689146)
    Date Created: April 3rd, 2019

Script Purposes and goals:
    1) Input required number of trials
    2) Randomize initial conditions sensitive to uncertainity in testing
    3) Run thermodynamic and rocket equation models over number of trials
    4) Output final matrix of monte carlo results

Outputs: [ FOR EACH TRIAL - AND SEPERATED BY TESTING METHOD ]
    1) Maximum downrange distance (x)
    2) Maximum crossrange distance (y)
    3) Maximum height (z)

%}
tic

%% BASE FOR LOOP CONTROLS MONTE CARLO SIMULATION

% Preallocate final matrix
monteMat = zeros(numTrials,6);

for i = 1:numTrials

    %% Identify Constants from each model

        % Given Constants
        g=9.81;             %[m/s] % Gravity Constant
        Cd=.8;              % Dischart Coefficient
        rhoAtm=.961;        %[kg/m3] % Ambient Air Density
        rhoWater=1000;      %[kg/m3] % Density of Water
        volBottle=.002;     %[m^3] % Empty bottle Volume
        Patm=83427;         %[Pa]   % Atmospheric Pressure
        gamma=1.4;          % Ratio of Specific Heats for Air
        Dthroat=.021;       %[m] % Diameter of Throat
        Dbottle=.105;       %[m] % Diameter of Bottle
        R=287;              % Universal Gas Constant Air
        mBottle=.15;        %[kg] % Mass of Empty Bottle
        CDrag=.5;           % Drag Coefficient
        Pgage=345738;       %[Pa] % Intial Gauge Pressure of air in bottle
        volWater=.001;      %[m^3] % Initial Volume of Water
        Tair=300;           %[K] % Initital Temp of Air
        v0=0;               % [m/s] Initial Velocity
        theta=pi/4;         %[rad] % Initial Angle
        ls=.5;              %[m] % Length of Test Stand
        mProp = 1;          % [kg] [1 L = 1000g = 1kg]

        x0=0;   %[m] % Initial Horizontal Distance
        z0=.25; %[m] % Initial Vertical Distance
        y0 = 0; %[m] % Initial Lateral Distance

        % Calculated Constants

        volAir=volBottle-volWater;              % Volume of air in Bottle Initial
        mAir=((Pgage+Patm)*volAir)/(R*Tair);    % Mass of Air Initital
        Athroat=pi*(Dthroat/2)^2;               % Area of Throat
        Abottle=pi*(Dbottle/2)^2;               % Area of Bottle
        mAirFinal=((Patm)*volAir)/(R*Tair);     % Mass of Air Final
        mWater = rhoWater * volWater;           % Initial mass of water

        vx = 0; % Initial Velocity x
        vy = 0; % Initial Velocity y
        vz = 0; % Initial Velocity z

    %% Thermodynamic Model
    % Note: Thermodynamic & Isp model comments in their respsective files
    state1Vec = [volAir vx vz vy x0 z0 y0];
    tspan = 0:.00001:.5;
    optstage1 = odeset('Events',@opt1);
    [t, ds] = ode45(@(t,s) state1func(t,s,Cd,Athroat,gamma,...
    rhoAtm,Patm,Pgage,volAir,Abottle,mBottle,mAir,CDrag,volBottle,...
    rhoWater,vw(i,:)), tspan, state1Vec, optstage1);
    m2 = mAir;
    Pair = Pgage + Patm;
    state2Vec = [m2;ds(end,2);ds(end,3);ds(end,4);ds(end,5);ds(end,6);ds(end,7)];
    tspan2 = t(end):.000001:6;
    optstage2 = odeset('Events',@opt2); % This function modifies ode45 to stop
        [t2,ds2] = ode45(@(t2,s) state2func(t2,s,mAir,gamma,Pair,volBottle,R,Patm,...
        Athroat,Cd,Abottle,rhoAtm,volAir,CDrag,mBottle,vw(i,:)),tspan2, state2Vec,optstage2);

    % Remove values where z has gone below the surface of the Earth
    indices = find(ds2(:,6) < 0);
    ds2(indices,:) = [];

        % Concatenate two stages
        xpos=[ds(:,5);ds2(:,5)];
        zpos=[ds(:,6);ds2(:,6)];
        ypos=[ds(:,7);ds2(:,7)];

    % Put current max iterator values in final matrix
    monteMat(i,1) = max(xpos);
    monteMat(i,2) = min(ypos);
    monteMat(i,3) = max(zpos);

    %% Rocket ISP Model Datas
    % Note: Commented and formatted code is present in 'rocketMain_Isp.m'

    filename1 = 'LA8am_test3';
    data1 = load(filename1);
    data1 = data1(2215:3300,3); % Just the summed load data
    indicies = find(data1 <= 0);
    data1(indicies) = [];
    force1 = data1*4.44822; % Summed force converted to newtons
    frequency = 1.652 * 1000; % [Hz] Sampling Rate
    time = (1 / frequency) * linspace(0,length(data1),length(data1)); % time vector
    time = time'; % make it into a column vector
    isp = trapz(time,force1);
    isp = isp / (mProp*g);
    mInitial = mBottle + mProp + mAir; % Initial Mass
    mFinal = mInitial - mProp - mAir + mAirFinal; % Final Mass
    deltaV = isp*g*log(mInitial / mFinal); % [m/s] % DeltaV
    vx = deltaV*cos(theta); % Initial Velocity x
    vy = 0;                 % Initial Velocity y
    vz = deltaV*sin(theta); % Initial Velocity z
    x0 = 0;     %[m] % Initial Horizontal Distance
    z0 = .25;   %[m] % Initial Vertical Distance
    y0 = 0;     %[m] % Initial Lateral Distance
    stateVec = [vx;vz;vy;x0;z0;y0];
    tspan = 0:.00001:5; % Create a time span, 0 to 2 seconds
    optstage1 = odeset('Events',@opt1); % opt1.m outputs the cancel signal
    [t, ds] = ode45(@(t,s) state2funcROCKET(t,s,vw(i,:)), tspan, stateVec, optstage1);

    % Put current max iterator values in final matrix
    monteMat(i,4) = max(ds(:,4));
    monteMat(i,5) = min(ds(:,6));
    monteMat(i,6) = max(ds(:,5));
end
toc
end
