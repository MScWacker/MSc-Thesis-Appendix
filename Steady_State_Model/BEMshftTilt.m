function [AeroPower, AeroQ, AeroTx] = BEMshftTilt(n_b, l_b, l_c, pitch, beta, r_hub, h_hub, omega, v_w, psi)
%Rotor module using AeroDyn with wind shear and AeroDyn initialisation
% Instead of changing the yaw, the shftTilt could be changed to 'num2str(rad2deg(beta)) '
beta = rad2deg(beta);
rpm = round((omega*60)/(2*pi), 4);
%theta_0 is the starting azimuth postion of the rotor between 0 and 2pi 
%find time to move to azimuth and set time step

Tmax = 2*pi/n_b/omega;
dt = Tmax/250*n_b;   
%Tmax = dt; % uncomment if only one timestep should be looked at

% AeroDyn manual recommends at least 200 time steps per revolution
% AeroDyn will not run if time step is greater that simulation time

n_bn = 15; % number of blade nodes


%% Location of AeroDyn
AeroDyn = 'C:\Users\Jannis\AeroDyn\bin\AeroDyn_Driver_x64.exe ';
% Other files used are 
% 'ad_primary_Daisy.inp'
% 'ad_blade_Daisy.inp',
% 'ad_aerfoil_S809.inp'
% output file is named 
% 'AeroDynResults.1.out'

%% Generate and populate the AeroDyn Driver File
% create the AeroDyn driver file
DriverFile = 'C:\Users\Jannis\AeroDyn\ad_driver.inp';

% open the AeroDyn driver file
DriverFileID = fopen(DriverFile, "w");

% populate the AeroDyn driver file
fprintf(DriverFileID, [...
    '-----  AeroDyn Driver v15 Input File  -------------------------------------- \r\n'...
    'UAE Phase 3 turbine (BEM) \r\n' ...
    '-----  Input Configuration  ---------------------------------------------------- \r\n'...
    '       False               Echo            -  Echo input parameters to "<rootname>.ech"? \r\n'...
    '"ad_primary_Daisy.inp"     AD_InputFile    -  Name of the primary AeroDyn input file \r\n'...
    '-----  Turbine Data  -----------------------------------------------------------\r\n'...
    '      ' num2str(n_b) ' \t NumBlades \t - Number of blades (-) \r\n' ...
    '      ' num2str(r_hub) ' \t HubRad \t - Hub radius (m) \r\n' ...
    '      ' num2str(h_hub) ' \t HubHt \t - Hub height (m) \r\n' ...
    '      0 \t Overhang \t - Overhang (m) \r\n' ...
    '      ' num2str(beta) '   ShftTilt        - Shaft tilt (deg) \r\n' ...
    '      0                    Precone         - Blade precone (deg) \r\n' ...
    '-----  I/O Settings  ----------------------------------------------------------- \r\n' ...
    '       "AeroDynResults"       OutFileRoot     - Root name for any output files (use "" for .dvr rootname) (-) \r\n' ...
    '       True                TabDel          - When generating formatted output (OutForm=True), make output tab-delimited (fixed-width otherwise) (flag) \r\n' ...
    '       "ES10.3E2"          OutFmt          - Format used for text tabular output, excluding the time channel.  Resulting field should be 10 characters. (quoted string) \r\n' ...
    '       True                Beep            - Beep on exit (flag) \r\n' ...
    '-----  Combined-Case Analysis  ------------------------------------------------- \r\n' ...
    '     1      NumCases        - Number of cases to run \r\n']);

% Input the Case Variables
fprintf(DriverFileID,'%8s %8s %8s %8s %8s %8s %8s\r\n','WndSpeed','ShearExp','RotSpd','Pitch','Yaw','dT','Tmax');
fprintf(DriverFileID,'%8s %8s %8s %8s %8s %8s %8s\r\n','(m/s)','(-)','(rpm)','(deg)','(deg)','(s)','(s)');
fprintf(DriverFileID,'%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g \r\n',[v_w psi rpm pitch 0 dt Tmax]);

fclose(DriverFileID);

%% Generate and populate Blade Input file
BladeFile = 'C:\Users\Jannis\AeroDyn\ad_blade_Daisy.inp';
BladeFileID = fopen(BladeFile, "w"); 
BladeProp = zeros(7, n_bn); 
BladeProp(1,:) = 0:l_b/(n_bn-1):l_b;    %BlSpn
BladeProp(2,:) = zeros(n_bn, 1);    %BlCrvAC
BladeProp(3,:) = zeros(n_bn, 1);    %BlSwpAC
BladeProp(4,:) = zeros(n_bn, 1);    %BlCrvAng
BladeProp(5,:) = zeros(n_bn, 1);    %BlTwist
BladeProp(6,:) = ones(n_bn, 1)*l_c; %BlChord
BladeProp(7,:) = ones(n_bn, 1);     %BlAFID

fprintf(BladeFileID, [...
    '------- AERODYN v15.04.* BLADE DEFINITION INPUT FILE -------------------------------------\r\n'...
    'This file describes the blade shape and is changed according to l_b and l_c \r\n'...
    '======  Blade Properties =================================================================\r\n'...
    '\t' num2str(n_bn) '  NumBlNds    - Number of blade nodes used in the analysis (-)\r\n']);
   
fprintf(BladeFileID,'%8s %8s %8s %8s %8s %8s %8s\r\n','BlSpn','BlCrvAC','BlSwpAC','BlCrvAng','BlTwist','BlChord','BlAFID');
fprintf(BladeFileID,'%8s %8s %8s %8s %8s %8s %8s\r\n','(m)','(m)','(m)','(deg)','(deg)','(m)','(-)');
fprintf(BladeFileID,'%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g \r\n',BladeProp);

fclose(BladeFileID);

%% Run AeroDyn
Input = [AeroDyn, DriverFile];
system(Input);

%% Define outputs from AeroDyn Data
output = 'AeroDynResults.1.out';
Data = importdata(output, '\t');
%Time, Omega, Power, Thrust, Torque, B1Azimuth, Cp, Cq, Ct, RtArea, TSR, 

%AeroTSR = mean(Data.data(:, 3));
AeroPower = mean(Data.data(:, 4));
AeroTx = mean(Data.data(:, 5));
% AeroTy = mean(Data.data(:, 6));
% AeroTz = mean(Data.data(:, 7));
AeroQ = mean(Data.data(:, 8));
% AeroCp = mean(Data.data(:, 9))
% AeroCq = mean(Data.data(:, 10));
% AeroCt = mean(Data.data(:, 11));


end



