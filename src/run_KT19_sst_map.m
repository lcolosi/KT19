%% Sea Surface Temperature Map from MASS Flight
% Luke Colosi | lcolosi@ucsd.edu | October 17th, 2022

clc, clear, close all

% Change directory 
cd '/Users/lukecolosi/Desktop/projects/graduate_research/deployments/SMODE2022/src'

% Set paths for KT19, SPAN, and kml data directories and output mat files for SST
% and figures
dirinspan = '../data/SPAN/20221026_2/';
dirinkt19 = '../data/KT19/20221026_2/';
dirinkml = '../daa/kml/flight_lines/20221026_MASS_TO.kml';
dirinp = '../data/kml/SMODE_Polygon.kml';
dirout = '../data/KT19/processed/'; 

% Set other inputs to function
dt = 1;                                                                     % Averaging time period (units: seconds)
region = [-124.7, -123.8; 36.3, 37.1];                                      % Regional extent of flight lines 'region', region
c_range = [13, 18]; 
dT = 0.2; 

% Generate files and figures for SST map
[fileout] = KT19_sst_map(dirinspan,dirinkt19,dirout,'region', region, 'dirinp', dirinp, 'dt', dt, 'c_range', c_range, 'dT', dT);


 