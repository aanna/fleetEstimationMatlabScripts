close all; clear; clc;
%% Cost for traveling between stations
% cost of rebalancing, travel time is represented in secs
% tt_matrix stores travel time for all trips between i and j
% travel time is in seconds
% distance is in meters
% velocity is in meters per second

%% Read file
% station file
mac = false;
if (~mac)
    facilityFile = sprintf('/home/kasia/Dropbox/matlab/2016-03-Demand_generation/facility_location/stations_ecbd34.txt');
else
    facilityFile = sprintf('/Users/katarzyna/Dropbox/matlab/2016-03-Demand_generation/facility_location/stations_ecbd34.txt');
end
stationsData = dlmread(facilityFile, ' ', 0, 0);

f_ids = stationsData(:,1);
% station coordinates are in utm coord system in meters
stationX = stationsData(:,2);
stationY = stationsData(:,3);

dist_factor = 1.37; % adjustment for euclidean dist, based on Spieser2014
% distance matrix is in meters
Dist_betweenSt_inM = dist_factor * pdist2([stationX(:), stationY(:)], [stationX(:), stationY(:)]);

%% Save to file
fileTOSave = sprintf('distance%dStations.txt', length(stationX));
delimiter = ' ';
dlmwrite(fileTOSave, Dist_betweenSt_inM, delimiter);

%% Cost based on the travel time
% travel time = distance travelled / speed
% ave speed for arterial roads in Singapore is 28.9 km/h (8.02 m/s)
speed_kmh = 28.9; % km/h
speed_ms = speed_kmh * 1000/3600;
tt_inSeconds = zeros(length(Dist_betweenSt_inM));
tt_inMinutes = zeros(length(Dist_betweenSt_inM));
for i = 1: length(Dist_betweenSt_inM)
    for j = 1: length(Dist_betweenSt_inM)
        % travel time in seconds
        tt_inSeconds (i,j) = ceil (Dist_betweenSt_inM(i,j) / speed_ms);
        % travel time in minutes
        tt_inMinutes (i,j) = ceil (tt_inSeconds (i,j)/60);
    end
end

%% Save to file
fileTime = sprintf('RebTimeInSecs%dStations.txt', length(Dist_betweenSt_inM));
delimiter = ' ';
dlmwrite(fileTime, tt_inSeconds, delimiter);
