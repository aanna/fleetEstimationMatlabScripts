close all; clear; clc;
%% Cost for traveling between stations
% cost of rebalancing
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
stationX = stationsData(:,2);
stationY = stationsData(:,3);

dist_factor = 1.37; % adjustment for euclidean dist, based on Spieser2014
Dist_matrix = dist_factor * pdist2([stationX(:), stationY(:)], [stationX(:), stationY(:)]);

%% Save to file
fileTOSave = sprintf('distance%dStations.txt', length(stationX));
delimiter = ' ';
dlmwrite(fileTOSave, Dist_matrix, delimiter);

%% Cost based on the travel time
ave_tt = 32; % km/h
ave_tt_ms = ave_tt * 1000/3600;
tt_matrix = zeros(length(Dist_matrix));
for i = 1: length(Dist_matrix)
    for j = 1: length(Dist_matrix)
        dist_temp = Dist_matrix(i,j);
        tt_matrix(i,j) = round(Dist_matrix(i,j) / ave_tt_ms);
    end
end

%% Save to file
fileTime = sprintf('RebTime%dStations.txt', length(Dist_matrix));
delimiter = ' ';
dlmwrite(fileTime, tt_matrix, delimiter);
