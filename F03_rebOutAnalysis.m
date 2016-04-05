close all; clear; clc;

% have to provide as an input the number of stations and number of
% rebalancing intervals (it can't be recognized from the output file)
% new output file in format time, from_node, to_node, count

n_stations = 34;
n_reb_periods = 96;
rebalancing_interval = 24*60*60/n_reb_periods; % in seconds
vec_for_reb = 1:n_stations*n_stations;
reb_matrix = transpose(reshape(vec_for_reb, [n_stations, n_stations]));

%% Import output from Gurobi cpp
disp('1. Import output from Gurobi cpp...')
filename = '/home/kasia/Documents/rebalancing_cpp/rebalancingMethods/rebalancing_offline/solution_rebalancing.sol';
delimiter = ',';
startRow = 3;
formatSpec = '%s%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

var_name = dataArray{:, 1};
time_ = dataArray{:, 2};
depSt_orSt = dataArray{:, 3}; % depending on thr var_name, i.e., if reb_veh
%then this is departing station, if available_veh then @ which station
arrSt_orSt = dataArray{:, 4}; % same as above
quantity = dataArray{:, 5};

clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% Import stations
disp('2. Stations file...');
facilityFile = sprintf('stations_ecbd34.txt');

stationsData = dlmread(facilityFile, ' ', 0, 0);

f_ids = stationsData(:,1);
stationX = stationsData(:,2);
stationY = stationsData(:,3);

%% Generate 3 matrices
disp('2. Generate 3 matrices...')
% 1) vehicles available at each station for all time intervals
available_veh = (strcmp(var_name,'vhs_st_i'));
available_veh_m = zeros(n_reb_periods, n_stations);
% 2) rebalancing counts between stations for all time intervals
reb_veh = (strcmp(var_name,'nEmptyVhsTime'));
reb_veh_m = zeros(n_reb_periods, n_stations*n_stations);
% 3) moving vehicles for all time intervals
veh_in_transit = (strcmp(var_name,'in_transit'));
veh_in_transit_m = zeros(n_reb_periods, 1);

for i =1 : length (var_name)
    if (available_veh(i))
        available_veh_m(time_(i)+1, arrSt_orSt(i)+1) = quantity(i);
    end
    
    if (reb_veh(i))
        column_indx = reb_matrix(depSt_orSt(i)+1, arrSt_orSt(i)+1);
        reb_veh_m(time_(i)+1, column_indx) = quantity(i);
    end
    
    if (veh_in_transit(i))
        veh_in_transit_m(time_(i)+1,1) = quantity(i);
    end
end


%% Analysis
disp('3. Analysis...')
% check if the number of vehicles is constant over the simulation
total_vehicles = zeros(n_reb_periods, 1);
available_veh_per_interval = zeros(n_reb_periods, 1);
reb_veh_per_interval = zeros (n_reb_periods, 1);

for i =1 : n_reb_periods
    total_vehicles(i,1) = sum(available_veh_m(i,:)) + sum(veh_in_transit_m(i,:));
    available_veh_per_interval(i,1) = sum(available_veh_m(i,:));
    reb_veh_per_interval(i,1) = sum(reb_veh_m(i,:));
end

%% Reformat output
% to be in the format time, from, to, count
% this file will be the input file for offline rebalancing in
% amodController

reb_time = 0; % (is it beginning or end of the interval???)




%% Save to file
disp('4. Save rebalancing counts version 1...')
% the file will serve as an input for the simulation
filenameC = sprintf('rebalancingCounts_ecbd_per%d_st%d.txt', n_reb_periods, n_stations);
delimiter = ' ';
dlmwrite(filenameC, reb_veh_m,  delimiter);

disp('4. Save rebalancing counts version 1...')


disp('All done.')