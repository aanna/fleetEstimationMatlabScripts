close all; clear; clc;

% input files:
% new output file in format time, from_node, to_node, count

simple_model = false;

if (simple_model)
    n_stations = 3;
    n_reb_periods = 3;
    rebalancing_interval = 1;
    vec_for_reb = 1:n_stations*n_stations;
    reb_matrix = transpose(reshape(vec_for_reb, [n_stations, n_stations]));
else
    n_stations = 34;
    n_reb_periods = 96;
    rebalancing_interval = 24*60*60/n_reb_periods; % in seconds
    vec_for_reb = 1:n_stations*n_stations;
    reb_matrix = transpose(reshape(vec_for_reb, [n_stations, n_stations]));
end

% input files:
gurobi_out = sprintf('/home/kasia/Documents/rebalancing_cpp/rebalancingMethods/rebalancing_offline/rebalancing_solution_simple.sol');
originFile = sprintf('origCounts_rebEvery900_stations34.txt');
destFile = sprintf('destCounts_rebEvery900_stations34.txt');
facilityFile = sprintf('stations_ecbd34.txt');

%% Import output from Gurobi cpp
disp('1. Import output from Gurobi cpp...')
delimiter = {',',' '};
startRow = 3;
formatSpec = '%s%f%f%f%f%[^\n\r]';
fileID = fopen(gurobi_out,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

GRB_var_name = dataArray{:, 1};
GRB_time = dataArray{:, 2};
GRB_depSt_orSt = dataArray{:, 3}; % depending on thr var_name, i.e., if reb_veh
%then this is departing station, if available_veh then @ which station
GRB_arrSt_orSt = dataArray{:, 4}; % same as above
GRBSOL_quantity = dataArray{:, 5};

clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% Import stations
disp('2. Stations file...');
stationsData = dlmread(facilityFile, ' ', 0, 0);

f_ids = stationsData(:,1);
stationX = stationsData(:,2);
stationY = stationsData(:,3);

%% Import trips data - origin counts
disp('3. Origin and destination counts from files...');
% we do not know ahead how many columns (stations) and rows (intervals) is
% in the file

% origin counts
originFile_ = fopen(originFile);

% find number of columns in the file
tline = fgetl(originFile_);
fclose(originFile_); % have to close and open the file again to make sure that the first line is not skipped later
N = length(find(tline==' ')) + 1;

% Keep track of the current state outside the loop:
nrows = 0;     % Number of rows in current run
origin_counts = zeros(1000, N);  % Prealocated for speed improvement, later it will be trunckated to the correct number of rows
% reopen the file
originFile_ = fopen(originFile);

while ~feof(originFile_)
    current_line = fgetl(originFile_);
    % Parse row of data from line
    row = textscan( current_line, '%u', 'delimiter', ' ' );
    row_double = cell2mat(row)';
    % Append the row to the current run data
    nrows = nrows + 1; 
    origin_counts(nrows,:) = row_double;
end

% Remove all rows from nrows until the end of matrix
origin_counts = origin_counts(1:nrows,:);

%% Import trips data - destination counts
% dest_counts; the same procedure as for origin_counts
destFile_ = fopen(destFile);
% find number of columns in the file
tline = fgetl(destFile_);
fclose(destFile_); % have to close and open the file again to make sure that the first line is not skipped later
N = length(find(tline==' ')) + 1;

% Keep track of the current state outside the loop:
nrows = 0;     % Number of rows in current run
dest_counts = zeros(1000, N);  % Prealocated for speed improvement, later it will be trunckated to the correct number of rows
% reopen the file
destFile_ = fopen(destFile);

while ~feof(destFile_)
    current_line = fgetl(destFile_);
    % Parse row of data from line
    row = textscan( current_line, '%u', 'delimiter', ' ' );
    row_double = cell2mat(row)';
    % Append the row to the current run data
    nrows = nrows + 1; 
    dest_counts(nrows,:) = row_double;
end

% Remove all rows from nrows until the end of matrix
dest_counts = dest_counts(1:nrows,:);

%% Generate 3 matrices
disp('4. Generate 3 matrices...')
% 1) vehicles available at each station for all time intervals
available_veh = (strcmp(GRB_var_name,'v_ti'));
available_veh_m = zeros(n_reb_periods, n_stations);
% 2) rebalancing counts between stations for all time intervals
reb_veh = (strcmp(GRB_var_name,'r_tij'));
reb_veh_m = zeros(n_reb_periods, n_stations*n_stations);

for i =1 : length (GRB_var_name)
    if (available_veh(i))
        available_veh_m(GRB_time(i)+1, GRB_arrSt_orSt(i)+1) = GRBSOL_quantity(i);
    end
    
    if (reb_veh(i))
        column_indx = reb_matrix(GRB_depSt_orSt(i)+1, GRB_arrSt_orSt(i)+1);
        reb_veh_m(GRB_time(i)+1, column_indx) = GRBSOL_quantity(i);
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



disp('All done.')