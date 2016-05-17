close all; clear; clc;

% input files:
% gurobi output file, file automatically generated b Gurobi
gurobi_out = sprintf('/home/kasia/Documents/rebalancing_cpp/rebalancingMethods/rebalancing_offline/rebalancing_solution_simple.sol');
% counts of origins at each station at each rebalancing interval
% size: matrix n_rebalancing_intervals x nstations
originFile = sprintf('origCounts_rebEvery900_stations34.txt');
% destination counts at each station at each rebalancing interval
% size: matrix n_rebalancing_intervals x nstations
destFile = sprintf('destCounts_rebEvery900_stations34.txt');
% facility file: facility_id, posX, posY; every line is a new facility
facilityFile = sprintf('stations_ecbd34.txt');
% vehicles in transit, matrix n_rebalancing_intervals x nstations
intransitFile = sprintf('inTransit900_stations34.txt');

% output files:


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

clearvars gurobi_out delimiter startRow formatSpec fileID dataArray ans;

%% Import stations
% disp('2. Stations file...');
% stationsData = dlmread(facilityFile, ' ', 0, 0);
% 
% f_ids = stationsData(:,1);
% stationX = stationsData(:,2);
% stationY = stationsData(:,3);

%% Import trips data - origin counts
disp('3. Origin destination and in_transit counts from files...');
% we do not know ahead how many columns (stations) and rows (intervals) is
% in the file

% origin counts
originFile_ = fopen(originFile);

% find number of columns in the file
tline = fgetl(originFile_);
fclose(originFile_); % have to close and open the file again to make sure that the first line is not skipped later
N_stations_orig = length(find(tline==' ')) + 1;

% Keep track of the current state outside the loop:
nrows_orig = 0;     % Number of rows in current run
origin_counts = zeros(1000, N_stations_orig);  % Prealocated for speed improvement, later it will be trunckated to the correct number of rows
% reopen the file
originFile_ = fopen(originFile);

while ~feof(originFile_)
    current_line = fgetl(originFile_);
    % Parse row of data from line
    row = textscan( current_line, '%u', 'delimiter', ' ' );
    row_double = cell2mat(row)';
    % Append the row to the current run data
    nrows_orig = nrows_orig + 1; 
    origin_counts(nrows_orig,:) = row_double;
end

% Remove all rows from nrows until the end of matrix
origin_counts = origin_counts(1:nrows_orig,:);

%% Import trips data - destination counts
% dest_counts; the same procedure as for origin_counts
destFile_ = fopen(destFile);
% find number of columns in the file
tline = fgetl(destFile_);
fclose(destFile_); % have to close and open the file again to make sure that the first line is not skipped later
N_stations_dest = length(find(tline==' ')) + 1;

% Keep track of the current state outside the loop:
nrows_dest = 0;     % Number of rows in current run
dest_counts = zeros(1000, N_stations_dest);  % Prealocated for speed improvement, later it will be trunckated to the correct number of rows
% reopen the file
destFile_ = fopen(destFile);

while ~feof(destFile_)
    current_line = fgetl(destFile_);
    % Parse row of data from line
    row = textscan( current_line, '%u', 'delimiter', ' ' );
    row_double = cell2mat(row)';
    % Append the row to the current run data
    nrows_dest = nrows_dest + 1; 
    dest_counts(nrows_dest,:) = row_double;
end

% Remove all rows from nrows until the end of matrix
dest_counts = dest_counts(1:nrows_dest,:);

% check if N_stations_orig == N_stations_dest
if (N_stations_orig == N_stations_orig && nrows_orig == nrows_dest)
    disp('CORRECT: N_stations_orig == N_stations_dest && nrows_orig == nrows_dest')
    % create a general variable N_stations
    nstations =  N_stations_orig;
    nrows = nrows_orig;
    clearvars N_stations_orig N_stations_dest nrows_orig nrows_dest destFile_ destFile originFile_ originFile ;
else
    disp('WRONG: N_stations_orig != N_stations_dest or && nrows_orig != nrows_dest')
end

%% in transit counts
% dest_counts; the same procedure as for origin_counts
intransitFile_ = fopen(intransitFile);
% find number of columns in the file
tline = fgetl(intransitFile_);
fclose(intransitFile_); % have to close and open the file again to make sure that the first line is not skipped later
N_stations_intr = length(find(tline==' ')) + 1;

% Keep track of the current state outside the loop:
nrows_intr = 0;     % Number of rows in current run
intransit_counts = zeros(1000, N_stations_intr);  % Prealocated for speed improvement, later it will be trunckated to the correct number of rows
% reopen the file
intransitFile_ = fopen(intransitFile);

while ~feof(intransitFile_)
    current_line = fgetl(intransitFile_);
    % Parse row of data from line
    row = textscan( current_line, '%u', 'delimiter', ' ' );
    row_double = cell2mat(row)';
    % Append the row to the current run data
    nrows_intr = nrows_intr + 1; 
    intransit_counts(nrows_intr,:) = row_double;
end

% Remove all rows from nrows until the end of matrix
intransit_counts = intransit_counts(1:nrows_intr,:);

% check if N_stations_orig == N_stations_dest
if (N_stations_intr == nstations && nrows == nrows_intr)
    disp('CORRECT: N_stations_intr == nstations && nrows == nrows_intr')
    % create a general variable N_stations
    clearvars N_stations_intr nrows_intr intransitFile_ intransitFile;
else
    disp('WRONG: N_stations_intr != nstations && nrows == nrows_intr')
end

%% Read available vehicles and rebalancing counts
disp('4. Read available vehicles and rebalancing counts...')
% 1) vehicles available at each station for all time intervals
% (strcmp(GRB_var_name,'v_ti')) returns 1 if GRB_var_name is equal to the
% string, otherwise it returns zero
idle_veh_ix = (strcmp(GRB_var_name,'v_ti'));
available_veh_m = zeros(nrows, nstations);
% 2) rebalancing counts between stations for all time intervals
reb_veh = (strcmp(GRB_var_name,'r_tij'));
reb_veh_m = zeros(nrows, nstations*nstations);

% reb_matrix = transpose(reshape(vec_for_reb, [nstations, nstations]));

for i =1 : length (GRB_var_name)
    if (idle_veh_ix(i))
        available_veh_m(GRB_time(i)+1, GRB_arrSt_orSt(i)+1) = GRBSOL_quantity(i);
    end
    
%     if (reb_veh(i))
%         column_indx = reb_matrix(GRB_depSt_orSt(i)+1, GRB_arrSt_orSt(i)+1);
%         reb_veh_m(GRB_time(i)+1, column_indx) = GRBSOL_quantity(i);
%     end
end


%% Analysis
disp('3. Analysis...')
% check if the number of vehicles is constant over the simulation
total_vehicles = zeros(nrows, 1);
available_veh_per_interval = zeros(nrows, 1);
reb_veh_per_interval = zeros (nrows, 1);

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