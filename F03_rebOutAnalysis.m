close all; clear; clc;

simple_model = false;

% input files:
if (simple_model)
    % gurobi output file, file automatically generated b Gurobi
    gurobi_out = sprintf('/home/kasia/Documents/rebalancing_cpp/rebalancingMethods/rebalancing_offline/rebalancing_solution_simple.sol');
    % counts of origins at each station at each rebalancing interval
    % size: matrix n_rebalancing_intervals x nstations
    originFile = sprintf('sampleFiles/origCounts3x3.txt');
    % destination counts at each station at each rebalancing interval
    % size: matrix n_rebalancing_intervals x nstations
    destFile = sprintf('sampleFiles/destCounts3x3.txt');
    % facility file: facility_id, posX, posY; every line is a new facility
    facilityFile = sprintf('sampleFiles/stationsXY.txt');
    % vehicles in transit, matrix n_rebalancing_intervals x nstations
    intransitFile = sprintf('sampleFiles/inTransit3x3.txt');
    % estimated travel cost for a trip between stations
    travelcostFile = sprintf('sampleFiles/costM3x3.txt');
    % rebalancing interval in seconds
    reb_interval = 1; % 1 for sample files
    
else
    % simMobility model
    % gurobi output file, file automatically generated b Gurobi
    gurobi_out = sprintf('/home/kasia/Documents/rebalancing_cpp/rebalancingMethods/rebalancing_offline/rebalancing_solution.sol');
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
    % estimated travel cost for a trip between stations
    travelcostFile = sprintf('RebTimeInSecs34Stations.txt');
    % rebalancing interval in seconds
    reb_interval = 900; % in seconds
end


% output files:
% possible 3 different output files for rebalancing input:
% 1) time, from, to, count
% or each trip separately
% 2) time, from, to. From to can be as x and y or as a node_id

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
disp('2. Stations file...');
stationsData = dlmread(facilityFile, ' ', 0, 0);

f_ids = stationsData(:,1);
stationX = stationsData(:,2);
stationY = stationsData(:,3);

%% Import trips data - origin counts
disp('3. Origin counts, destination counts, in_transit counts and travel time from files...');
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

%% travel time between stations
ttFile_ = fopen(travelcostFile);
% find number of columns in the file
tline = fgetl(ttFile_);
fclose(ttFile_); % have to close and open the file again to make sure that the first line is not skipped later
N_stations_tt = length(find(tline==' ')) + 1;

% Keep track of the current state outside the loop:
nrows_tt = 0;     % Number of rows in current run
traveltimes = zeros(N_stations_tt, N_stations_tt);
% reopen the file
ttFile_ = fopen(travelcostFile);

while ~feof(ttFile_)
    current_line = fgetl(ttFile_);
    % Parse row of data from line
    row = textscan( current_line, '%u', 'delimiter', ' ' );
    row_double = cell2mat(row)';
    % Append the row to the current run data
    nrows_tt = nrows_tt + 1;
    traveltimes(nrows_tt,:) = row_double;
end

%% Sum available vehicles and rebalancing counts
disp('4. Sum available vehicles and rebalancing counts...')

% 1) vehicles available at each station for all time intervals
% (strcmp(GRB_var_name,'v_ti')) returns 1 if GRB_var_name is equal to the
% string, otherwise it returns zero
idle_veh_ix = (strcmp(GRB_var_name,'v_ti'));
available_veh_m = zeros(nrows, nstations);

% 2) rebalancing departing counts between stations for all time intervals
reb_veh = (strcmp(GRB_var_name,'r_tij'));
reb_dep_counts = zeros(nrows, nstations*nstations);
% rebalancing arriving counts between stations for all time intervals,
% which is rebalancing_dep + travel time
reb_arr_counts = zeros(nrows, nstations*nstations);
% to retrive and store rebalancing counts
vec_for_reb = 1:nstations*nstations;
reb_matrix = transpose(reshape(vec_for_reb, [nstations, nstations]));

for i =1 : length (GRB_var_name)
    if (idle_veh_ix(i))
        available_veh_m(GRB_time(i)+1, GRB_depSt_orSt(i)+1) = GRBSOL_quantity(i);
    end
    
    if (reb_veh(i))
        column_indx = reb_matrix(GRB_depSt_orSt(i)+1, GRB_arrSt_orSt(i)+1);
        reb_dep_counts(GRB_time(i) + 1, column_indx) = GRBSOL_quantity(i);
        if (GRBSOL_quantity(i) > 0)
            % shift the counts forward by travel time
            tt = traveltimes(GRB_depSt_orSt(i)+1, GRB_arrSt_orSt(i)+1);
            arr_time = rem (GRB_time(i) + 1 + tt, nrows*reb_interval);
            reb_arr_counts(arr_time, column_indx) = GRBSOL_quantity(i);
        end
    end
end

% Add travel time to the reb_departure time

%% Analysis
disp('5. Check if the number of vehicles is constant over the simulation...')
% check if the number of vehicles is constant over the simulation
total_vehicles = zeros(nrows, 1);
available_veh_per_interval = zeros(nrows, 1);
reb_veh_per_interval = zeros (nrows, 1);

for i = 1 : nrows
    total_vehicles(i,1) = sum(available_veh_m(i,:)) + sum(intransit_counts(i,:)) + sum(dest_counts(i,:)) + sum(reb_arr_counts(i,:));
    
    if (i > 1)
        if (total_vehicles(i) ~= total_vehicles(i - 1))
            X = sprintf('Number of vehicles NOT EQUAL at i = %d!==> %d != %d', i, total_vehicles(i), total_vehicles(i - 1));
            disp(X)
            %         else
            %             X = sprintf('EQUAL at i = %d!==> %d == %d', i, total_vehicles(i), total_vehicles(i - 1));
            %             disp(X)
        end
    end
end

% last interval comparison
if (total_vehicles(1) ~= total_vehicles(end))
    X = sprintf('Number of vehicles NOT EQUAL at i = %d!==> %d != %d', 1, total_vehicles(1), total_vehicles(end));
    disp(X)
    %         else
    %             X = sprintf('EQUAL at i = %d!==> %d == %d', i, total_vehicles(i), total_vehicles(i - 1));
    %             disp(X)
end

%% Reformat output
disp('6. Reformat output...')
% to be in the format time, from, to, count (time is the time when vehicle
% has to depart for destination)
% this file will be the inputed in amodController

rebalances_n2n = zeros(100000, 3); % preallocated space, line: time, node, node
rebalances_xy2xy = zeros(100000, 5); % preallocated space, line: time, originX, originY, destX, destY
rebalances_xy2xy_comb = zeros(100000, 6); % preallocated space, line: time, originX, originY, destX, destY, count
rebalances_n2n_comb = zeros(100000, 4); % preallocated space, line: time, node, node, count
index = 0;
index_c = 0;
% check for all rebalancing variables in the optimization output and write
% a new line for each rebalancing trip: time, origin, dest,
% version 2 : time, originX, originY, destX, destY
% version 3 : time, origin, dest, counts (there is no separate line for
% each trip)
for i = 1 : length (GRB_time)
    if (reb_veh(i))
        column_indx = reb_matrix(GRB_depSt_orSt(i)+1, GRB_arrSt_orSt(i)+1);
        % find origin station id
        origin_st = f_ids(GRB_depSt_orSt(i)+1);
        dest_st = f_ids(GRB_arrSt_orSt(i)+1);
        originX = stationX(GRB_depSt_orSt(i)+1);
        originY = stationY(GRB_depSt_orSt(i)+1);
        destX = stationX(GRB_arrSt_orSt(i)+1);
        destY = stationY(GRB_arrSt_orSt(i)+1);
        % time has to be converted to seconds based on rebalancing interval
        rebtime = (GRB_time(i)+1)*reb_interval;
        % all rebalancing trips are saved separately so for each count
        % prodice one trip
        for j = 1 : GRBSOL_quantity(i)
            index = index + 1;
            rebalances_n2n(index, :) = [rebtime; origin_st; dest_st];
            rebalances_xy2xy(index, :) = [rebtime; originX; originY; destX; destY];
        end
        index_c = index_c + 1;
        rebalances_xy2xy_comb(index_c, :) = [rebtime; originX; originY; destX; destY; GRBSOL_quantity(i)];
        rebalances_n2n_comb(index_c, :) = [rebtime; origin_st; dest_st; GRBSOL_quantity(i)];
    end
end

rebalances_n2n = rebalances_n2n(1:index, :);
rebalances_xy2xy = rebalances_xy2xy(1:index, :);
rebalances_xy2xy_comb = rebalances_xy2xy_comb(1:index_c, :);
rebalances_n2n_comb = rebalances_n2n_comb(1:index_c, :);

clearvars ans column_indx current_line dest_st destX destY facilityFile i j originX originY origin_st row tline X;

%% Save to file
disp('7. Save rebalancing counts version 1...')
% the file will serve as an input for the simulation
filename_out = sprintf('rebalancingCounts_sample_per%d_st%d.txt', nrows, nstations);
delimiter = ' ';
dlmwrite(filename_out, rebalances_n2n_comb,  delimiter); % or rebalances_xy2xy_comb or rebalances_xy2xy

disp('All done.')