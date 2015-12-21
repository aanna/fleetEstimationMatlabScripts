close all; clear; clc;
%% Arrival Time Estimation
% for the offile rebalancing policy
% given the demand, travel cost matrix and station locations, the script
% aims to find the extimated arrival time for each customer plus the extra
% time needed for the vehicle to return to the station

%% Read files
% booking file
bookingFile = sprintf('inputDemand/ecbd_sorted_bookings.txt');
bookingData = dlmread(bookingFile, ' ', 0, 0);

bookingTime = bookingData(:,2);
origX = bookingData(:,4);
origY = bookingData(:,5);
destX = bookingData(:,6);
destY = bookingData(:,7);

% station file
facilityFile = sprintf('inputDemand/ecbd_stations21.txt');
stationsData = dlmread(facilityFile, ' ', 0, 0);

f_ids = stationsData(:,1);
stationX = stationsData(:,2);
stationY = stationsData(:,3);

% nodes file
nodesFile = sprintf('inputDemand/ecbd_stations21.txt');
nodesData = dlmread(nodesFile, ' ', 0, 0);

node_ids = nodesData(:,1);
nodeX = nodesData(:,2);
nodeY = nodesData(:,3);
%% for each booking, find distance from origin to destination
O2Ddistance = zeros(length(origX),1); % in m
arrivalTime = zeros(length(origX),1); % in secs

euclideanDist = true;
if (euclideanDist)
    % do the euclidian distances between origin and destination
    % distances are increased by 20% due to euclidean distance
    correction_factor = 1.2;
    for i=1:length(origX)
        O2Ddistance(i) = correction_factor*pdist2([origX(i), origY(i)], [destX(i), destY(i)]);
    end
else
    %find the nearest node for each origin and destination and find the
    %distance between nodes from the shortest_distances output file
    Dist_matrix_origin = pdist2([origX(:), origY(:)], [nodeX(:), nodeY(:)]);
    Dist_matrix_dest = pdist2([destX(:), destY(:)], [nodeX(:), nodeY(:)]);
    
    minDist_orig = zeros(length(origX),1);
    for ii = 1:length(origX)
        [minDist_, ind] = min(Dist_matrix_orig(ii,:));
        minDist_orig(ii) = ind;
    end
    
    minDist_dest = zeros(length(destX),1);
    for ii = 1:length(destX)
        [minDist_, ind] = min(Dist_matrix_dest(ii,:));
        minDist_dest(ii) = ind;
    end
end

%% Calculate the arrival time = request time + average waiting time +
% distance travelled*ave_speed
% waiting time is normally distributed with mean = 6 mins (360secs) and
% standard deviation = 5 mins (300 secs)
% speed is normally distributed with mean = 22 km/h (6 m/s) and
% st.dev = 30 km/h (8.33 m/s) and minimum speed is 2m/s (7.2 km/h)

for i = 1:length(origX)
    
    waiting_time_distr = normrnd(360, 300);
    waiting_time = max(0, waiting_time_distr);
    speed_distr = normrnd(6, 8.33);
    speed = max(2, speed_distr);
    arrivalTime(i) = bookingTime(i) + O2Ddistance(i)/speed + waiting_time;
    
end

new_bookings = [origX, origY, destX, destY, bookingTime, arrivalTime];
%% Save to file

fileTOSave = sprintf('arrivalTimes.txt');
delimiter = ' ';
dlmwrite(fileTOSave, new_bookings, delimiter);
