close all; clear; clc;
%% Arrival Time Estimation
% for the offile rebalancing policy
% given the demand, travel cost matrix and station locations, the script
% aims to find the extimated arrival time for each customer plus the extra
% time needed for the vehicle to return to the station
% each line in the output file is in the format:
% origX, origY, destX, destY, bookingTime, arrivalTime

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
traveltime = zeros(length(origX),1); % in secs

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

%% Calculate travel time and arrival time 
% arrival time = request time + average waiting time + travel time
% travel time = distance travelled*ave_speed
% waiting time is normally distributed with mean = 6 mins (360secs) and
% standard deviation = 5 mins (300 secs)
% speed is normally distributed with mean = 22 km/h (6 m/s) and
% st.dev = 10 km/h (2.78 m/s) and minimum speed is 2m/s (7.2 km/h)

for i = 1:length(origX)
    
    waiting_time_distr = normrnd(360, 300);
    waiting_time = max(0, waiting_time_distr);
    speed_distr = normrnd(6, 2.78);
    speed = max(2, speed_distr);
    arrivalTime(i) = bookingTime(i) + O2Ddistance(i)/speed + waiting_time;
    traveltime(i) = O2Ddistance(i)/speed;
end

bookings_wArrivalT = [origX, origY, destX, destY, int32(bookingTime), int32(arrivalTime)];
bookings_wTT = [origX, origY, destX, destY, int32(bookingTime), int32(traveltime)];
%% Save to file

fileTOSave = sprintf('arrivalTimes.txt');
delimiter = ' ';
dlmwrite(fileTOSave, bookings_wArrivalT, delimiter);

TTfileTOSave = sprintf('trips_withTTfor%dstations.txt', length(f_ids));
delimiter = ' ';
dlmwrite(TTfileTOSave, bookings_wTT, delimiter);

%% find the nearest station for each origin and destiantion
closestStO = zeros(length(origX), 1);
closestStD = zeros(length(destX), 1);
for i = 1:length(origX)
    [o_distance, o_closestStationID ] = min(pdist2([origX(i,:), origY(i,:)], [stationX(:), stationY(:)]));
    [d_distance, d_closestStationID ] = min(pdist2([destX(i,:), destY(i,:)], [stationX(:), stationY(:)]));
    closestStO(i) = o_closestStationID;
    closestStD(i) = d_closestStationID;
end

new_bookings = [closestStO, closestStD, bookingTime, int32(arrivalTime)];

newBooking_f = sprintf('tripsBetween%dStation.txt', length(f_ids));
delimiter = ' ';
dlmwrite(newBooking_f, new_bookings, delimiter);
