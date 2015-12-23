close all; clear; clc;

%% Count number of people waiting at each station for every interval
% Given the trip data file consisting on trips origin, destination and time
% of the beginning of the trip find how many cars do we need to run
% autonomous mobility on demand system.
%
% To run this script, first you should run the arrivalTimeEstimation.m to
% get the estimated arrival time of each trip
%% Read files
% booking file
bookingFile = sprintf('trips_withTTfor21stations.txt');
bookingData = dlmread(bookingFile, ' ', 0, 0);

origX = bookingData(:,1);
origY = bookingData(:,2);
bookingTime = bookingData(:,5);

% station file
facilityFile = sprintf('inputDemand/ecbd_stations21.txt');

stationsData = dlmread(facilityFile, ' ', 0, 0);

f_ids = stationsData(:,1);
stationX = stationsData(:,2);
stationY = stationsData(:,3);

%% show the graph for all the trips every deltaT (histogram)

figure()
hist(bookingTime, 24*4);
figure()
hist(bookingTime, 24);

maxN_cust15mins = max(hist(bookingTime, 24*4));
maxN_cust1hr = max(hist(bookingTime, 24));

%% Outgoing trips
rebalancing_period = 15*60; % #mins in seconds
dayLength = 60*60*24; % 24hrs in seconds
n_periods = ceil(dayLength/rebalancing_period); %number of reb periods
% trips count at each station and every time interval
trips_count = zeros(n_periods, length(stationX));

% find the nearest station for every origin
Dist_matrix_orig = pdist2([origX(:), origY(:)], [stationX(:), stationY(:)]);
minDist_orig = zeros(length(origX),1);
for ii = 1:length(origX)
    [minDist_, ind] = min(Dist_matrix_orig(ii,:));
    minDist_orig(ii) = ind;
end

% for each time interval find how many trips is within the interval
trips_assigned = [];
counter_origin = zeros(n_periods, length(f_ids));
for i = 1: n_periods
    for j = 1 : length(bookingTime)
        if(bookingTime(j) <= i*rebalancing_period)
            % check which station the trip belongs to
            station_ = minDist_orig(j);
            % add it to the station counter
            counter_origin(i, station_) = counter_origin(i, station_) + 1;
            trips_assigned = [trips_assigned j];
        end % end if
    end % end for each customer

    % delete trips from the array
    bookingTime(trips_assigned) = [];
    origX(trips_assigned) = [];
    origY(trips_assigned) = [];
    minDist_orig(trips_assigned) = [];
    trips_assigned = [];
end

%% Incoming trips
arrivalTime = bookingData(:,6);
destX = bookingData(:,3);
destY = bookingData(:,4);

incoming_trips = [arrivalTime, destX, destY];
incoming_trips_sorted = sortrows(incoming_trips, 1);

arrivalTime = incoming_trips_sorted(:,1);
destX = incoming_trips_sorted(:,2);
destY = incoming_trips_sorted(:,3);

% find the nearest station for every destination
Dist_matrix_dest = pdist2([destX(:), destY(:)], [stationX(:), stationY(:)]);
minDist_dest = zeros(length(destX),1);
for ii = 1:length(destX)
    [minDist_, ind] = min(Dist_matrix_dest(ii,:));
    minDist_dest(ii) = ind;
end

% for each time interval find how many trips is within the interval
counter_dest = zeros(n_periods, length(f_ids));
for i = 1: n_periods
    for j = 1 : length(arrivalTime)
        if(arrivalTime(j) <= i*rebalancing_period)
            % check which station the trip belongs to
            station_ = minDist_dest(j);
            % add it to the station counter
            counter_dest(i, station_) = counter_dest(i, station_) + 1;
            trips_assigned = [trips_assigned j];
        end % end if
    end % end for each customer

    % delete trips from the array
    arrivalTime(trips_assigned) = [];
    destX(trips_assigned) = [];
    destY(trips_assigned) = [];
    minDist_dest(trips_assigned) = [];
    trips_assigned = [];
end

%% Save to file

fileTOSave_orig = sprintf('origCounts_reb%d_stations%d_updated.txt', rebalancing_period, length(stationX));
delimiter = ' ';
dlmwrite(fileTOSave_orig, counter_origin, delimiter);

fileTOSave_dest = sprintf('destCounts_reb%d_stations%d_updated.txt', rebalancing_period, length(stationX));
delimiter = ' ';
dlmwrite(fileTOSave_dest, counter_dest, delimiter);

%% generate a matrix of trips for each interval
% change coordinates for the station number

% find the nearest station for each origin and destiantion
origX = bookingData(:,1);
origY = bookingData(:,2);
destX = bookingData(:,3);
destY = bookingData(:,4);
bookingTime = bookingData(:,5);
arrivalTime = bookingData(:,6);

closestStO = zeros(length(origX), 1);
closestStD = zeros(length(destX), 1);
for i = 1:length(origX)
    [o_distance, o_closestStationID ] = min(pdist2([origX(i,:), origY(i,:)], [stationX(:), stationY(:)]));
    [d_distance, d_closestStationID ] = min(pdist2([destX(i,:), destY(i,:)], [stationX(:), stationY(:)]));
    closestStO(i) = o_closestStationID;
    closestStD(i) = d_closestStationID;
end

new_bookings = [closestStO, closestStD, bookingTime, int32(arrivalTime)];

newBooking_f = sprintf('tripsAggregatedBy%dStation.txt', length(f_ids));
delimiter = ' ';
dlmwrite(newBooking_f, new_bookings, delimiter);