close all; clear; clc;

%% Count number of people waiting at each station for every interval
% Given the trip data file consisting on trips origin, destination and time
% of the beginning of the trip find how many cars do we need to run
% autonomous mobility on demand system.
%
% input1: trips between station in the following format:
% request_time_sec, origin_station (node_id), destination_station (node_id), 
% traveltime, arrivalTime

% input2: stations in the format:
% node_id, pos_x, pos_y

% output:

%% Read files
% booking file
disp('1. Import trips between stations...')
bookingFile = sprintf('tripsBetween34Station2016-03-21.txt');
bookingData = dlmread(bookingFile, ' ', 0, 0);

booking_time = bookingData(:,1);
origin_station = bookingData(:,2);
dest_station = bookingData(:,3);
travel_time = bookingData(:,4);

% station file
disp('2. Import stations...')
facilityFile = sprintf('/Users/katarzyna/Dropbox/matlab/2016-03-Demand_generation/facility_location/stations_ecbd34.txt');
stationsData = dlmread(facilityFile, ' ', 0, 0);

f_ids = stationsData(:,1);
stationX = stationsData(:,2);
stationY = stationsData(:,3);

%% show graph for all the trips every deltaT (histogram)
disp('3. Plot histogram for all trips...')
figure()
hist(booking_time, 24*4); % histogram every 15 mins
figure()
hist(booking_time, 24); % histogram every hour

% maxN_cust15mins = max(hist(bookingTime, 24*4));
% maxN_cust1hr = max(hist(bookingTime, 24));

%% Find station ids and substitute them by the sequence number
% sequence number is based on the order in the stations' file
% this is to allow counting of departures and arrivals for each station
indices = zeros(length(f_ids),1);
for i = 1:length(f_ids)
   indices(i) = i; 
end

origin_id = zeros(length(booking_time),1);
dest_id = zeros(length(booking_time),1);
for i = 1: length(booking_time)
    origin_id(i) = find(f_ids == origin_station(i));
    dest_id(i) = find(f_ids == dest_station(i));
end

%% Departures and arrivals at each station
rebalancing_period_start = 0;
rebalancing_period_end = 15*60; % #mins in seconds
reb_delta = rebalancing_period_end - rebalancing_period_start;
dayLength = 60*60*24; % 24hrs in seconds
n_periods = ceil(dayLength/rebalancing_period_end); %number of reb periods

% Convert travel time to adjust for periods
travel_time_ = zeros(length(booking_time),1);
for i = 1: length(booking_time)
    travel_time_(i) = floor(travel_time(i)/reb_delta) + 1;
end

% number of departures and arrivals at each station
current_period = 1;
counter_orig = zeros(n_periods, length(f_ids));
counter_dest = zeros(n_periods, length(f_ids));
for i = 1: length(booking_time) % 
    % we are in the rebalancing interval
    if ((booking_time(i) >= rebalancing_period_start) && (booking_time(i) < rebalancing_period_end))   
        counter_orig (current_period, origin_id(i)) = counter_orig (current_period, origin_id(i)) + 1;
        if (rem(current_period + travel_time_(i), n_periods) ~= 0)
            counter_dest (rem(current_period + travel_time_(i), n_periods), dest_id(i)) = counter_dest (rem(current_period + travel_time_(i), n_periods), dest_id(i)) + 1; 
        else
            counter_dest (n_periods, dest_id(i)) = counter_dest (n_periods, dest_id(i)) + 1;
        end
    end
    
    if (booking_time(i) > rebalancing_period_end)
        rebalancing_period_start = rebalancing_period_end;
        rebalancing_period_end = rebalancing_period_end + reb_delta;
        current_period = floor(booking_time(i)/reb_delta) + 1;    
    end
end

%% Save to file

fileTOSave_orig = sprintf('origCounts_rebEvery%d_stations%d.txt', reb_delta, length(stationX));
delimiter = ' ';
dlmwrite(fileTOSave_orig, counter_orig,  delimiter);

fileTOSave_dest = sprintf('destCounts_rebEvery%d_stations%d.txt', reb_delta, length(stationX));
delimiter = ' ';
dlmwrite(fileTOSave_dest, counter_dest, delimiter);
