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
% file consisting of the number of trips originating and arriving at all
% stations

%% Read files
% booking file
disp('1. Import trips between stations...')
bookingFile = sprintf('tripsBetween34Station2016-05-09.txt');
bookingData = dlmread(bookingFile, ' ', 0, 0);

booking_time = bookingData(:,1);
origin_station = bookingData(:,2);
dest_station = bookingData(:,3);
travel_time = bookingData(:,4);

% station file
disp('2. Import stations...')
facilityFile = sprintf('stations_ecbd34.txt');
stationsData = dlmread(facilityFile, ' ', 0, 0);

f_ids = stationsData(:,1);

%% show graph for all the trips every deltaT (histogram)
disp('3. Plot histogram for all trips...')
figure()
hist(booking_time, 24*4); % histogram every 15 mins
% figure()
% hist(booking_time, 24); % histogram every hour

% maxN_cust15mins = max(hist(bookingTime, 24*4));
% maxN_cust1hr = max(hist(bookingTime, 24));

%% Find station ids and substitute them by the sequence number
disp('4. Find station ids...')
% sequence number is based on the order in the stations' file
% this is to allow counting of departures and arrivals for each station

origin_id = zeros(length(booking_time),1);
dest_id = zeros(length(booking_time),1);
for i = 1: length(booking_time)
    origin_id(i) = find(f_ids == origin_station(i));
    dest_id(i) = find(f_ids == dest_station(i));
end

%% Departures and arrivals at each station
disp('5. Number of departures and arrivals at each station...')
rebalancing_period_start = 0;
rebalancing_period_end = 15*60; % #mins in seconds
reb_delta = rebalancing_period_end - rebalancing_period_start;
dayLength = 60*60*24; % 24hrs in seconds
n_periods = ceil(dayLength/rebalancing_period_end); %number of reb periods

% Convert travel time to adjust for periods
travel_time_ = zeros(length(booking_time),1);
for i = 1: length(booking_time)
    travel_time_(i) = ceil(travel_time(i)/reb_delta);
end

% number of departures and arrivals at each station
current_period = 1;
counter_orig = zeros(n_periods, length(f_ids));
counter_dest = zeros(n_periods, length(f_ids));

for i = 1: length(booking_time) %
    % check if we are in the rebalancing interval
    if ((booking_time(i) >= rebalancing_period_start) && (booking_time(i) < rebalancing_period_end))
        counter_orig (current_period, origin_id(i)) = counter_orig (current_period, origin_id(i)) + 1;
        
        % rem(current_period + travel_time_(i), n_periods) returns the time
        % period when the vehicle will arrive at destination. It is
        % reminder because, i.e., when the trip starts at the end of the
        % day then it may arrive at the beginning of the next day. If the
        % rem is zero then we are at the last period.
        
        if (rem(current_period + travel_time_(i), n_periods) ~= 0)
            counter_dest (rem(current_period + travel_time_(i), n_periods), dest_id(i)) = counter_dest (rem(current_period + travel_time_(i), n_periods), dest_id(i)) + 1;
        else
            % we are arriving in the last interval
            counter_dest (n_periods, dest_id(i)) = counter_dest (n_periods, dest_id(i)) + 1;
        end
        
        % elseif handles the line which does not satisfy the condition
        % (otherwise it would be neglected)
    elseif (booking_time(i) >= rebalancing_period_end)   
        % this is to prevent skipping the first line which does not
        % satisfy the above condition
        counter_orig (current_period + 1, origin_id(i)) = counter_orig (current_period + 1, origin_id(i)) + 1;
        
        if (rem(current_period + 1 + travel_time_(i), n_periods) ~= 0)
            counter_dest (rem(current_period + 1 + travel_time_(i), n_periods), dest_id(i)) = counter_dest (rem(current_period + 1 + travel_time_(i), n_periods), dest_id(i)) + 1;
        else
            counter_dest (n_periods, dest_id(i)) = counter_dest (n_periods, dest_id(i)) + 1;
        end
        
        % update the indices
        rebalancing_period_start = rebalancing_period_end;
        rebalancing_period_end = rebalancing_period_end + reb_delta;
        current_period = current_period + 1;
    end
    
end

%% check if the sum of origins and destinations is correct
sum_origins = sum(counter_orig(:));
sum_dest = sum(counter_dest(:));

if (sum_origins == sum_dest && sum_origins == length(booking_time))
    disp('CORRECT: sum_origins == sum_dest && sum_origins == length(booking_time)')
else
    disp('WRONG: sum_origins ~= sum_dest || sum_origins ~= length(booking_time)')
end

%% Save to file
disp('6. Save to file...')
fileTOSave_orig = sprintf('origCounts_rebEvery%d_stations%d.txt', reb_delta, length(f_ids));
delimiter = ' ';
dlmwrite(fileTOSave_orig, counter_orig,  delimiter);

fileTOSave_dest = sprintf('destCounts_rebEvery%d_stations%d.txt', reb_delta, length(f_ids));
delimiter = ' ';
dlmwrite(fileTOSave_dest, counter_dest, delimiter);

%% number of vehicles in transit at each period of time
disp('7. Number of vehicles in transit at each period of time...')
current_period = 1;
rebalancing_period_start = 0;
rebalancing_period_end = 15*60; % #mins in seconds
in_transit = zeros(n_periods, length(f_ids));

for i = 1: length(booking_time) %
    
    if (travel_time_(i) > 1) % at least 2 time intervals so the vehicle is missing in the simulation for at least one period
        for j = 1 : (travel_time_(i) - 1) % minus one because we do not count the period when vehicle arrives in destination

        end
    end
    
    
end

%% sum of arrivals + departures + in_transit at each time step
disp('8. Sum of arrivals + departures + in_transit at each time step...')
total_vehicles = zeros(n_periods, 1);
for i = 1: length(total_vehicles)
    
    total_vehicles(i) = sum(in_transit(i,:)) + sum(counter_orig(i,:)) + sum(counter_dest(i,:));
    
end

% check
disp('9. Checking number of vehicles at every interval...')

for i = 1: length(total_vehicles)
    
    if (i < length(total_vehicles))

        if (total_vehicles(i) ~= total_vehicles(i + 1))
            X = sprintf('NOT EQUAL!: %d != %d', total_vehicles(i), total_vehicles(i + 1));
            disp(X)
        else
            X = sprintf('EQUAL!: %d == %d', total_vehicles(i), total_vehicles(i + 1));
            disp(X)
        end
    else

        if (total_vehicles(i) ~= total_vehicles(1))
            X = sprintf('NOT EQUAL!: %d != %d', total_vehicles(i), total_vehicles(1));
            disp(X)
        else
            X = sprintf('EQUAL!: %d == %d', total_vehicles(i), total_vehicles(1));
            disp(X)
        end
    end
end

%% Save to file
disp('10. Saving to file...')
fileTOSave_transit = sprintf('inTransitReb%dStations%d%.txt', reb_delta, length(f_ids));
delimiter = ' ';
dlmwrite(fileTOSave_transit, in_transit, delimiter);

disp('All done.')
