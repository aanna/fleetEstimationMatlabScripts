close all; clear; clc;

%% Count number of people waiting at each station for every interval
% Given the trip data file consisting on trips origin, destination, booking
% time and arrival time, under the assumption that everyone is served
% within the interval when booking arrived
%
% input: trips between station in the following format:
% request_time_sec, origin_station (node_id), destination_station (node_id),
% traveltime, arrivalTime

% output:
% file consisting of the number of vehicles travelling now to each station

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
    travel_time_(i) = floor(travel_time(i)/reb_delta) + 1;
end

% number of departures and arrivals at each station
current_period = 1;
arriving_dest = zeros(n_periods, length(f_ids));
in_transit = zeros(n_periods, length(f_ids));

for i = 1: length(booking_time) %
    % check if we are in the rebalancing interval
    % this causes problem because if this is not true then we skip that
    % line which should not happen
    if ((booking_time(i) >= rebalancing_period_start) && (booking_time(i) < rebalancing_period_end))
        
        if (travel_time_(i) > 1)
            for j = 1 : travel_time_(i)
                in_transit (current_period + j - 1, dest_id(i)) = in_transit (current_period + j - 1, dest_id(i)) + 1;
            end
        else
            % do nothing
        end
        
    elseif (booking_time(i) >= rebalancing_period_end)
        
        if (travel_time_(i) > 1)
            for j = 1 : travel_time_(i)
                in_transit (current_period + j, dest_id(i)) = in_transit (current_period + j, dest_id(i)) + 1;
            end
        else
            % do nothing
        end
        
        % update the indices
        rebalancing_period_start = rebalancing_period_end;
        rebalancing_period_end = rebalancing_period_end + reb_delta;
        current_period = current_period + 1;
    end
    
end

still_in_transit = in_transit - arriving_dest;

%% Save to file

fileTOSave_dest = sprintf('inTransit%d_stations%d.txt', reb_delta, length(f_ids));
delimiter = ' ';
dlmwrite(fileTOSave_dest, arriving_dest, delimiter);
