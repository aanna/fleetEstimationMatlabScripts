close all; clear; clc;
%% Arrival Time Estimation
% for the offile rebalancing policy

% given the demand, travel cost matrix and station locations, the script
% aims to find the estimated arrival time for each customer plus the extra
% time needed for the vehicle to return to the station

% input: booking file in format: 
% booking_id, booking_time, customer_id, origin_x, origin_y, destination_x, destination_y 
% intput: station file in format: station_node_id, station_x, station_y
% input: nodes in the zone: node_id, pos_x, pos_y

% output1 is in the format:
% origX, origY, destX, destY, bookingTime, travel_time, arrival_time
% output2 is in the format:
% bookingTime, station_origin, station_dest, traveltime, arrivalTime

%% Read files
% booking file
disp('1. Import bookings...')
bookingFile = sprintf('/Users/katarzyna/Dropbox/matlab/2016-03-Demand_generation/boookings_ecbd_808356.txt');
bookingData = dlmread(bookingFile, ' ', 0, 0);

bookingTime = bookingData(:,2);
origX = bookingData(:,4);
origY = bookingData(:,5);
destX = bookingData(:,6);
destY = bookingData(:,7);

% station file
disp('2. Import stations...')
facilityFile = sprintf('/Users/katarzyna/Dropbox/matlab/2016-03-Demand_generation/facility_location/stations_ecbd34.txt');
stationsData = dlmread(facilityFile, ' ', 0, 0);

f_ids = stationsData(:,1);
stationX = stationsData(:,2);
stationY = stationsData(:,3);

% nodes file
disp('3. Import nodes...')
nodesFile = sprintf('/Users/katarzyna/Dropbox/matlab/2016-03-Demand_generation/facility_location/input/ecbd_nodes.csv');
nodesData = dlmread(nodesFile, ',', 0, 0);

node_ids = nodesData(:,1);
nodeX = nodesData(:,2);
nodeY = nodesData(:,3);
%% for each booking, find distance from origin to destination
disp('4. Find distance from origins to destinations...')
O2Ddistance = zeros(length(origX),1); % in m
arrivalTime = zeros(length(origX),1); % in secs
traveltime = zeros(length(origX),1); % in secs

euclideanDist = true;
if (euclideanDist)
    % do the euclidian distances between origin and destination
    % distances are increased by 20% due to euclidean distance
    correction_factor = 1.2; % should be a random value from an interval
    for i=1:length(origX)
        O2Ddistance(i) = correction_factor*pdist2([origX(i), origY(i)], [destX(i), destY(i)]);
    end
else
    %find the nearest node for each origin and destination and find the
    %distance between nodes from the shortest_distances output file
    node_origin = pdist2([origX(:), origY(:)], [nodeX(:), nodeY(:)]);
    node_dest = pdist2([destX(:), destY(:)], [nodeX(:), nodeY(:)]);
    
    costestNode_orig = zeros(length(origX),1);
    for ii = 1:length(origX)
        [minDist_, ind] = min(node_origin(ii,:));
        costestNode_orig(ii) = ind;
    end
    
    closestNode_dest = zeros(length(destX),1);
    for ii = 1:length(destX)
        [minDist_, ind] = min(node_dest(ii,:));
        closestNode_dest(ii) = ind;
    end
end

%% Calculate travel time and arrival time
% arrival time = request time + average waiting time + travel time
% travel time = distance travelled*ave_speed
% waiting time is normally distributed with mean = 6 mins (360secs) and
% standard deviation = 5 mins (300 secs)
% speed is normally distributed with mean = 22 km/h (6 m/s) and
% st.dev = 10 km/h (2.78 m/s) and minimum speed is 2m/s (7.2 km/h)

disp('5. Calculate travel time...')
for i = 1:length(origX)
    
    waiting_time_distr = normrnd(360, 300);
    waiting_time = max(0, waiting_time_distr);
    speed_distr = normrnd(6, 2.78);
    speed = max(2, speed_distr);
    arrivalTime(i) = bookingTime(i) + O2Ddistance(i)/speed + waiting_time;
    traveltime(i) = O2Ddistance(i)/speed;
end

bookings_wATT = [int32(bookingTime), origX, origY, destX, destY, int32(traveltime), int32(arrivalTime)];
%% Save to file
disp('6. SAve bookings with travel time...')
delimiter = ' ';
trips_ArrivalAndTT = sprintf('trips_ArrivalAndTT%dstations.txt', length(f_ids));
dlmwrite(trips_ArrivalAndTT, bookings_wATT, delimiter);

%% trips between the stations only
% find the nearest station for each origin and destiantion
disp('7. Find the nearest station for each origin and destiantion...')
closestStO = [];
closestStD = [];
ind_delete = [];
for i = 1:length(origX)
    [o_distance, o_closestStationID ] = min(pdist2([origX(i,:), origY(i,:)], [stationX(:), stationY(:)]));
    [d_distance, d_closestStationID ] = min(pdist2([destX(i,:), destY(i,:)], [stationX(:), stationY(:)]));
    if (o_closestStationID ==  d_closestStationID)
        % ind to be removed from booking
        ind_delete = [ind_delete; i];
    else
        closestStO = [closestStO; o_closestStationID];
        closestStD = [closestStD; d_closestStationID];
    end
end



bookingTime(ind_delete) = [];
traveltime(ind_delete) = [];
arrivalTime(ind_delete) = [];

disp('8. Save bookings between stations...')
new_bookings = [ bookingTime, closestStO, closestStD, int32(traveltime), int32(arrivalTime)];
newBooking_f = sprintf('tripsBetween%dStationV2.txt', length(f_ids));
delimiter = ' ';
dlmwrite(newBooking_f, new_bookings, delimiter);

disp('All done.')