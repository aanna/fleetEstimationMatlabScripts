close all; clear; clc;
%% Arrival Time Estimation
% for the offile rebalancing policy

% given the demand, travel cost matrix and station locations, the script
% aims to find the estimated arrival time for each customer

% the file estimates the number of trips between stations and travel
% time of these trips 
% The waiting time is considered in calculation of travel time.
% Additionally if we consider trips as journeys between stations only then
% we add a 1.2 factor to account for the last leg of the trip which is
% destination to carpark.

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
% bookingFile = sprintf('/Users/katarzyna/Dropbox/matlab/2016-03-Demand_generation/boookings_ecbd_808356.txt');
bookingFile = sprintf('/home/kasia/Dropbox/matlab/2016-03-Demand_generation/boookings_ecbd_sorted_808356.txt');

bookingData = dlmread(bookingFile, ' ', 0, 0);

bookingTime = bookingData(:,2);
origX = bookingData(:,4);
origY = bookingData(:,5);
destX = bookingData(:,6);
destY = bookingData(:,7);

% station file
disp('2. Import stations...')
% facilityFile = sprintf('/Users/katarzyna/Dropbox/matlab/2016-03-Demand_generation/facility_location/stations_ecbd34.txt');
facilityFile = sprintf('/home/kasia/Dropbox/matlab/2016-03-Demand_generation/facility_location/stations_ecbd34.txt');
stationsData = dlmread(facilityFile, ' ', 0, 0);

f_ids = stationsData(:,1);
stationX = stationsData(:,2);
stationY = stationsData(:,3);

% nodes file
disp('3. Import nodes...')
%nodesFile = sprintf('/Users/katarzyna/Dropbox/matlab/2016-03-Demand_generation/facility_location/input/ecbd_nodes.csv');
nodesFile = sprintf('/home/kasia/Dropbox/matlab/2016-03-Demand_generation/input/ecbd_nodes.csv');
nodesData = dlmread(nodesFile, ',', 0, 0);

node_ids = nodesData(:,1);
nodeX = nodesData(:,2);
nodeY = nodesData(:,3);
%% for each booking, find distance from origin to destination
disp('4. Find distance from origins to destinations...')
O2Ddistance = zeros(length(origX),1); % in m
arrivalTime = int32(zeros(length(origX),1)); % in secs
traveltime = int32(zeros(length(origX),1)); % in secs

% we calculate distance in euclidian space, not based on simmobility
% network, so euclideanDist = true; false if based on simmobility network
euclideanDist = true;
if (euclideanDist)
    % do the euclidian distances between origin and destination
    % distances are increased by 37% due to euclidean distance (factor
    % based on Emilio's papar)
    correction_factor = 1.37; % maybe it should be a random value from an interval
    for i=1:length(origX)
        O2Ddistance(i) = correction_factor*pdist2([origX(i) origY(i)], [destX(i) destY(i)]);
    end
else
    % exectuted only if the distance is based on simmobility network
    %find the nearest node for each origin and destination and find the
    %distance between nodes from the shortest_distances output file
%     node_origin = pdist2([origX(:), origY(:)], [nodeX(:), nodeY(:)]);
%     node_dest = pdist2([destX(:), destY(:)], [nodeX(:), nodeY(:)]);
%     
%     costestNode_orig = zeros(length(origX),1);
%     for ii = 1:length(origX)
%         [minDist_, ind] = min(node_origin(ii,:));
%         costestNode_orig(ii) = ind;
%     end
%     
%     closestNode_dest = zeros(length(destX),1);
%     for ii = 1:length(destX)
%         [minDist_, ind] = min(node_dest(ii,:));
%         closestNode_dest(ii) = ind;
%     end
end

%% Calculate travel time and arrival time
% arrival time = request time + waiting time + travel time
% travel time = distance travelled*ave_speed

% waiting time is gamma distributed with a = 2 and b = 4;
% R = gamrnd(A,B) generates random numbers from the gamma distribution 
% with shape parameters in A and scale parameters in B. 
% If A is restricted to integers, the gamma distribution is referred to as 
% the Erlang distribution used in queueing theory. A = 1 -> equivalent to
% exponential distribution
% B controls the scale of the data. When it becomes large, the gamma 
% distribution approaches the normal distribution.

% speed is normally distributed with mean = 28.9 km/h (8.02 m/s) and
% st.dev = 10 km/h (2.78 m/s) and minimum speed is 3m/s (10.08 km/h)

disp('5. Calculate travel time...')
for i = 1:length(origX)
    
    waiting_time = gamrnd(2,4);
    speed_distr = normrnd(8.02, 2.78);
    speed = max(3, speed_distr);
    arrivalTime(i) = bookingTime(i) + O2Ddistance(i)/speed + waiting_time;
    traveltime(i) = O2Ddistance(i)/speed;
end

bookings_wATT = [int32(bookingTime), origX, origY, destX, destY, int32(traveltime), int32(arrivalTime)];
%% Save to file
disp('6. Save bookings with travel time...')

trips_ArrivalAndTT = sprintf('trips_ArrivalAndTT%dstations2016-03-21.txt', length(f_ids));
fileArrivals = fopen(trips_ArrivalAndTT,'w');

for j = 1:length(bookingTime)
    fprintf(fileArrivals,'%0u %f %f %f %f %0u %0u\n', bookingTime(j), origX(j), origY(j), destX(j), destY(j), traveltime(j), arrivalTime(j));
end
fclose(fileArrivals);

disp('All done.')

%% trips between the stations only
% find the nearest station for each origin and destiantion
disp('7. Find the nearest station for each origin and destiantion...')
closestStO = zeros(length(origX),1);
closestStD = zeros(length(origX),1);
ind_delete = zeros(length(origX),1);
for i = 1:length(origX)
    [o_distance, o_closestStationID ] = min(pdist2([origX(i) origY(i)], [stationX(:) stationY(:)]));
    [d_distance, d_closestStationID ] = min(pdist2([destX(i) destY(i)], [stationX(:) stationY(:)]));
    
    if (o_closestStationID ~=  d_closestStationID)
        closestStO(i) = o_closestStationID; % this gives the indices
        closestStD(i) = d_closestStationID; % this gives the indices
    else
        ind_delete(i) = i;
    end
end

ind_delete = find(ind_delete); % non zero indx 
% delete empty rows
closestStO(ind_delete) = []; 
closestStD(ind_delete) = [];
bookingTime(ind_delete) = [];
traveltime(ind_delete) = [];
arrivalTime(ind_delete) = [];

%% find node_id for each station

station_id_o = int32(zeros(length(closestStO),1));
station_id_d = int32(zeros(length(closestStD),1));

f_ids_ = int32(f_ids);

for i = 1:length(closestStO)
    station_id_o(i) = f_ids_(closestStO(i));
    station_id_d(i) = f_ids_(closestStD(i));
end

%% Save bookings
disp('8. Save bookings between stations...')
% new_bookings = [bookingTime, station_id_o, station_id_d, traveltime, arrivalTime];
% newBooking_f = sprintf('tripsBetween%dStation2016-03-21.txt', length(f_ids));
% delimiter = ' ';
% dlmwrite(newBooking_f, new_bookings, 'precision','%0u', delimiter);

newBooking_f = sprintf('tripsBetween%dStation2016-03-21.txt', length(f_ids));
fileBookings = fopen(newBooking_f,'w');

for j = 1:length(bookingTime)
    fprintf(fileBookings,'%0u %0u %0u %0u %0u\n', bookingTime(j), station_id_o(j), station_id_d(j), traveltime(j), arrivalTime(j));
end
fclose(fileBookings);

disp('All done.')