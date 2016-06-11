close all; clear; clc;
%% Convert rebalancing start time

sim_start_time = 6 *3600; % start which we want to have, in seconds
% default start time of the output file starts at 0.
end_day = 24*3600; % how many simulation hours

disp('1. Import rebalancing trips...')
bookingFile = sprintf('rebalancingCounts_withInit_per96_st34.txt');
bkData = dlmread(bookingFile, ' ', 0, 0);

reb_time = bkData(:,1);
fromNode = bkData(:,2);
toNode = bkData(:,3);
nOfTrips = bkData(:,4);

clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% find new time for each trip
disp('2. Find new time for each trip')
new_time = rem(reb_time + end_day - sim_start_time, end_day);
reb_interval = reb_time(1);

% shift all vectors circularly
index = -1;
for i = 1 : length(reb_time)
    if (reb_time(i) >= sim_start_time)
        index = i - 1;
        break;
    end
end

% sort rows
new_time = circshift(new_time, -index);
fromNode = circshift(fromNode, -index);
toNode = circshift(toNode, -index);
nOfTrips = circshift(nOfTrips, -index);

%% save to file with the a new start time
disp('3. Save to file.')
filename = sprintf('rebalancingCounts_withInit_per96_st34_start%dam.txt', sim_start_time/3600);
file = fopen(filename,'w');

% we add 900 seconds because first rebalancing is not at time zero, but 1
% rebalancing interval later
for j = 1:length(reb_time)
    fprintf(file,'%0u %0u %0u %0u\n', new_time(j) + reb_interval, fromNode(j),toNode(j), nOfTrips(j));
end
fclose(file);

disp('Done.')