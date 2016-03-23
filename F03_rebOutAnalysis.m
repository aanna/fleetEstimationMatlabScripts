close all; clear; clc;

% have to provide as an input the number of stations and number of
% rebalancing intervals (it can't be recognized from the output file)
n_stations = 34;
n_reb_periods = 96;
vec_for_reb = 1:n_stations*n_stations;
reb_matrix = transpose(reshape(vec_for_reb, [n_stations, n_stations]));

%% Import list of nodes within the analysed zone
disp('1. Import output from Gurobi cpp...')
filename = '/home/kasia/Documents/rebalancing_cpp/rebalancingMethods/rebalancing_offline/solution_rebalancing.sol';
delimiter = ',';
startRow = 3;
formatSpec = '%s%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

var_name = dataArray{:, 1};
time_ = dataArray{:, 2};
depSt_orSt = dataArray{:, 3}; % depending on thr var_name, i.e., if reb_veh 
%then this is departing station, if available_veh then @ which station
arrSt_orSt = dataArray{:, 4}; % same as above
quantity = dataArray{:, 5};

clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% Generate 3 matrices:
disp('2. Generate 3 matrices...')
% 1) vehicles available at each station for all time intervals
available_veh = (strcmp(var_name,'vhs_st_i'));
available_veh_m = zeros(n_reb_periods, n_stations);
% 2) rebalancing counts between stations for all time intervals
reb_veh = (strcmp(var_name,'nEmptyVhsTime'));
reb_veh_m = zeros(n_reb_periods, n_stations*n_stations);
% 3) moving vehicles for all time intervals
veh_in_transit = (strcmp(var_name,'in_transit'));
veh_in_transit_m = zeros(n_reb_periods, 1);

for i =1 : length (var_name)
   if (available_veh(i))
       available_veh_m(time_(i)+1, arrSt_orSt(i)+1) = quantity(i);
   end
   
   if (reb_veh(i))
       column_indx = reb_matrix(depSt_orSt(i)+1, arrSt_orSt(i)+1);
       reb_veh_m(time_(i)+1, column_indx) = quantity(i);
   end  
   
   if (veh_in_transit(i))
       veh_in_transit_m(time_(i)+1,1) = quantity(i);
   end 
end

%% save to file
disp('9. Save customer file...')
% filenameC = sprintf('customers_ecbd_%d.txt', length(origin_x));
% fileCustomers = fopen(filenameC,'w');
%
% for j = 1:length(origin_x)
%     fprintf(fileCustomers,'%0u %0f %0f\n', j, origin_x(j), origin_y(j));
% end
% fclose(fileCustomers);
%
% disp('10. Save booking file...')
% amod_mode = 1; % mode = 1 if this is amod trip
% filenameB = sprintf('boookings_ecbd_%d.txt', length(origin_x));
% fileBookings = fopen(filenameB,'w');
%
% for j = 1:length(origin_x)
%     fprintf(fileBookings,'%0u %0u %0u %0f %0f %0f %0f %0u\n', j, time_sec(j), j, origin_x(j), origin_y(j), dest_x(j), dest_y(j), amod_mode);
% end
% fclose(fileBookings);

disp('All done.')