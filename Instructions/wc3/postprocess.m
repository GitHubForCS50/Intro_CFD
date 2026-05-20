% evaluate_hx.m
% Purpose: Read output.dat and calculate heat exchanger performance indicators

% 1. Load the data 
% readtable uses the headers we added to create an easily searchable table
data = readtable('output.dat');

% 2. Identify the inflow and outflow X-coordinates
% We use a small tolerance (1e-6) to avoid floating-point rounding errors
x_in = min(data.x);
x_out = max(data.x);

inflow_idx  = abs(data.x - x_in) < 1e-6;
outflow_idx = abs(data.x - x_out) < 1e-6;

% 3. Extract the variables at the inflow and outflow planes
T_in  = data.T(inflow_idx);
p_in  = data.p(inflow_idx);

T_out = data.T(outflow_idx);
p_out = data.p(outflow_idx);

% 4. Calculate the simple arithmetic averages
T_in_avg  = mean(T_in);
T_out_avg = mean(T_out);

p_in_avg  = mean(p_in);
p_out_avg = mean(p_out);

% 5. Calculate the performance indicators
% Note: Pressure drop is usually (In - Out), Temp increase is (Out - In)
delta_T = mean(T_out - T_in); 
delta_p = mean(p_out - p_in); 

% 6. Print the results to the console
fprintf('\n================================================\n');
fprintf('       HEAT EXCHANGER PERFORMANCE REPORT        \n');
fprintf('================================================\n');
fprintf('Temperature (T):\n');
fprintf('  Average Inflow:  %8.2f K\n', T_in_avg);
fprintf('  Average Outflow: %8.2f K\n', T_out_avg);
fprintf('  Temp Increase:   %8.2f K   <-- (Higher is better)\n', delta_T);
fprintf('  Improvement:     %8.2f %%\n\n', (delta_T/59.4167 - 1)*100);

fprintf('Relative Pressure (p):\n');
fprintf('  Average Inflow:  %8.2e Pa\n', p_in_avg);
fprintf('  Average Outflow: %8.2e Pa\n', p_out_avg);
fprintf('  Pressure Drop:   %8.2e Pa  <-- (Lower is better)\n', abs(delta_p));
fprintf('  Worse by:        %8.2f %%\n', (abs(delta_p)/8.4475e-04 - 1)*100);
fprintf('================================================\n\n');

% figure;
% plot(p_out);