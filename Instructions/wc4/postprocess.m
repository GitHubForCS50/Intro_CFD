osc_data_name = fullfile('outputs','2026-05-10_16-48-11','osc.dat');
data = readtable(osc_data_name);

figure;
yyaxis left;
plot(data.time,data.T,'yellow');
ylabel('Temperature');
yyaxis right;
plot(data.time,data.u,'r',data.time,data.v,'b');
ylabel('Velocities');
xlabel('Time');
legend(["T", "u", "v"]);
grid on;