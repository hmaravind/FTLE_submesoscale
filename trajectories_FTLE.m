clear; close all; clc;
%% Start and end times, depth
d = 2;                  % depth
tMin = 40; tMax = 80;   % minimum and maximum time (hours) for which data is available
nt = 6;                 % number of time windows with distinct tStart
flag = 0;               % 0 to use particle grid defined in function; 1 to define another

tStart = [linspace(tMin,tMax,nt+1) linspace(tMax,tMin,nt+1)];
tStart = tStart([1:end/2-1, end/2+1:end-1]);
tStop = [tMax*ones(1,nt) tMin*ones(1,nt)];
dt = [(tStart(2)-tStart(1))*ones(1,nt), (tStart(nt+2)-tStart(nt+1))*ones(1,nt)];

tStart = [tMin:0.25:tMax]; nt = length(tStart); tStart = [tStart fliplr(tStart)];
dt = [10*ones(1,nt) -10*ones(1,nt)];
tStop = tStart+dt;

dt = [1 2 3 4 5 6 -1 -2 -3 -4 -5 -6 -15];
tStop = [ones(1,length(dt)-1)*80 70];
tStart = tStop-dt;

tStart = 80; tStop = 70; dt = -10; d = 2;

tStep = 300*((tStop-tStart)/abs(tStop-tStart));
disp([tStart; tStop; dt; tStep]);

%% Trajectory extraction and FTLE computation
for i = 1:length(tStart)
    display(['tStart = ', num2str(tStart(i)), ', tStop = ', num2str(tStop(i)), ', d = ', num2str(d)]);
    trajectory_calculation_periodic(tStart(i),tStop(i),d,dt(i),tStep(i),flag);     % Trajectory extraction
    compute_FTLE(tStart(i),tStop(i),d,tStep(i),0);                        % FTLE computation
end

%% FTLE computation
% trajectoryFiles = dir('Trajectories_FTLE/depth*.mat');
% parfor i = 1:size(trajectoryFiles,1)
% %     display(['tStart = ', num2str(tStart(i)), ', tStop = ', num2str(tStop(i)), ', d = ', num2str(d)]);
%     fileName = trajectoryFiles(i).name; fileName = fileName(1:end-4)
%     j = 0; d = str2num(fileName(6)); 
%     if d == 1; j = 1; d = str2num(fileName(6:7)); end
%     tStart = str2num(fileName([11:12]+j)); tStop = str2num(fileName([14:15]+j)); 
%     compute_FTLE(tStart,tStop,d);
% end
