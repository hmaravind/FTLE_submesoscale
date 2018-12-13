clear; close all; clc;
%% Initializations
% Time for which data is available
% Data is currently available for y = 2018; m = 9; d = 1:30; H = 0:3:21;
tMin = datetime([2018,09,01,00,00,00]);
tMax = datetime([2018,09,30,21,00,00]);   

% Input time-window for analysis, start time (tStart) and end time (tStop).
% tStart and tStop - format: [yyyy,mm,dd,HH,MM,SS]; should lie between tMin and tMax
timeWindow = -days(5);
% tStart = datetime([2018,09,10,00,00,00]);

y = 2018; m = 9; d = 10:1:20; H = 0;
tStart = repmat(datetime([0,0,0]),(length(y)*length(m)*length(d)*length(H)),1);
for i = 1:length(d)
    tStart(i) = datetime([y,m,d(i),H,0,0]);
end
tStop = tStart + timeWindow;
i = 1; tStart = tStart(i); tStop = tStop(i);
if(min(min([tStart;tStop]>tMax)) || min(min([tStart;tStop]<tMin)))
    disp('No data for the times specified'); return
else
    disp([tStart,tStop]')
end

flag = 0;               % 0 to use particle grid defined in function; 1 to define another
tStep = seconds(hours(3))/30*((tStop-tStart)/abs(tStop-tStart));

%% Trajectory extraction and FTLE computation
for i = 1:length(tStart)
    disp(['tStart = ',datestr(tStart(i),0), ', tStop = ', datestr(tStop(i),0)]);
    trajectory_calculation_periodic(tStart(i),tStop(i),timeWindow,tStep(i),flag);      % Trajectory extraction
    compute_FTLE(tStart(i),tStop(i),tStep,0);                                    % FTLE computation
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
