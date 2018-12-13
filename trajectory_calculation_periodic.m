function [xArray, yArray, tArray] = trajectory_calculation_periodic(t0i,tFi,dt,tStep,flagParticles,flagAux,xParticles,yParticles,tStepOutput)
% flagParicles      0 to use grid defined in the function (for FTLE calculations)
%                   1 to input xParticles, yParticles to function

%% Select case, load time and grid data
folderLabel = '2018_09';
folderName = ['data_', folderLabel];
load(fullfile(folderName,'time_data.mat'));
load(fullfile(folderName,'grid_data.mat'));
tData = tData;
[xData, yData] = meshgrid(lon_uv,lat_uv);

%% Select time steps for advection
tStart = seconds(t0i-datetime(datetime(1968,5,23)));
tStop = seconds(tFi-datetime(datetime(1968,5,23)));
timeWindow = seconds(dt);

%% Select time step and output
% tStep = seconds(hours(3))/30*((tStop-tStart)/abs(tStop-tStart)); % adjusts tStep based on whether it is forward- or backward-time FTLE
if ~exist('tStepOutput','var'); tStepOutput = timeWindow; end % needs to be a multiple of tStep

%% Select particle grid
if flagParticles == 0
    xP = xData;    yP = yData;
%     dxA = min(xP(1,2)-xP(1,1),yP(2,1)-yP(1,1))/100;
%     dyA = dxA;
    XP = xP; YP = yP;
end

dxA = 1e-4;     dyA = 1e-4;
if ~exist('flagAux','var') 
    if flagParticles == 0; flagAux = 1;
    else; flagAux = 0; end
end

if ~exist('xParticles','var')
    if flagParticles == 0; xParticles = repmat(XP(:)',5,1);
    else; disp('Error: Initial particle locations not specified'); return; end
    if flagAux == 1
        xParticles(2,:) = xParticles(2,:) - dxA;
        xParticles(3,:) = xParticles(3,:) + dxA;
    end
end
if ~exist('yParticles','var')
    if flagParticles == 0; yParticles = repmat(YP(:)',5,1);
    else; disp('Error: Initial particle locations not specified'); return; end
    if flagAux == 1
        yParticles(4,:) = yParticles(4,:) - dyA;
        yParticles(5,:) = yParticles(5,:) + dyA;
    end
end
tic;
[xArray, yArray, tArray] = trajectory_folder(folderLabel, folderName, tStart, tStop, tStep, tStepOutput, xParticles, yParticles, flagAux, tData, xData, yData);
toc;

end

% Function that does what
function [xArray, yArray, tArray] = trajectory_folder(folderLabel, folderName, tStart, tStop, tStep, tStepOutput, xParticles, yParticles, flagAux, tData, xData, yData)
    
%     fprintf('Case %s started ...\n', folderLabel)
    
    if tStart<tStop; rt = 1; else rt = -1; end
    % rt = 1 for forward in time, -1 for backward in time
    destName = 'Trajectories';
    
    %% Initial Data Files
    
    d = dir(fullfile(folderName,'2018*.mat')); %%%%%%%%%%%
    folder = fullfile(pwd,folderName); %%%%%%%%%%%%%

    %% Skip files if we need it

%     fileSkipNumber = round(skipNumber/100);
%     d = d(1:fileSkipNumber:end);

    %% Set the time range and steps of interest
    
    stepSkip = tStepOutput/tStep; % so this is an integer
    tParticles = tStart:tStep:tStop;

    if rt == 1
        f2 = find(tData<tStart, 1, 'last');
    else
        f2 = find(tData>tStart, 1, 'first');
    end

    %% Set the advection particle distribution       
    
    xParticles = xParticles(:);
    yParticles = yParticles(:);

    tPos = zeros(1, 1+(length(tParticles)-1)/stepSkip);
    xPos = zeros(length(xParticles), length(tPos));
    yPos = zeros(length(yParticles), length(tPos));
    
    % Initial configuration
    xPos(:,1) = xParticles;
    yPos(:,1) = yParticles;
    tPos(1) = tStart;
    

    %% Load initial data distribution

    fprintf('t = %g / %g \n', abs(tParticles(1)-tStart), abs(tStart-tStop))
    
    %%%% coordinates
    x0 = 0; % min(unique(xData)); % get left value (periodic BC)
    Lx = 0; % max(unique(xData)) - x0; % get period
    R = 6.371e6; % Radius of the earth
    
    load(fullfile(folderName,d(f2-1*rt).name));
%     u(isnan(u)) = 0; v(isnan(v)) = 0;
    u1Data = u'./(R*cosd(yData));
    v1Data = v'./R;

    load(fullfile(folderName,d(f2).name));
%     u(isnan(u)) = 0; v(isnan(v)) = 0;
    u2Data = u'./(R*cosd(yData));
    v2Data = v'./R;

    load(fullfile(folderName,d(f2+1*rt).name));
%     u(isnan(u)) = 0; v(isnan(v)) = 0;
    u3Data = u'./(R*cosd(yData));
    v3Data = v'./R;
    tBreak = tData(f2+1*rt);

    load(fullfile(folderName,d(f2+2*rt).name));
%     u(isnan(u)) = 0; v(isnan(v)) = 0;
    u4Data = u'./(R*cosd(yData));
    v4Data = v'./R;

    lastInd = f2+2*rt;
    xCurrent = xPos(:,1);
    yCurrent = yPos(:,1);
    
    %% Loop through time steps advecting particles forward

    for i = 1:length(tParticles)-1
        
        fprintf('t = %g / %g \n',abs(tParticles(i+1)-tStart), abs(tStart-tStop))

        % Update data as necessary

        if tParticles(i)^rt > tBreak^rt

            lastInd = lastInd + 1*rt;

            u1Data = u2Data;
            v1Data = v2Data;
            
            u2Data = u3Data;
            v2Data = v3Data;
            
            u3Data = u4Data;
            v3Data = v4Data;

            load(fullfile(folder,d(lastInd).name));
%             u(isnan(u)) = 0; v(isnan(v)) = 0;
            u4Data = u'./(R*cosd(yData));
            v4Data = v'./R;
            tBreak = tData(lastInd-1*rt);
            
        end

        %% Calculate the velocity at the current, half step and next step - Lagranage Interpolating Polynomial order 3

        % Current step

        coef1 = (tParticles(i) - tData(lastInd-2*rt))*(tParticles(i) - tData(lastInd-1*rt))*(tParticles(i) - tData(lastInd))  /((tData(lastInd-3*rt) - tData(lastInd-2*rt))*(tData(lastInd-3*rt) - tData(lastInd-1*rt))*(tData(lastInd-3*rt) - tData(lastInd)))  ;
        coef2 = (tParticles(i) - tData(lastInd-3*rt))*(tParticles(i) - tData(lastInd-1*rt))*(tParticles(i) - tData(lastInd))  /((tData(lastInd-2*rt) - tData(lastInd-3*rt))*(tData(lastInd-2*rt) - tData(lastInd-1*rt))*(tData(lastInd-2*rt) - tData(lastInd)))  ;
        coef3 = (tParticles(i) - tData(lastInd-3*rt))*(tParticles(i) - tData(lastInd-2*rt))*(tParticles(i) - tData(lastInd))  /((tData(lastInd-1*rt) - tData(lastInd-3*rt))*(tData(lastInd-1*rt) - tData(lastInd-2*rt))*(tData(lastInd-1*rt) - tData(lastInd)))  ;
        coef4 = (tParticles(i) - tData(lastInd-3*rt))*(tParticles(i) - tData(lastInd-2*rt))*(tParticles(i) - tData(lastInd-1*rt))/((tData(lastInd)   - tData(lastInd-3*rt))*(tData(lastInd)   - tData(lastInd-2*rt))*(tData(lastInd)   - tData(lastInd-1*rt)));

        uCurrent = coef1*u1Data + coef2*u2Data + coef3*u3Data + coef4*u4Data;
        vCurrent = coef1*v1Data + coef2*v2Data + coef3*v3Data + coef4*v4Data;

        % Half step

        coef1 = (tParticles(i)+tStep/2 - tData(lastInd-2*rt))*(tParticles(i)+tStep/2 - tData(lastInd-1*rt))*(tParticles(i)+tStep/2 - tData(lastInd))  /((tData(lastInd-3*rt) - tData(lastInd-2*rt))*(tData(lastInd-3*rt) - tData(lastInd-1*rt))*(tData(lastInd-3*rt) - tData(lastInd)))  ;
        coef2 = (tParticles(i)+tStep/2 - tData(lastInd-3*rt))*(tParticles(i)+tStep/2 - tData(lastInd-1*rt))*(tParticles(i)+tStep/2 - tData(lastInd))  /((tData(lastInd-2*rt) - tData(lastInd-3*rt))*(tData(lastInd-2*rt) - tData(lastInd-1*rt))*(tData(lastInd-2*rt) - tData(lastInd)))  ;
        coef3 = (tParticles(i)+tStep/2 - tData(lastInd-3*rt))*(tParticles(i)+tStep/2 - tData(lastInd-2*rt))*(tParticles(i)+tStep/2 - tData(lastInd))  /((tData(lastInd-1*rt) - tData(lastInd-3*rt))*(tData(lastInd-1*rt) - tData(lastInd-2*rt))*(tData(lastInd-1*rt) - tData(lastInd)))  ;
        coef4 = (tParticles(i)+tStep/2 - tData(lastInd-3*rt))*(tParticles(i)+tStep/2 - tData(lastInd-2*rt))*(tParticles(i)+tStep/2 - tData(lastInd-1*rt))/((tData(lastInd)   - tData(lastInd-3*rt))*(tData(lastInd)   - tData(lastInd-2*rt))*(tData(lastInd)   - tData(lastInd-1*rt)));

        uHalf = coef1*u1Data + coef2*u2Data + coef3*u3Data + coef4*u4Data;
        vHalf = coef1*v1Data + coef2*v2Data + coef3*v3Data + coef4*v4Data;

        % Next step

        coef1 = (tParticles(i)+tStep - tData(lastInd-2*rt))*(tParticles(i)+tStep - tData(lastInd-1*rt))*(tParticles(i)+tStep - tData(lastInd))  /((tData(lastInd-3*rt) - tData(lastInd-2*rt))*(tData(lastInd-3*rt) - tData(lastInd-1*rt))*(tData(lastInd-3*rt) - tData(lastInd)))  ;
        coef2 = (tParticles(i)+tStep - tData(lastInd-3*rt))*(tParticles(i)+tStep - tData(lastInd-1*rt))*(tParticles(i)+tStep - tData(lastInd))  /((tData(lastInd-2*rt) - tData(lastInd-3*rt))*(tData(lastInd-2*rt) - tData(lastInd-1*rt))*(tData(lastInd-2*rt) - tData(lastInd)))  ;
        coef3 = (tParticles(i)+tStep - tData(lastInd-3*rt))*(tParticles(i)+tStep - tData(lastInd-2*rt))*(tParticles(i)+tStep - tData(lastInd))  /((tData(lastInd-1*rt) - tData(lastInd-3*rt))*(tData(lastInd-1*rt) - tData(lastInd-2*rt))*(tData(lastInd-1*rt) - tData(lastInd)))  ;
        coef4 = (tParticles(i)+tStep - tData(lastInd-3*rt))*(tParticles(i)+tStep - tData(lastInd-2*rt))*(tParticles(i)+tStep - tData(lastInd-1*rt))/((tData(lastInd)   - tData(lastInd-3*rt))*(tData(lastInd)   - tData(lastInd-2*rt))*(tData(lastInd)   - tData(lastInd-1*rt)));

        uNext = coef1*u1Data + coef2*u2Data + coef3*u3Data + coef4*u4Data;
        vNext = coef1*v1Data + coef2*v2Data + coef3*v3Data + coef4*v4Data;

        
        %% Interpolate data to get the velocity for RK4
        % select if 'linear', 'nearest', 'cubic', 'makima', or 'spline'
        method = 'linear';
% keyboard

        uK1 = interp2(xData,yData,uCurrent, x0 + mod((xCurrent-x0),Lx)            ,  yCurrent,              method);
        vK1 = interp2(xData,yData,vCurrent, x0 + mod((xCurrent-x0),Lx)            ,  yCurrent,              method);
        uK2 = interp2(xData,yData,uHalf   , x0 + mod((xCurrent+uK1*tStep/2-x0),Lx),  yCurrent+vK1*tStep/2,  method);
        vK2 = interp2(xData,yData,vHalf   , x0 + mod((xCurrent+uK1*tStep/2-x0),Lx),  yCurrent+vK1*tStep/2,  method);
        uK3 = interp2(xData,yData,uHalf   , x0 + mod((xCurrent+uK2*tStep/2-x0),Lx),  yCurrent+vK2*tStep/2,  method);
        vK3 = interp2(xData,yData,vHalf   , x0 + mod((xCurrent+uK2*tStep/2-x0),Lx),  yCurrent+vK2*tStep/2,  method);    
        uK4 = interp2(xData,yData,uNext   , x0 + mod((xCurrent+uK3*tStep-x0),Lx)  ,  yCurrent+vK3*tStep,    method);
        vK4 = interp2(xData,yData,vNext   , x0 + mod((xCurrent+uK3*tStep-x0),Lx)  ,  yCurrent+vK3*tStep,    method);

        xCurrent = xCurrent + tStep/6*(uK1 + 2*uK2 + 2*uK3 + uK4);
        yCurrent = yCurrent + tStep/6*(vK1 + 2*vK2 + 2*vK3 + vK4);

        uFinal = interp2(xData,yData,uNext   , x0 + mod((xCurrent-x0),Lx)  ,  yCurrent,    method);
        vFinal = interp2(xData,yData,vNext   , x0 + mod((xCurrent-x0),Lx)  ,  yCurrent,    method);
    
        %%% Add value to the output xPos, yPos, tPos;
        if mod(i,stepSkip)==0
            xPos(:,ceil((i)/stepSkip)+1) = xCurrent;
            yPos(:,ceil((i)/stepSkip)+1) = yCurrent;
            tPos(ceil((i)/stepSkip)+1) = tParticles(i+1);
        end

    end

    %% Remove points that went out of the domain
    indx = find(isnan(xPos(:,end)));
    indy = find(isnan(yPos(:,end)));
    indu = find(isnan(uFinal));
    indv = find(isnan(vFinal));

    fprintf('%.0f NaNs obtained\n',length(indu))
%     ind = find(~isnan(xPos(:,end)));
%     xPos = xPos(ind,:);
%     yPos = yPos(ind,:);
%     

    if ~exist(destName, 'dir'); mkdir(destName); end
    
    xArray = xPos';
    yArray = yPos';
    tArray = tPos';
    xLeft = x0;
    
    if flagAux == 1
        for i = 2:size(tArray,1)
        x = [xArray(1,:); xArray(i,:)]; y = [yArray(1,:); yArray(i,:)]; t = [tArray(1); tArray(i)];
        t1 = datetime(datetime(1968,5,23)+seconds(tArray(1)));
        t2 = datetime(datetime(1968,5,23)+seconds(tArray(i)));
        name1 = [num2str(year(t1)), sprintf('%02d',month(t1)), sprintf('%02d',day(t1)), sprintf('%02d',hour(t1))];
        name2 = [num2str(year(t2)), sprintf('%02d',month(t2)), sprintf('%02d',day(t2)), sprintf('%02d',hour(t2))];
        save(fullfile(destName,[folderLabel,'_',char('f'*(tStart<tStop)+'b'*(tStart>tStop)), '_',name1,'_',name2,'_',num2str(tStep),'.mat']),'x','y','t','xLeft','Lx','uFinal','vFinal','-v7.3');
        end
    else
        x = xArray; y = yArray; t = tArray;
        t1 = datetime(datetime(1968,5,23)+seconds(tStart));
        t2 = datetime(datetime(1968,5,23)+seconds(tStop));
        name1 = [num2str(year(t1)), sprintf('%02d',month(t1)), sprintf('%02d',day(t1)), sprintf('%02d',hour(t1))];
        name2 = [num2str(year(t2)), sprintf('%02d',month(t2)), sprintf('%02d',day(t2)), sprintf('%02d',hour(t2))];
        save(fullfile(destName,[folderLabel,'_',name1,'_',name2,'_',num2str(tStep),'.mat']),'x','y','t','xLeft','Lx','-v7.3');
    end
    disp('Trajectory data saved successfully.');
end