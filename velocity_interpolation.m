% Script
function velocity_interpolation

%% Select case
folderLabel = 'depth2m';

%% Select time steps for advection
tS = 60*3600;
tF = 70*3600;

%% Select time step and output
tStep = 300;
tStepOutput = 300; % needs to be a multiple of tStep

%% Select grid
dxP = 2;
dyP = 2;
[xParticles, yParticles] = meshgrid(0:dxP:600, 0:dyP:600);

tic;
velocity_folder(folderLabel, tS, tF, tStep, tStepOutput, xParticles, yParticles);
toc;

end

% Function that does what
function velocity_folder(folderLabel, tS, tF, tStep, tStepOutput, xParticles, yParticles)  
    
    fprintf('Case %s started ...\n', folderLabel)
    
    folderName = sprintf('%s_vel',folderLabel);
    destName = sprintf('%s_vel_interp',folderLabel);
    
    %% Set the parameters for the files to be used

%     dataFileTimeStep = 5e-4;
%     skipNumber = 100;

    %% Initial Data Files

    % load tData
    load(fullfile(folderName,'coordinates_time.mat'), 'tData');
    
    d = dir(fullfile(folderName,'file_*.mat')); %%%%%%%%%%%
    folder = fullfile(pwd,folderName); %%%%%%%%%%%%%

%     step_S = str2double(d(1).name(end-9:end-4));
%     step_F = str2double(d(end).name(end-9:end-4));
% 
%     tData = step_S:skipNumber:step_F;
%     tData = tData * dataFileTimeStep;

    %% Skip files if we need it

%     fileSkipNumber = round(skipNumber/100);
%     d = d(1:fileSkipNumber:end);

    %% Set the time range and steps of interest
    
    stepSkip = tStepOutput/tStep; % so this is an integer
    tParticles = tS:tStep:tF;

    f2 = find(tData<tS, 1, 'last');

    %% Set the advection particle distribution       

    tPos = zeros(1, 1+(length(tParticles)-1)/stepSkip);
%     xPos = zeros(length(xParticles), length(tPos));
%     yPos = zeros(length(yParticles), length(tPos));
    
    % Initial configuration
    xPos = xParticles(:);
    yPos = yParticles(:);
    tPos(1) = tS;
    

    %% Load initial data distribution

    fprintf('t = %g / %g \n', tParticles(1), tF)
    
    %%%% coordinates
    load(fullfile(folderName,'coordinates_time.mat'));
    xData = x1G;
    yData = x2G;
    x0 = min(unique(xData)); % get left value (periodic BC)
    Lx = max(unique(xData)) - x0; % get period

    load(fullfile(folderName,d(f2-1).name));
    u1Data = u1G;
    v1Data = u2G;
%     w1Data = u3G;

    load(fullfile(folderName,d(f2).name));
    u2Data = u1G;
    v2Data = u2G;
%     w2Data = u3G;

    load(fullfile(folderName,d(f2+1).name));
    u3Data = u1G;
    v3Data = u2G;
%     w3Data = u3G;
    tBreak = tData(f2+1);

    load(fullfile(folderName,d(f2+2).name));
    u4Data = u1G;
    v4Data = u2G;
%     w4Data = u3G;

    lastInd = f2+2;
    xCurrent = xPos;
    yCurrent = yPos;
    

    %% Loop through time steps advecting particles forward

    for i = 1:length(tParticles)-1
        
        fprintf('t = %g / %g \n',tParticles(i+1), tF)

        % Update data as necessary

        if tParticles(i) > tBreak % tStep needs to be small so that update is precise

            lastInd = lastInd + 1;

            u1Data = u2Data;
            v1Data = v2Data;
%             w1Data = w2Data;
            
            u2Data = u3Data;
            v2Data = v3Data;
%             w2Data = w3Data;
            
            u3Data = u4Data;
            v3Data = v4Data;
%             w3Data = w4Data;

            load(fullfile(folder,d(lastInd).name));
            u4Data = u1G;
            v4Data = u2G;
%             w4Data = u3G;
            tBreak = tData(lastInd-1); %%%% THIS LINE WAS WRONG, I HAD lastInd-2 INSTEAD OF lastInd-1
            
        end
        
        if mod(i,stepSkip)~=0
            continue
        end

        %% Calculate the velocity at the current, half step and next step

%         % Current step
% 
%         coef1 = (tParticles(i) - tData(lastInd-2))*(tParticles(i) - tData(lastInd-1))*(tParticles(i) - tData(lastInd))  /((tData(lastInd-3) - tData(lastInd-2))*(tData(lastInd-3) - tData(lastInd-1))*(tData(lastInd-3) - tData(lastInd)))  ;
%         coef2 = (tParticles(i) - tData(lastInd-3))*(tParticles(i) - tData(lastInd-1))*(tParticles(i) - tData(lastInd))  /((tData(lastInd-2) - tData(lastInd-3))*(tData(lastInd-2) - tData(lastInd-1))*(tData(lastInd-2) - tData(lastInd)))  ;
%         coef3 = (tParticles(i) - tData(lastInd-3))*(tParticles(i) - tData(lastInd-2))*(tParticles(i) - tData(lastInd))  /((tData(lastInd-1) - tData(lastInd-3))*(tData(lastInd-1) - tData(lastInd-2))*(tData(lastInd-1) - tData(lastInd)))  ;
%         coef4 = (tParticles(i) - tData(lastInd-3))*(tParticles(i) - tData(lastInd-2))*(tParticles(i) - tData(lastInd-1))/((tData(lastInd)   - tData(lastInd-3))*(tData(lastInd)   - tData(lastInd-2))*(tData(lastInd)   - tData(lastInd-1)));
% 
%         uCurrent = coef1*u1Data + coef2*u2Data + coef3*u3Data + coef4*u4Data;
%         vCurrent = coef1*v1Data + coef2*v2Data + coef3*v3Data + coef4*v4Data;
%         wCurrent = coef1*w1Data + coef2*w2Data + coef3*w3Data + coef4*w4Data;

%         % Half step
% 
%         coef1 = (tParticles(i)+tStep/2 - tData(lastInd-2))*(tParticles(i)+tStep/2 - tData(lastInd-1))*(tParticles(i)+tStep/2 - tData(lastInd))  /((tData(lastInd-3) - tData(lastInd-2))*(tData(lastInd-3) - tData(lastInd-1))*(tData(lastInd-3) - tData(lastInd)))  ;
%         coef2 = (tParticles(i)+tStep/2 - tData(lastInd-3))*(tParticles(i)+tStep/2 - tData(lastInd-1))*(tParticles(i)+tStep/2 - tData(lastInd))  /((tData(lastInd-2) - tData(lastInd-3))*(tData(lastInd-2) - tData(lastInd-1))*(tData(lastInd-2) - tData(lastInd)))  ;
%         coef3 = (tParticles(i)+tStep/2 - tData(lastInd-3))*(tParticles(i)+tStep/2 - tData(lastInd-2))*(tParticles(i)+tStep/2 - tData(lastInd))  /((tData(lastInd-1) - tData(lastInd-3))*(tData(lastInd-1) - tData(lastInd-2))*(tData(lastInd-1) - tData(lastInd)))  ;
%         coef4 = (tParticles(i)+tStep/2 - tData(lastInd-3))*(tParticles(i)+tStep/2 - tData(lastInd-2))*(tParticles(i)+tStep/2 - tData(lastInd-1))/((tData(lastInd)   - tData(lastInd-3))*(tData(lastInd)   - tData(lastInd-2))*(tData(lastInd)   - tData(lastInd-1)));
% 
%         uHalf = coef1*u1Data + coef2*u2Data + coef3*u3Data + coef4*u4Data;
%         vHalf = coef1*v1Data + coef2*v2Data + coef3*v3Data + coef4*v4Data;
% 
%         % Next step
% 
        coef1 = (tParticles(i)+tStep - tData(lastInd-2))*(tParticles(i)+tStep - tData(lastInd-1))*(tParticles(i)+tStep - tData(lastInd))  /((tData(lastInd-3) - tData(lastInd-2))*(tData(lastInd-3) - tData(lastInd-1))*(tData(lastInd-3) - tData(lastInd)))  ;
        coef2 = (tParticles(i)+tStep - tData(lastInd-3))*(tParticles(i)+tStep - tData(lastInd-1))*(tParticles(i)+tStep - tData(lastInd))  /((tData(lastInd-2) - tData(lastInd-3))*(tData(lastInd-2) - tData(lastInd-1))*(tData(lastInd-2) - tData(lastInd)))  ;
        coef3 = (tParticles(i)+tStep - tData(lastInd-3))*(tParticles(i)+tStep - tData(lastInd-2))*(tParticles(i)+tStep - tData(lastInd))  /((tData(lastInd-1) - tData(lastInd-3))*(tData(lastInd-1) - tData(lastInd-2))*(tData(lastInd-1) - tData(lastInd)))  ;
        coef4 = (tParticles(i)+tStep - tData(lastInd-3))*(tParticles(i)+tStep - tData(lastInd-2))*(tParticles(i)+tStep - tData(lastInd-1))/((tData(lastInd)   - tData(lastInd-3))*(tData(lastInd)   - tData(lastInd-2))*(tData(lastInd)   - tData(lastInd-1)));
% 
        uNext = coef1*u1Data + coef2*u2Data + coef3*u3Data + coef4*u4Data;
        vNext = coef1*v1Data + coef2*v2Data + coef3*v3Data + coef4*v4Data;
%         wNext = coef1*w1Data + coef2*w2Data + coef3*w3Data + coef4*w4Data;


        % Linear in time
        
%         fprintf('tNow is %g, t2 and t3 are : %g and %g\n', tParticles(i) , tData(lastInd-2), tData(lastInd-1) )
%         
%         coef = (tParticles(i) - tData(lastInd-2))/(tData(lastInd-1) - tData(lastInd-2));
%         uCurrent = u2Data + coef*(u3Data-u2Data);
%         vCurrent = v2Data + coef*(v3Data-v2Data);
%         wCurrent = w2Data + coef*(w3Data-w2Data);
        
        
        %% Interpolate data to get the velocity for RK4
        % select if 'linear', 'nearest', 'cubic', 'makima', or 'spline'
        method = 'cubic';

        uK1 = interp2(xData,yData,uNext, x0 + mod((xCurrent-x0),Lx)            ,  yCurrent,              method);
        vK1 = interp2(xData,yData,vNext, x0 + mod((xCurrent-x0),Lx)            ,  yCurrent,              method);
%         wK1 = interp2(xData,yData,wNext, x0 + mod((xCurrent-x0),Lx)            ,  yCurrent,              method);
%         uK2 = interp2(xData,yData,uHalf   , x0 + mod((xCurrent+uK1*tStep/2-x0),Lx),  yCurrent+vK1*tStep/2,  method);
%         vK2 = interp2(xData,yData,vHalf   , x0 + mod((xCurrent+uK1*tStep/2-x0),Lx),  yCurrent+vK1*tStep/2,  method);
%         uK3 = interp2(xData,yData,uHalf   , x0 + mod((xCurrent+uK2*tStep/2-x0),Lx),  yCurrent+vK2*tStep/2,  method);
%         vK3 = interp2(xData,yData,vHalf   , x0 + mod((xCurrent+uK2*tStep/2-x0),Lx),  yCurrent+vK2*tStep/2,  method);    
%         uK4 = interp2(xData,yData,uNext   , x0 + mod((xCurrent+uK3*tStep-x0),Lx)  ,  yCurrent+vK3*tStep,    method);
%         vK4 = interp2(xData,yData,vNext   , x0 + mod((xCurrent+uK3*tStep-x0),Lx)  ,  yCurrent+vK3*tStep,    method);
% 
%         xCurrent = xCurrent + tStep/6*(uK1 + 2*uK2 + 2*uK3 + uK4);
%         yCurrent = yCurrent + tStep/6*(vK1 + 2*vK2 + 2*vK3 + vK4);

        

        %%% Add value to the output xPos, yPos, tPos;
        if mod(i,stepSkip)==0
%             xPos(:,ceil((i)/stepSkip)+1) = xCurrent;
%             yPos(:,ceil((i)/stepSkip)+1) = yCurrent;
            uPos(:,ceil((i)/stepSkip)+1) = uK1;
            vPos(:,ceil((i)/stepSkip)+1) = vK1;
%             wPos(:,ceil((i)/stepSkip)+1) = wK1;
            tPos(ceil((i)/stepSkip)+1) = tParticles(i+1);
            
%             fprintf('max u3 is %.5f, after interpolation: %.5f\n', max(abs(wCurrent(:))) , max(abs(wK1)) )
            
        end

    end

    %% Remove points that went out of the domain
    ind = find(isnan(xPos(:,end)));
    fprintf('%.0f NaNs obtained\n',length(ind))
    % ind = find(~isnan(xPos(:,end)));
    % xPos = xPos(ind,:);
    % yPos = yPos(ind,:);

    if ~exist(destName, 'dir'); mkdir(destName); end
    
    t = tPos';
    xLeft = x0;
    
    xGrid = reshape(xPos,size(xParticles));
    yGrid = reshape(yPos,size(xParticles));
    uGrid = reshape(uPos,[size(xParticles), length(t)]);
    vGrid = reshape(vPos,[size(xParticles), length(t)]);
%     wGrid = reshape(wPos,[size(xParticles), length(t)]);
    
%     save(fullfile(destName,'interpolated_velocity.mat'),'t','xGrid','yGrid','uGrid','vGrid','wGrid','xLeft','Lx','-v7.3');
    save(fullfile(destName,'interpolated_velocity.mat'),'t','xGrid','yGrid','uGrid','vGrid','xLeft','Lx','-v7.3');

end
