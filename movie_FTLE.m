% Movie of FTLE fields from t1 to t2 - both forward and backward in time
% with partilcles initialized at regions with high forward in time FTLE
% values at t1 advected with the flow till t2.

% d    - depth
% t1   - start time
% t2   - end time
% dt   - width of time window
function movie_FTLE(d,t1,t2,dt,percentageMax)
    if t1>t2; temp = t2; t2 = t1; t1 = temp; end   % ensure t1<t2
    tArray = t1:0.25:t2;
    
    % Start and end times for forward and backward time calculations
    t0Fwd = tArray; tfFwd = tArray + dt;
    t0Bwd = tArray; tfBwd = tArray - dt;
    
    name = ['depth',num2str(d),'m_',num2str(t1),'_',num2str(t2)];
    vidObj = VideoWriter(['Movies/FTLE_', name, '_', num2str(percentageMax), '_movie.mp4'],'MPEG-4');
    vidObj.Quality = 100; vidObj.FrameRate = 4; open(vidObj);
    
    % Load the grid data
    gridData = load(['depth',num2str(d),'m_vel/coordinates_time.mat']);
        
    for i = 1:length(tArray)
        display(['Iteration ', num2str(i), ' of ', num2str(length(tArray))]);
        
        % Load forward and backward time FTLE data, if they exist
        fileName_f = ['FTLE/depth',num2str(d),'m_','f','_',num2str(t0Fwd(i)),'_',num2str(tfFwd(i)),'.mat'];
        fileName_b = ['FTLE/depth',num2str(d),'m_','b','_',num2str(t0Bwd(i)),'_',num2str(tfBwd(i)),'.mat'];

        if(exist(fileName_f,'file')~=2 || exist(fileName_b,'file')~=2)
            disp('Error: FTLE data not found for the specified time window');
%             return;
        end
        fData = load(fileName_f);    bData = load(fileName_b);

        % Get velocity data using interpolation, calculate vorticity
        
        %% Find particles to be advected, calculate trajectories
        if i == 1
%             percentageMax = 0.6;
            xP = fData.XP0(fData.FTLE>=percentageMax*max(max(fData.FTLE)))'; % x positions of particles to be advected
            yP = fData.YP0(fData.FTLE>=percentageMax*max(max(fData.FTLE)))'; % y positions of particles to be advected
            trajectoryFile = ['Trajectories/depth',num2str(d),'m_',num2str(t1),'_',num2str(t2),'_',num2str(percentageMax),'.mat'];
            if ~exist(trajectoryFile,'file')
                disp('Trajectory file not found. Computing required trajectories ...')
                trajectory_calculation_periodic(t1,t2,d,tArray(2)-tArray(1),1,xP,yP);
                load(['Trajectories/depth',num2str(d),'m_',num2str(t1),'_',num2str(t2),'.mat']);
                save(trajectoryFile,'Lx','lxP','lyP','t','x','y','xLeft','percentageMax');
            end
            load(trajectoryFile); xArray = x; yArray = y;
            x0 = repmat(xArray(1,:),length(tArray),1);
            xArray = xLeft + mod((xArray-xLeft),Lx);
        end
        
        %% Plotting figures
        fig = figure; fig.Visible = 'off';
        hold on; colormap(fig,flipud(bone));
        fig.Position([3,4]) = fig.Position([3,4])*2;
        
        subplot(2,1,1); hold on;
        imagesc(fData.XP0(1,:),fData.YP0(:,1),fData.FTLE); 
        set(gca,'YDir','normal'); caxis([0 0.8]);
        c1 = colorbar; title(c1,'FTLE','Interpreter','Latex');
        ax1 = gca; ax2 = axes('Parent', fig, 'Position', ax1.Position);
        s=scatter(xArray(i,:), yArray(i,:),10,1:size(xArray,2),'filled');
        colormap(ax2,jet); s.MarkerFaceAlpha = 0.2;
        linkaxes([ax1,ax2]); axis([0 4098 -1000 1000]);
        ax2.Visible = 'off'; c = colorbar; c.Visible = 'off';
        xlabel(ax1,'$x_1$','Interpreter','Latex'); ylabel(ax1,'$x_2$','Interpreter','Latex');
        title(ax1,['Forward time: $t_0 = ',num2str(t0Fwd(i),'%0.2f'),'$, $t_f = ',num2str(tfFwd(i),'%0.2f'),'$, $d = ',num2str(d),'m$'], 'Interpreter','Latex');
        
        subplot(2,1,2); hold on;
        imagesc(bData.XP0(1,:),bData.YP0(:,1),bData.FTLE); 
        set(gca,'YDir','normal'); caxis([0 0.8]);
        c1 = colorbar; title(c1,'FTLE','Interpreter','Latex');
        ax1 = gca; ax2 = axes('Parent', fig, 'Position', ax1.Position);
        s=scatter(xArray(i,:), yArray(i,:),10,1:size(xArray,2),'filled');
        colormap(ax2,jet); s.MarkerFaceAlpha = 0.2;
        linkaxes([ax1,ax2]); axis([0 4098 -1000 1000]);
        ax2.Visible = 'off'; c = colorbar; c.Visible = 'off';
        xlabel(ax1,'$x_1$','Interpreter','Latex'); ylabel(ax1,'$x_2$','Interpreter','Latex');
        title(ax1,['Backward time: $t_0 = ',num2str(t0Bwd(i),'%0.2f'),'$, $t_f = ',num2str(tfBwd(i),'%0.2f'),'$, $d = ',num2str(d),'m$'], 'Interpreter','Latex');
        
        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',15);
        print(fig,['Movies/movie_FTLE_figures/',name,'_',num2str(percentageMax),'_',num2str(i),'.png'],'-dpng','-r400');
        if i == length(tArray); writeVideo(vidObj,getframe(fig)); writeVideo(vidObj,getframe(fig)); end
        writeVideo(vidObj,getframe(fig)); close(fig);
    end
    close(vidObj);
end