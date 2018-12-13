function FTLE = compute_FTLE(tStart,tStop,d,tStep,flagPlot)
    %% Loading up the data and reshaping for our needs
%     tStart = 42; tStop = 52; d = 2;
    fileName = ['depth',num2str(d),'m_',char('f'*(tStart<tStop)+'b'*(tStart>tStop)),'_',num2str(tStart),'_',num2str(tStop),'_',num2str(tStep)];
    load(['Trajectories_FTLE/',fileName, '.mat']);
    
    % Particles and auxiliary points
    x0 = reshape(x(1,:),5,size(x,2)/5);
    y0 = reshape(y(1,:),5,size(y,2)/5);
    xf = reshape(x(end,:),5,size(x,2)/5);
    yf = reshape(y(end,:),5,size(y,2)/5);
    uf = reshape(uFinal,5,size(uFinal,2)/5);
    vf = reshape(vFinal,5,size(vFinal,2)/5);
    
    % Remove points with NaNs in final x, y, u and v - points that crossed the domain
    for i = 1:5
        indTraj = unique([find((~isnan(xf(i,:)))), find((~isnan(yf(i,:)))), ...
            find((~isnan(uf(i,:)))), find((~isnan(vf(i,:))))], 'sorted');
        x0 = x0(:,indTraj);     xf = xf(:,indTraj);
        y0 = y0(:,indTraj);     yf = yf(:,indTraj);
    end
    fprintf('Number of points that crossed the domain: %.0f out of %.0f\n',length(x(1,:))/5-length(xf),length(x(1,:))/5);
    
    % Particle locations
    XP0 = x0(1,:);
    YP0 = y0(1,:);
    
    % Finite differencing for flow map gradient
    dx0 = x0(3,:) - x0(2,:);
    dy0 = y0(5,:) - y0(4,:);
    dxfh = xf(3,:) - xf(2,:);
    dxfv = xf(5,:) - xf(4,:);
    dyfh = yf(3,:) - yf(2,:);
    dyfv = yf(5,:) - yf(4,:);
    
    % Flow map gradient
    f11 = dxfh./dx0; f12 = dxfv./dy0;
    f21 = dyfh./dx0; f22 = dyfv./dy0;
    
    % Cauchy-Green deformation tensor
    c11 = f11.*f11 + f21.*f21;
    c12 = f11.*f12 + f21.*f22;
    c21 = f11.*f12 + f21.*f22; % c21 = c12
    c22 = f12.*f12 + f22.*f22;
    
    % Eigenvalues and FTLE
    timeInterval = tStop-tStart; eigs = zeros(2,length(c11));
    %% Eigenvalues from closed from expression
    %{
    eig1 = reshape(((c11+c22)+sqrt((c11+c22).^2-4*(c11.*c22-c12.*c21)))/2,lyP,lxP);
    eig2 = reshape(((c11+c22)-sqrt((c11+c22).^2-4*(c11.*c22-c12.*c21)))/2,lyP,lxP);
    %}
    %% Eigenvalues using a for loop
    for i = 1:length(c11)
        eigs(1:2,i) = eig([c11(i), c12(i); c21(i), c22(i)]);
    end
    %% Eigenvalues using a cell array
    %{
    len = length(c11);
    A(1,1,1:len) = c11; A(1,2,1:len) = c12;
    A(2,1,1:len) = c21; A(2,2,1:len) = c22;
    c = arrayfun(@(k){A(:,:,k)}, 1:len);
    eig_cell = cellfun(@eig,c,'un',0);
    eigs = cell2mat(eig_cell);
    %}
    %%
    % Identify points with non-positive eigenvalues - as a result of numerial errors
    eig1 = eigs(1,:);           eig2 = eigs(2,:);
    indEig = find(eig1>0 & eig2>0);
    eig1pos = eig1(indEig);    eig2pos = eig2(indEig);
    fprintf('Number of points with non-positive eigenvalues: %.0f out of %.0f\n',length(eig1)-length(indEig),length(eig1));
    
    % FTLE and dilation rate, and corresponding grids
    FTLE = log(max(eig1,eig2))/(2*abs(timeInterval));
    DilR = log(eig1pos.*eig2pos)/abs(timeInterval);
    xFTLE = XP0;                yFTLE = YP0;
    xDilR = XP0(indEig);        yDilR = YP0(indEig);
    
    % Plot data, and save to disk if needed
    dirData = 'FTLE'; dirFig = 'Figures';
    if ~exist(dirData, 'dir'); mkdir(dirData); end
    if ~exist(dirFig, 'dir'); mkdir(dirFig); end
    
    if flagPlot == 1
    fig1 = figure; scatter(xFTLE,yFTLE,5,FTLE,'filled'); colorbar; caxis([0 max(max(FTLE))]);
    axis tight; xlabel('$x_1$','Interpreter','Latex'); ylabel('$x_2$','Interpreter','Latex');
    title(['FTLE: $t_0 = ',num2str(tStart),'$, $t_f = ',num2str(tStop),'$, $d = ',num2str(d),'m$'], 'Interpreter','Latex');

    fig2 = figure; scatter(xDilR,yDilR,5,DilR,'filled'); colorbar; colormap('bluewhitered');
    axis tight; xlabel('$x_1$','Interpreter','Latex'); ylabel('$x_2$','Interpreter','Latex');
    title(['DilR: $t_0 = ',num2str(tStart),'$, $t_f = ',num2str(tStop),'$, $d = ',num2str(d),'m$'], 'Interpreter','Latex');
    
    print(fig1,'-r400','-dpng',fullfile(dirFig, [fileName, '.png'])); close(fig1);
    print(fig2,'-r400','-dpng',fullfile(dirFig, [fileName, '.png'])); close(fig2);
    end
    save(fullfile(dirData, [fileName, '.mat']),'XP0','YP0','xFTLE','yFTLE','xDilR','yDilR','eig1','eig2','FTLE','DilR','timeInterval');
    disp('FTLE data saved successfully.');
end