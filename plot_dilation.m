t0 = 70; tf = 80;
FTLE = load(['C:\Aravind_Research\UCSD_2018_08_31\FTLE\depth2m_b_',num2str(tf),'_',num2str(t0),'.mat']);
Traj = load(['C:\Aravind_Research\UCSD_2018_08_31\Trajectories_FTLE\depth2m_f_',num2str(t0),'_',num2str(tf),'.mat']);

badEig1Ind = FTLE.eig1<=0; badEig2Ind = FTLE.eig2<=0;
badEig1 = FTLE.eig1(badEig1Ind); badEig2 = FTLE.eig2(badEig2Ind);
nBadEig1 = length(badEig1); nBadEig2 = length(badEig2);
FTLE.eig1(badEig1Ind) = NaN; FTLE.eig2(badEig2Ind) = NaN;
Traj.xNew = Traj.xLeft + mod((Traj.x-Traj.xLeft),Traj.Lx);
dilation = log(FTLE.eig1.*FTLE.eig2)/abs(FTLE.timeInterval);

fig = figure; fig.Position([3 4]) = fig.Position([3 4])*2.5; fig.Position([1 2]) = 0;

subplot(2,1,1); ax1 = gca; pcolor(FTLE.XP0,FTLE.YP0,FTLE.FTLE); shading interp;
% subplot(2,1,1); ax1 = gca; pcolor(FTLE.XP0,FTLE.YP0,dilation); shading interp;
colormap(ax1,bluewhitered); colorbar;
% caxis([-2 1]); colormap(ax1,bluewhitered); colorbar;
% colormap(ax1,flipud(gray)); %caxis([0 max(max(dilation))]);
% hold on; scatter(FTLE.XP0((badEig1Ind+badEig2Ind)>0),FTLE.YP0((badEig1Ind+badEig2Ind)>0),10,'b');
hold on; s = scatter(Traj.xNew(2,:),Traj.y(2,:),5,'k.'); s.MarkerEdgeAlpha = 0.01;
% xlabel('x'); ylabel('y'); title(['Dilation for backward calculation from t = ',num2str(tf),' to ',num2str(t0), ' with points that give negative eigenvalues for CG tensor (o)']);
xlabel('x'); ylabel('y'); title(['Dilation for backward calculation from t = ',num2str(tf),' to ',num2str(t0)]);
xlabel('x'); ylabel('y'); title(['FTLE for backward calculation from t = ',num2str(tf),' to ',num2str(t0)]);
xlabel('x'); ylabel('y'); title(['FTLE (backward) and trajectories (forward) at t = ',num2str(tf)]);


subplot(2,1,2); ax2 = gca; pcolor(FTLE.XP0,FTLE.YP0,dilation); shading interp;
% subplot(2,1,2); ax1 = gca; pcolor(FTLE.XP0,FTLE.YP0,FTLE.FTLE); shading interp;
caxis([-2 1]); colormap(ax2,bluewhitered); colorbar; 
% colormap(ax1,bluewhitered); colorbar;

% colormap(ax2,flipud(gray)); %caxis([0 max(max(dilation))]);
hold on; s = scatter(Traj.xNew(2,:),Traj.y(2,:),5,'k.'); s.MarkerEdgeAlpha = 0.01;
% xlabel('x'); ylabel('y'); title(['FTLE (backward) and trajectories (forward) at t = ',num2str(tf)]);
xlabel('x'); ylabel('y'); title(['Dilation (backward) and trajectories (forward) at t = ',num2str(tf)]);