
%% Preliminary tests
% Sample figure 1 - sea suface height for the entire domain;
varIndices = [1 4 5 8];
d = 1; NetCDF_direct_download_saveas_mat; a = out;
fig = figure; pcolor(a.lon_uv,a.lat_uv,a.zeta(:,:,1)'); shading interp;
colorbar; title('Sample - SSH');
xlabel('longitude'); ylabel('latitude');
print(fig,'Figures/sample_SSH.png','-dpng','-r400');

% ocean_time, after data for three consecutive days have been loaded into a, b and c
varIndices = [1];
d = 2; NetCDF_direct_download_saveas_mat; a = out;
d = 3; NetCDF_direct_download_saveas_mat; b = out;
d = 4; NetCDF_direct_download_saveas_mat; c = out;

fig = figure; hold on; plot(1:25,a,1:25,b,1:25,c);
plot(1:25,a(end)*ones(1,25),'k--',1:25,b(1)*ones(1,25),'k--',1:25,b(end)*ones(1,25),'k--',1:25,c(1)*ones(1,25),'k--');
xlabel('Index'); ylabel('ocean\_time (hours since 1968-05-23 00:00:00 GMT)');
title('SOCIB data: ocean\_time - 25 data points every 3 hours'); 
l = legend('2018/01/02','2018/01/03','2018/01/04'); l.Location = 'northwest';
print(fig,'Figures/ocean_time.png','-dpng','-r400');

% difference in vbar between data for the same time instances, but stored for different days - a is the day after b
varIndices = [1 12];
d = 4; NetCDF_direct_download_saveas_mat; a = out;
d = 3; NetCDF_direct_download_saveas_mat; b = out;

for k = 0:16; A(k+1) = min(min(a.vbar(:,:,1+k) == b.vbar(:,:,9+k))); B(k+1) = a.ocean_time(1+k)==b.ocean_time(9+k); end; [A;B]
for k = 0:16; figure; pcolor((a.vbar(:,:,1+k)-b.vbar(:,:,1+8+k))'); shading interp; colorbar; colormap('bluewhitered'); title(['k = ',num2str(k)]); end
a.vbar(isnan(a.vbar(:))) = 0; b.vbar(isnan(b.vbar(:))) = 0;
for k = 0:16; A(k+1) = min(min(a.vbar(:,:,1+k) == b.vbar(:,:,9+k))); B(k+1) = a.ocean_time(1+k)==b.ocean_time(9+k); end; [A;B]

%% Sample figure - sea suface height for the required domain
load('C:\Aravind_Research\SOCIB_Alboran Sea\data\2018090100.mat')
load('C:\Aravind_Research\SOCIB_Alboran Sea\data\grid_data.mat')
fig = figure; pcolor(lon_uv,lat_uv,zeta'); shading interp;
colorbar; title(['Sample - SSH on ', datestr(datetime(1968,5,23)+seconds(ocean_time),0)]);
xlabel('longitude'); ylabel('latitude');
print(fig,'Figures/sample_SSH.png','-dpng','-r400');

%% Sample figure - trajectories
load('C:\Aravind_Research\SOCIB_Alboran Sea\data_2018_09\grid_data.mat')
load('C:\Aravind_Research\SOCIB_Alboran Sea\Trajectories_FTLE\2018_09_b_2018091000_2018090500.mat')
fig1 = figure; scatter(x(1,:),y(1,:),5,'k.'); axis([min(lon_uv) max(lon_uv) min(lat_uv) max(lat_uv)]);
xlabel('longitude'); ylabel('latitude'); title('Initial particle distribution on 2018/09/10');
fig2 = figure; scatter(x(2,:),y(2,:),5,'k.'); axis([min(lon_uv) max(lon_uv) min(lat_uv) max(lat_uv)]);
xlabel('longitude'); ylabel('latitude'); title('Final particle distribution on 2018/09/5');
% print(fig1,'Figures/Particles_initial.png','-dpng','-r400');
% print(fig2,'Figures/Particles_final.png','-dpng','-r400');