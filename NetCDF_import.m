clear; close all; clc;

allMonths = [1 2 3 4 5 6 7 8 9 10 11 12];
maxDays = [31 28 31 30 31 30 31 31 30 31 30 31];    % days in each month of the year specified above

year = 2018; months = 9;
for i = 1:length(months); days{months(i)} = 1:maxDays(months(i)); end
varIndices = [1 6:12]; 
xIndex = 5; yIndex = 4;
st = 8; sx = 250; sy = 180; dt = 3;                 % sizes of ocean_time, lon_uv and lat_uv required, tStep between two data points in hours
folderName = 'data'; infoFile = fullfile(folderName,'sample_info.txt');

%% Create a folder, save grid data, and a sample_info.txt file
if ~exist(folderName,'dir'); mkdir(folderName); end
month = 1; day = 1; nc_fileName = ['http://thredds.socib.es/thredds/dodsC/operational_models/oceanographical/hydrodynamics/wmop_surface/',num2str(year),'/',sprintf('%02d',month),'/roms_wmop_surface_',num2str(year), sprintf('%02d',month), sprintf('%02d',day),'.nc'];
diary(infoFile); ncdisp(nc_fileName); diary off;
nc_info = ncinfo(nc_fileName); xName = nc_info.Variables(xIndex).Name; yName = nc_info.Variables(yIndex).Name;
data_xy.(xName) = ncread(nc_fileName,xName); data_xy.(yName) = ncread(nc_fileName,yName);
data_xy.(xName) = data_xy.(xName)(1:sx); data_xy.(yName) = data_xy.(yName)(1:sy);
[data_xy.xGrid, data_xy.yGrid, data_xy.xUnif, data_xy.yUnif] = xyfromlatlon(data_xy.(yName),data_xy.(xName));

save(fullfile(folderName,'grid_data.mat'),'-struct','data_xy');
disp('Saved grid_data.mat');

%% Downloading and saving data
parObj = parpool(8);                                % number of cores to use - only if parfor is used in line 18
for month = months
    parfor day = 1:maxDays(month) %days{month}      % use parfor while downloading data for more than a few days - change lines 9 and 36       
        % Open NetCDF file
        nc_fileName = ['http://thredds.socib.es/thredds/dodsC/operational_models/oceanographical/hydrodynamics/wmop_surface/',num2str(year),'/',sprintf('%02d',month),'/roms_wmop_surface_',num2str(year), sprintf('%02d',month), sprintf('%02d',day),'.nc'];
        ncid = netcdf.open(nc_fileName,'NOWRITE');
        
        % Import data, reshape as per sizes specified and save
        disp(['Importing data for ', num2str(year), '/', sprintf('%02d',month), '/', sprintf('%02d',day), ' ...'])
        saveData(nc_fileName,folderName,varIndices,st,sx,sy,year,month,day,dt);
        
        netcdf.close(ncid); 
    end
end
delete(parObj)                                      % only if parfor is used in line 18

%% Saving time data in tData.mat
i = 0;
for month = months
    for day = days{month}
        for t = 1:st
            fileName = [num2str(year), sprintf('%02d',month), sprintf('%02d',day), sprintf('%02d',(t-1)*dt), '.mat'];
            load(fullfile(folderName,fileName),'ocean_time');
            i = i+1; tData(i) = ocean_time;
        end
    end
end
save(fullfile(folderName,'time_data.mat'),'tData');
disp('Saved time_data.mat');

%% Function to import data into a structure containing the same variables as in the NetCDF file, reshape it as per requirements and then save as .mat
function saveData(nc_fileName,folderName,varIndices,st,sx,sy,y,m,d,dt)
    nc_info = ncinfo(nc_fileName);
    varNames = cell(max(varIndices),1);
    for k = varIndices
        varNames{k} = nc_info.Variables(k).Name;
        out.(varNames{k}) = ncread(nc_fileName,varNames{k});
    end
    for t = 1:st
        for k = varIndices
            switch k
                case 1
                    out.(['data',num2str(t)]).(varNames{k}) = out.(varNames{k})(t);
                otherwise
                    out.(['data',num2str(t)]).(varNames{k}) = out.(varNames{k})(1:sx,1:sy,t);
            end
        end
        structSave = out.(['data',num2str(t)]);
        fileName = [num2str(y), sprintf('%02d',m), sprintf('%02d',d), sprintf('%02d',(t-1)*dt)];
        save(fullfile(folderName,fileName),'-struct','structSave');
        disp(['Data saved for ', fileName]);
    end
end

function [X, Y, uniX, uniY] = xyfromlatlon(lat,lon)
    % Actual grid
    [LON LAT] = meshgrid(lon,lat);
    latlon1 = [LAT(:) LON(:)];
    latlon2x = [LAT(:), ones(size(latlon1,1),1)*min(lon)];
    latlon2y = [ones(size(latlon1,1),1)*min(lat), LON(:)];

    DistX = lldistkm_edited(latlon1,latlon2x)*1000;
    DistY = lldistkm_edited(latlon1,latlon2y)*1000;
    X = reshape(DistX,size(LON));
    Y = reshape(DistY,size(LAT));
    
    % Cartesian grid - has an error by definition
    uniX = repmat(X(1,:),size(X,1),1);
    uniY = Y;
end

% dont do the above; convert velocities to \deg/second