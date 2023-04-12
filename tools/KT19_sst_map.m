function [fileout_mat] = KT19_sst_map(dirinspan,dirinkt19,dirout,varargin)

    %%%%
    %
    % [fileout_mat] = KT19_sst_map(dirinspan,dirinkt19,dirout,varargin)
    % 
    % Function for creating sea surface temperature map from a MASS Flight
    % and saving data into a .mat/kml files. 
    %
    % Written by Nick Statom with adaptions by Luke Colosi. 
    %
    %   Parameters
    %   ----------
    %   dirinspan : Path to SPAN (Novatel IMU) data for the trajectory of
    %               the aircraft during flight. File formats are ASCII
    %               files. 
    %
    %   dirinkt19 : Path to KT19 (point measurement pyrometer) data for the
    %               sea surface temperature observations. File formats are
    %               text files. 
    %
    %   dirout : Path to directory where processed KT19 data (.mat file) 
    %            and figures will be saved to.  
    %
    %   Varargin (Optional parameters) :
    %       time_offset : Specifies the time offset for the KT19 time stamp
    %                     to account for the time drift in the internal 
    %                     computer clock if it is not synced. Keyword 
    %                     argument is a scalar value in units of seconds. 
    %                     Default value is 0 seconds. 
    %       region : Latitude and longitude ranges saved in a matrix: 
    %                   region = [lon_min, lon_max; 
    %                             lat_min, lat_max]
    %                Default value is an empty matrix (i.e. plot entire
    %                flight path). 
    %       dirinkml : Path to kml file for flight lines including file
    %                  name. Keyword argument is a string. Default value is
    %                  an empty matrix so that nothing gets plotted. 
    %       dirinp : Path to a ploygon kml file. Default value is
    %                an empty matrix so that nothing gets plotted.
    %       dt : Time period to average temperature data. Units: seconds. 
    %       c_range : SST range for colorbar saved in a matrix: 
    %                   c_range = [cmin, cmax]
    %                 Default value is [13, 19]. Units: celcius.
    %       ds : Along track increments for kml file. Integer value that
    %            specifies the increment of SST data points saved in the
    %            kml. Default value is 10. 
    %       alt_max : Maximum altitude threshold value. Removes all data
    %                 in the SST data from altitudes higher than the
    %                 threshold value. Units: meters. Default is an empty
    %                 array. (i.e. all altitudes are plotted) 
    %       dT : Temperature offset for rough atmospheric correction. 
    %            Default value is zero. Units: celcius.   
    % 
    %   Returns
    %   -------
    %   Plots Aircraft trajectory and altitude, and KT19 sea
    %   surface temperacture. Saves data within area specified by region in
    %   a mat and kml file. 
    %
    %   Notes
    %   -----
    %   (1) To-do list: 
    %        (a) Add an objective mapping capabilities. 
    %        (b) Figure out how to isolate clouds contaminated data from the sea
    %            surface data.
    %    
    %   !!!!!!!!!!!!!!!!!!!!!!!!!!! Important  !!!!!!!!!!!!!!!!!!!!!!!!!!!
    %
    %   (2) The coordinates of the experimental site are HARD CODED into 
    %       the function!!!!! Make sure to change variable region_exp to
    %       the coordinates of your experiment!!!!!
    % 
    %   (3) The required toolboxes for this function include
    %         (a) m_map :  Geospatial Mapping toolbox found online at 
    %                      https://www.eoas.ubc.ca/~rich/map.html#download
    %         (b) gescatter : Function for saving KML files that render as
    %                         scatter plots in google earth. Luc sent me
    %                         this function, but it should be available 
    %                         through MatLab's plotting toolbox. 
    %         (c) kml2str : Found online at https://www.mathworks.com/matlabcentral/fileexchange/35642-kml2struct
    % 
    %   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %   
    %%%%

    clc, close all

    % Set Varargins
    p = inputParser();
    p.CaseSensitive = true;
    p.addOptional('time_offset',0);
    p.addOptional('region',[]);
    p.addOptional('dirinkml',[]);
    p.addOptional('dirinp',[]);
    p.addOptional('dt',0.5);
    p.addOptional('c_range',[13, 19]);
    p.addOptional('ds',10);
    p.addOptional('alt_max',[]);
    p.addOptional('dT',0);
    p.parse(varargin{:});
    time_offset = p.Results.time_offset;
    region = p.Results.region;
    dirinkml = p.Results.dirinkml;
    dirinp = p.Results.dirinp;
    dt = p.Results.dt;
    c_range = p.Results.c_range;
    ds = p.Results.ds;
    alt_max = p.Results.alt_max;
    dT = p.Results.dT;

    % Set default interpreter to latex for all text in figures
    set(groot, 'DefaultTextInterpreter', 'latex')
    set(groot, 'DefaultLegendInterpreter', 'latex')
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 

    % Set date constants
    tcst = datenum(1980,1,6,0,0,0);                                         % gps zero time constant
    leap = 18;                                                              % current leap seconds since time zero (will change based on trajectory date, see wikipedia)
    
    %% Extract BESTPOS lines from ascii converted SPAN trajectory file to get relevant time/lat/lon
    
    % Set filename, counter, and waitbar
    input = dir([dirinspan '*.ASC']);
    step = 0;
    
    % Generate Waitbar
    f = waitbar(0,'Please wait...','Name','Compiling SPAN Trajectory'); 
    pause(1)

    % Set time stamp date string for files and plot titles
    KT19.timestamp_plot = datestr(input.datenum,'mmm. dd, yyyy'); 
    KT19.timestamp_file = datestr(input.datenum,'yyyymmdd');

    % Set path for output mat and kml files for SST with filename
    fileout_mat = [dirout 'KT19_'  KT19.timestamp_file '.mat'];
    fileout_kml = [dirout 'KT19_'  KT19.timestamp_file '.kml'];
    
    % Open file for reading
    waitbar(.33,f,'Reading SPAN ASCII file.');
    pause(1)
    fid = fopen([dirinspan input.name]);
               
    % Update wait bar
    waitbar(.67,f,'Extracting BESTPOS lines from ascii.');

    % Loop through while end of file indicator (feof) is false 
    while ~feof(fid)
    
        % Read line excluding newline character
        tline = fgetl(fid);
    
        % Extract data if BESTPOSA string is in line
        if tline(2:9)=='BESTPOSA'
    
            % Set counter
            step=step+1;        
    
            % Find commas in sting
            ind = strfind(tline,',');
    
            % Obtain the week
            week = str2num(tline(ind(5)+1:ind(6)-1));
    
            % Save the time, latitude, and longitude from the SPAN
            KT19.time_span(step) = str2num(tline(ind(6)+1:ind(7)-1))/(60*60*24) + tcst - leap/(60*60*24) + week*7;
            KT19.lat(step) = str2num(tline(ind(11)+1:ind(12)-1));
            KT19.lon(step) = str2num(tline(ind(12)+1:ind(13)-1));
            KT19.alt(step) = str2num(tline(ind(13)+1:ind(14)-1));
    
        end
    end

    % Close SPAN file and waiting bar
    fclose(fid);
    waitbar(1,f,'Finishing');
    pause(1)
    close(f)

    %% Open KT19 files and extract time and temperature information
    
    % Set filenames and counter
    input=dir([dirinkt19 '*.txt']);
    step=0;
    
    % Generate Waitbar
    f = waitbar(0,'Please wait...','Name','Compiling KT19 SST Data');
    
    % Loop through files 
    for i=1:length(input)
    
        % Open ith file for reading
        fid=fopen([dirinkt19 input(i).name]);
    
        % Loop through while end of file indicator (feof) is false 
        while ~feof(fid)
    
            % Read line excluding newline character
            tline = fgetl(fid);
    
            % Set counter
            step=step+1;
    
            % Find commas in sting
            ind = strfind(tline,',');
    
            % Extract the time and temperature from the KT19
            time = tline(1:ind(1)-1);
            KT19.time_kt19(step) = datenum(time,'yyyymmdd_HHMMSSFFF');
            KT19.temp(step) = str2num(tline(ind(1)+1:end-2));
    
        end
    
        % Close ith KT19 file
        fclose(fid);
    
        % Update waitbar
        waitbar(i/length(input),f,...
            [num2str(round((i/length(input))*100)) '$\%$ complete...']);
    
    end
    
    % Close Waitbar
    close(f)

    %% Process data 

    % Shift time where KT19 timestamp is incorrect
    if time_offset
        KT19.time_kt19 = KT19.time_kt19 - time_offset/(24*60*60);
    end
    
    % Generate Waitbar
    f = waitbar(0,'Please wait...','Name','Computing dt second mean and std');
    
    % Find all kt19 data +/- dt second from bestpos span time and store in a standard deviation and mean vector
    for i=1:length(KT19.time_span)
    
        % Find indices for dt second average
        ind=find(KT19.time_kt19 < KT19.time_span(i) + (dt/2)/(60*60*24) & KT19.time_kt19 >= KT19.time_span(i) - (dt/2)/(60*60*24));
    
        % Compute 1 minute mean and std
        KT19.temp_std(i) = nanstd(KT19.temp(ind));
        KT19.temp_m(i) = nanmean(KT19.temp(ind));
    
        % Update waitbar
        waitbar(i/length(KT19.time_span),f,...
            [num2str(round((i/length(KT19.time_span))*100)) '$\%$ complete...']);
    end
    
    % Close Waitbar
    close(f)

    if ~isempty(region)
    
        % Find the region where data is collected in the Experimental site
        ind_lat = KT19.lat > region(2,1) & KT19.lat < region(2,2);
        ind_lon = KT19.lon > region(1,1) & KT19.lon < region(1,2);
        ind = logical(ind_lat.*ind_lon); 

        % Extract data from the region 
        KT19_r.time = KT19.time_span(ind);
        KT19_r.lat = KT19.lat(ind);
        KT19_r.lon = KT19.lon(ind);
        KT19_r.alt = KT19.alt(ind);
        KT19_r.temp_m = KT19.temp_m(ind);
        KT19_r.temp_std = KT19.temp_std(ind);

    else

        % Set SST to KT19  
        KT19_r.time = KT19.time_span;
        KT19_r.lat = KT19.lat;
        KT19_r.lon = KT19.lon;
        KT19_r.alt = KT19.alt;
        KT19_r.temp_m = KT19.temp_m;
        KT19_r.temp_std = KT19.temp_std;
    
    end

    if ~isempty(alt_max)

        % Find the altitudes less than or equal to alt_max where data is collected 
        ind_alt = SST.alt <= alt_max;

        % Extract data from specofied altitude 
        SST.time = KT19_r.time(ind_alt);
        SST.lat = KT19_r.lat(ind_alt);
        SST.lon = KT19_r.lon(ind_alt);
        SST.alt = KT19_r.alt(ind_alt);
        SST.temp = KT19_r.temp_m(ind_alt) + dT;                               % Add atmospheric correction
        SST.temp_std = KT19_r.temp_std(ind_alt);

    else 

        % Set SST to KT19  
        SST.time = KT19_r.time;
        SST.lat = KT19_r.lat;
        SST.lon = KT19_r.lon;
        SST.alt = KT19_r.alt;
        SST.temp = KT19_r.temp_m + dT;                                        % Add atmospheric correction
        SST.temp_std = KT19_r.temp_std;

    end

    
    % Save data as mat 
    save(fileout_mat,'SST')

    % Save data as a kml file
    gescatter(fileout_kml,SST.lon(1:ds:end),SST.lat(1:ds:end),SST.temp(1:ds:end),'scale',0.8)

    %% Plot Trajectories, Time series, and SST maps

    clc, close all 

    % Process kml lines file if specified in varargin
    if ~isempty(dirinkml)

        % Read in kml file for flight lines and set lat and lon arrays
        kmlstruct = kml2struct(dirinkml);
        lat = []; lon = [];

        % Loop through lines
        for ilines = 1:length(kmlstruct)

            % Extract filght lines
            lat = [lat fliplr(kmlstruct(ilines).Lat')  NaN];
            lon = [lon fliplr(kmlstruct(ilines).Lon')  NaN];

        end

        % Save longitude and latitude end points in a coordinates matrix
        coordinates = [lat; lon];

    end

    % Process kml ploygon file if specified in varargin
    if ~isempty(dirinp)

        % Read in kml file for flight lines and set lat and lon arrays
        kmlstructp = kml2struct(dirinp);
        lat_p = []; lon_p = []; 

        % Loop through lines
        for ivertex = 1:length(kmlstructp.Lat)

            if ivertex ~= length(kmlstructp.Lat)

                % Extract sides of polygon and save as lines
                lat_p = [lat_p kmlstructp.Lat(ivertex) kmlstructp.Lat(ivertex+1)];
                lon_p = [lon_p kmlstructp.Lon(ivertex) kmlstructp.Lon(ivertex+1)];

            elseif ivertex == length(kmlstructp.Lat)

                % Extract sides of polygon and save as lines
                lat_p = [lat_p kmlstructp.Lat(ivertex) kmlstructp.Lat(1)];
                lon_p = [lon_p kmlstructp.Lon(ivertex) kmlstructp.Lon(1)];

            end

        end

        % Save longitude and latitude end points in a coordinates matrix
        coordinates_p = [lat_p; lon_p];

    end 

    % Set plotting parameters
    t_ticks = datenum((floor(SST.time(1)*24)*(1/24)):hours(1):(ceil(SST.time(end)*24)*(1/24))); % Time tick marks
    region_exp = [35.5 38.5; -127, -121]; %!!!!!!!!! Hard Code !!!!!!!!                         % Longitude and latitude extent for m_map 
    if ~isempty(dirinkml)
        region_flight = [round(min([SST.lat, coordinates(1,:)],[],'all')-0.2,1)...              % Flight over flight lines extent
                         round(max([SST.lat, coordinates(1,:)],[],'all')+0.2,1);...
                         round(min([SST.lon, coordinates(2,:)],[],'all')-0.2,1)...
                         round(max([SST.lon, coordinates(2,:)],[],'all')+0.2,1)];
    else 
        region_flight = [round(min(SST.lat,[],'all')-0.2,1)...                                  % Flight over flight lines extent
                         round(max(SST.lat,[],'all')+0.2,1);...
                         round(min(SST.lon,[],'all')-0.2,1)...
                         round(max(SST.lon,[],'all')+0.2,1)];
    end
    land_c = [0.5, 0.5, 0.5];                                                                   % Color of land for m_map
    fontsize = 15;

    %---------------- SPAN Trajectory ----------------% 
    % Create figure
    figure('Name', 'SPAN Trajectory', 'units','normalized','outerposition',[0 0 0.5 0.8])
    
    % Generate map projection, coastlines, and grid
    m_proj('Mercator','long',region_exp(2,:),'lat',region_exp(1,:));
    m_gshhs_f('patch',land_c ,'edgecolor','k');
    m_grid('linest','none', 'tickstyle', 'dd', 'tickdir','out','box','fancy','fontsize',fontsize);

    % Plot trajectory from SPAN onto projection 
    if ~isempty(dirinkml) && ~isempty(dirinp)
        hold on 
            m_plot(coordinates(2,:), coordinates(1,:),'k-','Linewidth',0.1);
            m_plot(coordinates_p(2,:), coordinates_p(1,:),'k-','Linewidth',1);
            m_scatter(KT19.lon(KT19.lon < 0),KT19.lat(KT19.lon < 0),5,KT19.alt(KT19.lon < 0))
    elseif ~isempty(dirinkml) && isempty(dirinp)
        hold on 
            m_plot(coordinates(2,:), coordinates(1,:),'k-','Linewidth',0.1);
            m_scatter(KT19.lon(KT19.lon < 0),KT19.lat(KT19.lon < 0),5,KT19.alt(KT19.lon < 0))
    elseif isempty(dirinkml) && ~isempty(dirinp)
        hold on 
            m_plot(coordinates_p(2,:), coordinates_p(1,:),'k-','Linewidth',1);
            m_scatter(KT19.lon(KT19.lon < 0),KT19.lat(KT19.lon < 0),5,KT19.alt(KT19.lon < 0))
    elseif isempty(dirinkml) && isempty(dirinp)
        hold on
            m_scatter(KT19.lon(KT19.lon < 0),KT19.lat(KT19.lon < 0),5,KT19.alt(KT19.lon < 0))
    end
    
    % Set figure attributes
    title('Aircraft Trajectory and Altitude')
    set(gca,'FontSize',fontsize)

    % Set colorbar
    cb = colorbar;
    colormap(flipud(cbrewer2('Blues')))
    set(gca,'ColorScale','log')
    caxis([10^3, 3*10^3]);
    cb.Label.Interpreter = 'Latex';
    cb.Label.String = 'Altitude (km)';
    cb.TickDirection = 'out';
    cb.TickLength = 0.015;
    cb.FontSize = fontsize;
    cb.TickLabelInterpreter = 'latex';

    % Save Figure
    saveas(gcf, [dirout 'SPAN_'  KT19.timestamp_file '_traj_alt.png'])
    
    %---------------- Temperature Histogram ----------------% 
    
    % Find time indices where aircraft is in experiment site specified by region
    ind_time = KT19.time_kt19 >= min(SST.time) & KT19.time_kt19 <= max(SST.time);

    % Create figure
    figure('Name', 'KT19 SST histogram Raw', 'units','normalized','outerposition',[0 0 0.5 0.8]);
    
    % Plot histogram
    h = histogram(KT19.temp(ind_time),min(KT19.temp(ind_time)):0.08:max(KT19.temp(ind_time))); 
    
    % Set figure Attributes
    h.Normalization = 'pdf';
    % h.EdgeColor = 'none';
    xlabel('Temperature ($^\circ$C)')
    ylabel('Probability Density')

    % Save Figure
    saveas(gcf, [dirout 'KT19_'  KT19.timestamp_file '_sst_hist.png'])
    
    %---------------- KT19 SST Map dt second Average ----------------% 
    % Create figure
    figure('Name', ['KT19 SST Map ' num2str(dt) ' second Average'], 'units','normalized','outerposition',[0 0 1 1])
    
    %------------- Subplot 1: Time -------------%
    ax1 = subplot(1,2,1);
    
    % Generate map projection, coastlines, and grid
    m_proj('Mercator','long',region_flight(2,:),'lat',region_flight(1,:));
    m_gshhs_f('patch',land_c ,'edgecolor','k');
    m_grid('linest',':', 'tickstyle', 'dd', 'tickdir','out','box','fancy','fontsize',fontsize);
    
    % Plot  
    if ~isempty(dirinkml)
        hold on 
            m_plot(coordinates(2,:), coordinates(1,:),'k-','Linewidth',0.1);
            m_scatter(SST.lon,SST.lat,5,SST.time);
    else
        hold on
            m_scatter(SST.lon,SST.lat,5,SST.time);
    end

    % Set Figure Attributes
    title(['UTC Time on ' KT19.timestamp_plot ]);
    set(gca,'FontSize',fontsize)

    % Set colorbar attributes
    cb = colorbar;
    colormap(ax1, flipud(cbrewer2('Spectral')))
    cb.Ticks = t_ticks;
    caxis([t_ticks(1) t_ticks(end)])
    cb.TickLabels = datestr(cb.Ticks, 'HH');
    cb.Label.String = ['UTC Time on ' KT19.timestamp_plot];
    cb.Label.Interpreter = 'Latex';
    cb.Label.String = 'hours';
    cb.TickDirection = 'out';
    cb.TickLength = 0.015;
    cb.FontSize = fontsize;
    cb.TickLabelInterpreter = 'latex';
    
    %------------- Subplot 2: Mean Temp -------------%
    ax2 = subplot(1,2,2);

    % Generate map projection, coastlines, and grid
    m_proj('Mercator','long',region_flight(2,:),'lat',region_flight(1,:));
    m_gshhs_f('patch',land_c ,'edgecolor','k');
    m_grid('linest',':', 'tickstyle', 'dd', 'tickdir','out','box','fancy','fontsize',fontsize);
    
    % Plot  
    if ~isempty(dirinkml)
        hold on 
            m_plot(coordinates(2,:), coordinates(1,:),'k-','Linewidth',0.1);
            m_scatter(SST.lon,SST.lat,5,SST.temp);
    else
        hold on
            m_scatter(SST.lon,SST.lat,5,SST.temp);
    end

    
    % Set Figure Attributes
    title(['KT19 SST ',  num2str(dt), ' second Mean']);
    set(gca,'FontSize',fontsize)

    % Set colorbar attributes
    cb = colorbar;
    colormap(ax2, flipud(cbrewer2('RdYlBu')))
    caxis([c_range(1) c_range(2)])
    cb.Label.Interpreter = 'Latex';
    cb.Label.String = '$^\circ$C';
    cb.TickDirection = 'out';
    cb.TickLength = 0.015;
    cb.FontSize = fontsize;
    cb.TickLabelInterpreter = 'latex';

    % Save Figure
    saveas(gcf, [dirout 'KT19_'  KT19.timestamp_file '_sst_map.png'])
    
