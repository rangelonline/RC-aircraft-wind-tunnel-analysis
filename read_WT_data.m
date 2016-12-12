%{
Function to read data files provided by the Oran W. Nicks Low Speed Wind 
Tunnel in College Station, Texas. The provided .txt files are formatted 
for Texas A&M University's aero 402 aerospace vehicle design II course. 

"Raw" data files contain both aerodynamic and inertial force and moment 
measurements. "Stat" data files contain only aerodynamic force and moment 
measurements. 

The Oran W. Nicks LSWT provides psuedo-coefficients of aerodynamic force or
moment divided by dynamic pressure. Useful aerodynamic coefficients are
calculated by normalizing these psuedo-coefficients by a user-specified
reference area or reference length. 

VARIABLES:
filename:               data file name. 
filetype:               specify whether filename is a "raw" or "stat" file.
Sref:                   vehicle reference area. 
cref:                   vehicle reference chord. 
data:                   struct to be passed back to main function or input
                        deck. Struct fields are headers. Each struct field
                        contains its data stored as 1D arrays. 

CALLED FUNCTIONS:
1)   write_struct_field_names(field_names)

NOTE: headers must be in order from left to right. 

By John Rangel 11/29/16
%}

function data = read_WT_data(filename, filetype, Sref, cref)
    
    %-----------------------------------------------------------------------
    %                      wind tunnel data file headers
    %-----------------------------------------------------------------------
    
    % environment, model orientation headers:
    headers = {'point'
               'dynamic_pressure' %Qact
               'velocity' %Uact
               'gage_total_pressure' %Pt
               'gage_static_pressure' %Ps
               'gage_barometric_pressure' %Baro
               'temperature' %Temp
               'phi' %commanded model roll angle
               'theta' %commanded model pitch angle
               'psi' %commanded model pitch angle
               'alpha' %actual model pitch angle
               'beta' %actual model yaw angle
               'gamma'}; %actual model roll angle
           
    % data without static tare reduction:
    raw_headers = {'axial_force' %AF
                    'side_force' %SF
                    'normal_force' %NF
                    'rolling_moment' %RM
                    'pitching_moment' %PM
                    'yawing_moment'}; %YM
                
    % data with static tare reduction:
    % NOTE: coefficients with a 'B' are in body frame. All others are in 
    %       wind frame. 
    stat_headers = {'CD' 
                    'CY'
                    'CL'
                    'Cl'
                    'Cm'
                    'Cn'
                    'CAB' 
                    'CYB'
                    'CNB'
                    'ClB'
                    'CmB'
                    'CnB'};
                
    raw_headers = [headers; raw_headers];
    stat_headers = [headers; stat_headers];
    
    %-----------------------------------------------------------------------
    %                         open data file
    %-----------------------------------------------------------------------
    
    fileID = fopen(filename);
    raw_data = fgetl(fileID);
    
    % ignore data file notes, headers:
    comment = raw_data(1);
    row_number = 1;
    while comment == '#'
        raw_data = fgetl(fileID);
        comment = raw_data(1);
        row_number = row_number+1;
    end %end while loop
    
    start_of_data = row_number-1;
    
    %-----------------------------------------------------------------------
    %                         read raw data
    %-----------------------------------------------------------------------
        
    % read numeric data:
    if strcmp(filetype, 'raw') == 1
        raw_data = dlmread(filename, '', start_of_data, 0);
        
        % check header, data match:
        number_of_file_headers = size(raw_data, 2);
        number_of_function_headers = size(raw_headers);
        if number_of_file_headers ~= number_of_function_headers
            disp('ERROR: file, function headers do not match.')
        end
        
        % create data struct:
        raw_data_struct = write_struct_field_names(raw_headers);
        
        % add data to struct:
        fields = fieldnames(raw_data_struct);
        for i = 1:number_of_file_headers
            raw_data_struct.(fields{i}) = raw_data(:, i);
        end %end i
        
        % trim first, last data points to account for system zeroing:
        trimmed_raw_data_struct = write_struct_field_names(raw_headers); 
        number_of_file_rows = size(raw_data, 1);
        
        for i = 1:number_of_file_headers
            count = 1;
            
            % store trimmed data in new struct:
            for j = 2:(number_of_file_rows-1)
                trimmed_raw_data_struct.(fields{i})(count) = raw_data_struct.(fields{i})(j);
                count = count+1;                
            end %end j
            
        end %end i
        
        % pass back trimmed raw data:
        data = trimmed_raw_data_struct;
    
    %-----------------------------------------------------------------------
    %                         read stat data
    %-----------------------------------------------------------------------
        
    elseif strcmp(filetype, 'stat') == 1
        stat_data = dlmread(filename, '', start_of_data, 0);
        
        % check header, data match:
        number_of_file_headers = size(stat_data, 2);
        number_of_function_headers = size(stat_headers);
        if number_of_file_headers ~= number_of_function_headers
            disp('ERROR: file, function headers do not match')
        end %end if
        
        % create data struct:
        stat_data_struct = write_struct_field_names(stat_headers);
        
        % add data to struct:
        fields = fieldnames(stat_data_struct);
        for i = 1:number_of_file_headers
            stat_data_struct.(fields{i}) = stat_data(:, i);
        end %end for loop
        
        % transform wind-frame force psuedo-coefficients:
        stat_data_struct.CD = stat_data_struct.CD/Sref;
        stat_data_struct.CL = stat_data_struct.CL/Sref;
        stat_data_struct.CY = stat_data_struct.CY/Sref;
        
        % transform wind-frame moment psuedo-coefficients:
        stat_data_struct.Cl = stat_data_struct.Cl/(Sref*cref);
        stat_data_struct.Cm = stat_data_struct.Cm/(Sref*cref);
        stat_data_struct.Cn = stat_data_struct.Cn/(Sref*cref);
        
        % transform body-frame force psuedo-coefficients:
        stat_data_struct.CAB = stat_data_struct.CAB/Sref;
        stat_data_struct.CYB = stat_data_struct.CYB/Sref;
        stat_data_struct.CNB = stat_data_struct.CNB/Sref;
        
        % transform body-frame moment psuedo-coefficients:
        stat_data_struct.ClB = stat_data_struct.ClB/(Sref*cref);
        stat_data_struct.CmB = stat_data_struct.CmB/(Sref*cref);
        stat_data_struct.CnB = stat_data_struct.CnB/(Sref*cref);
        
        % trim first, last data points to account for system zeroing:
        trimmed_stat_data_struct = write_struct_field_names(stat_headers); 
        number_of_file_rows = size(stat_data, 1);
        
        for i = 1:number_of_file_headers
            count = 1;
            
            % store trimmed data in new struct:
            for j = 2:(number_of_file_rows-1)
                trimmed_stat_data_struct.(fields{i})(count) = stat_data_struct.(fields{i})(j);
                count = count+1;                
            end %end j
            
        end %end i     
        
        % pass back trimmed stat data:
        data = trimmed_stat_data_struct;
        
    %-----------------------------------------------------------------------
    %                           error messages
    %-----------------------------------------------------------------------  
    
    else 
        disp('ERROR: data file type not supported. Choose (raw) or (stat)')
        data = 0.0;
                
    end % end if elseif else statement. 
    

    
