function data = load_STM_dat(file_path)
    % Load 2D data from a Scienta text file and return an OxA_CUT object

    % Open the file and read all lines
    fileID = fopen(file_path);
    raw_lines = textscan(fileID, '%s', 'Delimiter', '\n');
    fclose(fileID);
    lines = raw_lines{1};

    % read position
    x_line = lines{find(startsWith(lines, 'X (m)'), 1)};
    x = sscanf(x_line(6:end), '%f')';

    y_line = lines{find(startsWith(lines, 'Y (m)'), 1)};
    y = sscanf(y_line(6:end), '%f')';

    z_line = lines{find(startsWith(lines, 'Z (m)'), 1)};
    z = sscanf(z_line(6:end), '%f')';


    % Find the index of the line starting with '[Data', indicating the start of the data section
    data_start_idx = find(startsWith(lines, '[DATA]'), 1);

    % Initialize the 'value' matrix and read the data from the lines
    value = [];
    for i = 1:(length(lines) - data_start_idx - 1)
        data_line = lines{data_start_idx + i + 1};
        value(i, :) = sscanf(data_line, '%f');
    end

    Bias_V = value(:,1);
    Current_A = value(:,3);
    Current_bwd_A = value(:,21);
    Current_both_A = (Current_A+Current_bwd_A)/2;

    data.info.x = x*1E9;
    data.info.y = y*1E9;
    data.info.z = z*1E9;
    data.Bias_V = Bias_V;
    data.Current_A = Current_A;
    data.Current_bwd_A = Current_bwd_A;
    data.Current_both_A = Current_both_A;

end



