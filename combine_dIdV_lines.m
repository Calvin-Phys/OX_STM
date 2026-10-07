function CUT = combine_dIdV_lines(var_list)
    position = [];
    value = [];
    
    for i = 1:length(var_list)
        file_path = ['C:\Users\pengc\Dropbox\@Oxford\Group Meeting\1_Germanium\Ge STM data\second exp\20230727Ge_2\Bias-Spectroscopy00' sprintf('%03.0f', var_list(i)) '.dat'];
        data = load_STM_dat(file_path);
        position(i,1) = data.info.x;
        position(i,2) = data.info.y;
        position(i,3) = data.info.z;
        value(i,:) = data.Current_both_A;
    end
    
    % figure
    % scatter3(position(:,1),position(:,2),position(:,3));
    
    CUT = OxA_CUT(1:length(var_list),data.Bias_V,value);
    CUT.x_name = 'Points';
    CUT.x_unit = 'Idx';
    CUT.y_name = 'Bias Voltage';
    CUT.y_unit = 'Volt';
    CUT.info.position = position;

end