classdef OX_STM < OxArpes_2D_Data
    % OxA_CUT: A class for handling 2D ARPES data - CUT.
    % This class inherits from the OxArpes_2D_Data class and adds methods
    % for processing ARPES data, such as converting to k-space, smoothing,
    % and calculating second derivatives and curvature.
    
    properties

    end
    
    methods
        function obj = OX_STM(inputData, iCh, isSimulated)
            % OxA_CUT: Class constructor.
            % Initializes the properties of the class with default values.
            if nargin < 3
                isSimulated = false; % Default to real data if not specified
            end

            if isSimulated
                % If the input data is simulated, create an empty object
                obj.x = [];
                obj.y = [];
                obj.value = [];
                obj.x_name = 'Position X';
                obj.x_unit = 'nm';
                obj.y_name = 'Position Y';
                obj.y_unit = 'nm';
                obj.name = 'Simulated STM Data';
            else
                % If the input data is real, process it as before
                obj.x = 1E9.*( linspace(-inputData.header.scan_range(1)/2, inputData.header.scan_range(1)/2, inputData.header.scan_pixels(1)) + inputData.header.scan_offset(1));
                obj.y = 1E9.*( linspace(-inputData.header.scan_range(2)/2, inputData.header.scan_range(2)/2, inputData.header.scan_pixels(2)) + inputData.header.scan_offset(2));
                if iCh == 0
                    obj.value = flip(inputData.channels(1).data', 2) + flip(inputData.channels(2).data', 2);
                else
                    obj.value = flip(inputData.channels(iCh).data', 2);
                end
                obj.x_name = 'Position X';
                obj.x_unit = 'nm';
                obj.y_name = 'Position Y';
                obj.y_unit = 'nm';

                name = split(inputData.header.scan_file, '\');
                name = name{7};
                range1 = num2str(round(inputData.header.scan_range(1) * 10^9, 1));
                range2 = num2str(round(inputData.header.scan_range(2) * 10^9, 1));
                bias = num2str(inputData.header.bias);
                current = num2str(round(10^12 * str2num(inputData.header.current_current__a_)));
                obj.name = [name ' ' range1 'n*' range2 'n ' bias 'V ' current 'pA'];
            end
        end

        function CUT_out = remove_STM_background(obj)
            [X, Y] = meshgrid(obj.x, obj.y);
            sf = fit([X(:), Y(:)], obj.value(:), 'poly55');
            FIT_SURF = sf(X, Y);
               
            CUT_out = obj;
            CUT_out.value = obj.value - FIT_SURF;
        end

        function generate_fig(obj)

            set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'});

            h1 = figure('Name',obj.name,'Renderer', 'painters');
%             h1.Units = 'centimeters';
%             h1.Position = [3 3 15 15];
            ha1 = axes('parent',h1);
            % ha1.Units = 'centimeters';
            % ha1.Position = [2 2 3.5 3.2];

            set(ha1,'linewidth',1.5);
            fontname(ha1,"Arial");
            
            % imagesc meshc
            % hold on
            imagesc(ha1,obj.x,obj.y,obj.value');
%             xlabel(ha1,[obj.x_name ' (' obj.x_unit ')']);
%             ylabel(ha1,[obj.y_name ' (' obj.y_unit ')']);
%             title(ha1,append(obj.name, ': ', num2str(round(obj.info.photon_energy,1)), 'eV ', obj.info.polarization),'interpreter', 'none');
            set(ha1,'YDir','normal');

            load('oxa_colourmap.mat');
            colormap(ha1,oxa_blue);

%             clim(ha1,[6 *min(obj.value,[],"all") 1*max(obj.value,[],"all")]);

%             set(ha1,'TickDir','out');
            box(ha1,'off');
%             set(ha1,'XMinorTick','on','YMinorTick','on');
%             ha1.XAxis.MinorTickValues = interp1(1:length(ha1.XTick),ha1.XTick,0.5:0.5:(length(ha1.XTick)+0.5),'linear','extrap');
%             ha1.YAxis.MinorTickValues = interp1(1:length(ha1.YTick),ha1.YTick,0.5:0.5:(length(ha1.YTick)+0.5),'linear','extrap');
            % axis tight

            % This syntax just sets the axis limits to their current value
            xlim(ha1, xlim(ha1));
            ylim(ha1, ylim(ha1));
            
            % Set right and upper axis lines to same color as axes
%             xline(max(xlim(ha1)), 'k-', 'Color', ha1.XAxis.Color);
%             yline(max(ylim(ha1)), 'k-', 'Color', ha1.YAxis.Color);

            if strcmp(obj.x_unit,obj.y_unit)
                pbaspect([1 1 1]);
            end
            set(ha1,'fontsize',10);
            % set(findall(gcf,'-property','FontSize'),'FontSize',6);
        end

        function CUT_out = stm_cut_fft(obj)
            CUT_out = obj;
            CUT_out.value = abs(fftshift(fft2(obj.value)));
            CUT_out.x = CUT_out.x - mean(CUT_out.x);
            CUT_out.x = CUT_out.x / (CUT_out.x(2) - CUT_out.x(1)) / (max(CUT_out.x) - min(CUT_out.x));
        
            CUT_out.y = CUT_out.y - mean(CUT_out.y);
            CUT_out.y = CUT_out.y / (CUT_out.y(2) - CUT_out.y(1)) / (max(CUT_out.y) - min(CUT_out.y));
        
            CUT_out.x_name = 'KX';
            CUT_out.x_unit = 'nm^{-1}';
            CUT_out.y_name = 'KY';
            CUT_out.y_unit = 'nm^{-1}';
        end

        function CUT_out = stm_cut_2ft(obj, kx, ky)
            CUT_out = OxA_CUT();
        
            CUT_out.x = kx;
            CUT_out.y = ky;
            CUT_out.value = zeros(length(kx), length(ky));
        
            [Y, X] = meshgrid(obj.y, obj.x);
        
            mat = zeros(length(kx), length(ky));
            mm2 = obj.value(:);
            kx_l = 1:length(kx);
            ky_l = 1:length(ky);
            for r = kx_l
                for s = ky_l
                    mask = exp(-2 * pi * 1i * (kx(r) * X + ky(s) * Y));
                    mat(r, s) = dot(mask(:), mm2);
                end
            end
            CUT_out.value = abs(mat);
         
            CUT_out.x_name = 'KX';
            CUT_out.x_unit = 'nm^{-1}';
            CUT_out.y_name = 'KY';
            CUT_out.y_unit = 'nm^{-1}';
        end

        % Method to set the size and resolution of the simulated data
        function obj = setSimulatedAxes(obj, matrixSize, axisLength)
            obj.x = linspace(0, axisLength, matrixSize);
            obj.y = linspace(0, axisLength, matrixSize);
            obj.value = zeros(matrixSize, matrixSize);
        end

        % Method to generate atom peaks and apply Gaussian convolution
        function obj = generateSimulatedPeaks(obj, bondLength, gaussianSigma)
            % Ensure that the axes have been set
            if isempty(obj.x) || isempty(obj.y) || isempty(obj.value)
                error('Axes must be set before generating peaks.');
            end
            
            matrixSize = size(obj.value, 1);
            
            peakPeriod_x = floor(max(obj.x)/bondLength/3);
            peakPeriod_y = floor(max(obj.y)/bondLength/sqrt(3));
            

            % Create hexagonal lattice
            for i = 1:peakPeriod_x
                for j = 1:peakPeriod_y
%                     x_ind = round(bondLength *3 *i * length(obj.x) / max(obj.x));
%                     y_ind = round(bondLength *sqrt(3) *j * length(obj.y) / max(obj.y));
%                     obj.value(x_ind,y_ind) = 1;

                    x_ind = round(bondLength *(3 *i-1) * length(obj.x) / max(obj.x));
                    y_ind = round(bondLength *sqrt(3) *j * length(obj.y) / max(obj.y));
                    obj.value(x_ind,y_ind) = 1;

                    x_ind = round(bondLength *(3 *i-1.5) * length(obj.x) / max(obj.x));
                    y_ind = round(bondLength *sqrt(3) *(j-0.5) * length(obj.y) / max(obj.y));
                    obj.value(x_ind,y_ind) = 1;

%                     x_ind = round(bondLength *(3 *i-2.5) * length(obj.x) / max(obj.x));
%                     y_ind = round(bondLength *sqrt(3) *(j-0.5) * length(obj.y) / max(obj.y));
%                     obj.value(x_ind,y_ind) = 1;
                end
            end
            
            % Convolve with Gaussian
            gaussianFilter = fspecial('gaussian', [matrixSize, matrixSize], gaussianSigma/(obj.x(end)-obj.x(1))*length(obj.x));
            obj.value = conv2(obj.value, gaussianFilter, 'same');
            
            obj.name = 'Simulated STM Data';
        end
    
    
    end
end
