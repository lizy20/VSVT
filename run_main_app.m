classdef run_main_app < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        GridLayout                    matlab.ui.container.GridLayout
        Tree                          matlab.ui.container.CheckBoxTree
        AllstepsNode                  matlab.ui.container.TreeNode
        VariableDehomogenazationNode  matlab.ui.container.TreeNode
        VoronoiseedingNode            matlab.ui.container.TreeNode
        VoronoitessellationNode       matlab.ui.container.TreeNode
        TrussinfocomputationNode      matlab.ui.container.TreeNode
        BuildPolyshapeNode            matlab.ui.container.TreeNode
        BoundarySimplifyNode          matlab.ui.container.TreeNode
        GenerateSTLfileButton         matlab.ui.control.Button
        ComputegeometryButton         matlab.ui.control.Button
        BusyLamp                      matlab.ui.control.Lamp
        BusyLampLabel                 matlab.ui.control.Label
        ConsoleOutput                 matlab.ui.control.TextArea
        Panel                         matlab.ui.container.Panel
        IsotropicdesignLabel          matlab.ui.control.Label
        AnisotropicdesignLabel        matlab.ui.control.Label
        Switch                        matlab.ui.control.Switch
        VolumefractionEditField       matlab.ui.control.NumericEditField
        VolumeFractionEditFieldLabel  matlab.ui.control.Label
        R_minEditField                matlab.ui.control.NumericEditField
        R_minEditFieldLabel           matlab.ui.control.Label
        LocalporositySlider           matlab.ui.control.RangeSlider
        LocalRelativeDensityrhoLabel  matlab.ui.control.Label
        BoundaryConditionDropDown     matlab.ui.control.DropDown
        BoundaryconditionLabel        matlab.ui.control.Label
        ElementNumberHEditField       matlab.ui.control.NumericEditField
        ElementNumberHLabel           matlab.ui.control.Label
        ElementNumberLEditField       matlab.ui.control.NumericEditField
        ElementnumberLLabel           matlab.ui.control.Label
        DomainHeightEditField         matlab.ui.control.NumericEditField
        DomainHeightEditFieldLabel    matlab.ui.control.Label
        DomainLengthEditField         matlab.ui.control.NumericEditField
        DomainlengthLabel             matlab.ui.control.Label
        PauseResumeButton             matlab.ui.control.Button
        StartButton                   matlab.ui.control.Button
        LocalporosityMax              matlab.ui.control.NumericEditField
        LocalporosityMin              matlab.ui.control.NumericEditField
        Image                         matlab.ui.control.Image
        UIAxes5                       matlab.ui.control.UIAxes
        UIAxes3                       matlab.ui.control.UIAxes
        UIAxes4                       matlab.ui.control.UIAxes
        UIAxes2                       matlab.ui.control.UIAxes
        UIAxes                        matlab.ui.control.UIAxes
    end

    
      properties (Access = public)
        IsPaused =0; %     
    end
     properties (Access = public)
        Status =[1,1,1,1,1,1]; %     
     end

    properties (Access = private)
    BusyLampTimer           % 
    BusyLampBreathPhase = 0 % 
    end


    methods (Access = public)
        function status = getStepCheckVector(app)
             stepNodes = [
        app.VariableDehomogenazationNode;
        app.VoronoiseedingNode;
        app.VoronoitessellationNode;
        app.TrussinfocomputationNode;
        app.BuildPolyshapeNode 
        app.BoundarySimplifyNode  ];
      status = ismember(stepNodes, app.Tree.CheckedNodes);
        end
    end
    
    methods (Access = public)
      function startBusyLamp(app)
    if ~isempty(app.BusyLampTimer) && isvalid(app.BusyLampTimer)
        return;
    end
    app.BusyLampBreathPhase = 0;
    app.BusyLampTimer = timer( ...
        'ExecutionMode', 'fixedRate', ...
        'Period', 0.05, ...  % 
        'TimerFcn', @(~,~)updateLampColor(app) ...
    );
    start(app.BusyLampTimer);
     end
    end

    methods (Access = public)
    function updateLampColor(app)
    app.BusyLampBreathPhase = app.BusyLampBreathPhase + 0.05;
    alpha = (sin(app.BusyLampBreathPhase) + 1) / 2;  % 0 ~ 1
    darkRed = [0.55,0.98,0.82]; 
    brightRed = [0.39,0.83,0.07];
    currentColor = (1 - alpha) * darkRed + alpha * brightRed;
    app.BusyLamp.Color = currentColor;
    end
    end

   methods (Access = public)
      function stopBusyLamp(app)
            if ~isempty(app.BusyLampTimer) && isvalid(app.BusyLampTimer)
                stop(app.BusyLampTimer);
                delete(app.BusyLampTimer);
                app.BusyLampTimer = [];
            end
            app.BusyLamp.Color = [0.7, 0.7, 0.7];  %
        end
    end
        
      
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: StartButton
        function StartButtonPushed(app, event)
            
   %==================initialize===============
    cd(fileparts(mfilename('fullpath')));
    %=================================
    L = app.DomainLengthEditField.Value;
    W = app.DomainHeightEditField.Value;
    Nx = app.ElementNumberLEditField.Value;
    Ny = app.ElementNumberHEditField.Value;
    VolFrac = app.VolumefractionEditField.Value;
    rmin =app.R_minEditField.Value;
    app.IsPaused =0;
    app.PauseResumeButton.Text = 'Pause';
    LocalPorosityRange = app.LocalporositySlider.Value;
    LocalPorosityMin = LocalPorosityRange(1);
    LocalPorosityMax = LocalPorosityRange(2);
    Fkey = app.BoundaryConditionDropDown.Value;
    % Validate input parameters
    if L <= 0 || W <= 0 || Nx <= 0 || Ny <= 0 || VolFrac <= 0 || VolFrac > 1
        uialert(app.UIFigure, 'Invalid input parameters. Please check your values!', 'Error');
        return;
    end
    
    paths = genpath(pwd);
    addpath(paths);
    disp('All subfolders added to path.');
    myPrint(app, 'Computation started. Running optimization...');
    run_main(L, W,  Nx, Ny, VolFrac, LocalPorosityMin ,LocalPorosityMax,rmin,Fkey,app);
    % Display completion message
    uialert(app.UIFigure, ['Optimization completed successfully!' newline,...
        'Please go ahead, inverse-design geometry'],  'All Set', ...
    'Icon', 'success');
        end

        % Button pushed function: PauseResumeButton
        function PauseResumeButtonPushed(app, event)
               app.IsPaused = ~app.IsPaused;
          if app.IsPaused
           app.PauseResumeButton.Text = 'Resume';
           disp('Paused by user.');
           app.ConsoleOutput.Value = [ app.ConsoleOutput.Value; 'Paused by user.'];
          else
            app.PauseResumeButton.Text = 'Pause';
            disp('Resumed by user.');
             app.ConsoleOutput.Value = [app.ConsoleOutput.Value;'Resumed by user.'];
          end
                 scroll(app.ConsoleOutput, 'bottom')
                 drawnow;
        end

        % Value changed function: LocalporosityMin
        function LocalporosityMinValueChanged(app, event)
            minVal = app.LocalporosityMin.Value;
            maxVal = app.LocalporosityMax.Value;

            if minVal > maxVal
                minVal = maxVal;
                app.LocalporosityMin.Value = minVal;
            end

            minLimit = app.LocalporositySlider.Limits(1);
            if minVal < minLimit
                minVal = minLimit;
                app.LocalporosityMin.Value = minVal;
            end  
            app.LocalporositySlider.Value(1) = minVal;
        end

        % Value changed function: LocalporosityMax
        function LocalporosityMaxValueChanged(app, event)

        minVal = app.LocalporosityMin.Value;
        maxVal = app.LocalporosityMax.Value;
        if maxVal < minVal
            maxVal = minVal;
            app.LocalporosityMax.Value = maxVal;
        end
        maxLimit = app.LocalporositySlider.Limits(2);
        if maxVal > maxLimit
            maxVal = maxLimit;
            app.LocalporosityMax.Value = maxVal;
        end
        app.LocalporositySlider.Value(2) = maxVal;
        end

        % Value changed function: LocalporositySlider
        function LocalporositySliderValueChanged(app, event)
            val = app.LocalporositySlider.Value;
            app.LocalporosityMin.Value = val(1);
            app.LocalporosityMax.Value = val(2);
        end

        % Button pushed function: ComputegeometryButton
        function ComputegeometryButtonPushed(app, event)
            cd(fileparts(mfilename('fullpath')));
    %=================================
    L = app.DomainLengthEditField.Value;
    W = app.DomainHeightEditField.Value;
    Nx = app.ElementNumberLEditField.Value;
    Ny = app.ElementNumberHEditField.Value;
    VolFrac = app.VolumefractionEditField.Value;
    rmin =app.R_minEditField.Value;

    LocalPorosityRange = app.LocalporositySlider.Value;
    LocalPorosityMin = LocalPorosityRange(1);
    LocalPorosityMax = LocalPorosityRange(2);
    Fkey = app.BoundaryConditionDropDown.Value;
    % Validate input parameters
    if L <= 0 || W <= 0 || Nx <= 0 || Ny <= 0 || VolFrac <= 0 || VolFrac > 1
        uialert(app.UIFigure, 'Invalid input parameters. Please check your values!', 'Error');
        return;
    end
    paths = genpath(pwd);
    addpath(paths);
    disp('All subfolders added to path.');         
    myPrint(app, 'Building Geometry:...');
    run_geo_main(L, W,  Nx, Ny, VolFrac, LocalPorosityMin ,LocalPorosityMax,app);
    % Display completion message
    CStatus = app.getStepCheckVector();
    if  CStatus(6)==1 
    uialert(app.UIFigure, ['Geometry completed successfully!  ' newline ...
        'Please Using [Generate STL file] button for saving model.'], 'All Set', ...
    'Icon', 'success');
    else
       uialert(app.UIFigure, ['Selected Steps completed successfully!'], 'All Set', ...
    'Icon', 'success');
    end

        end

        % Button pushed function: GenerateSTLfileButton
        function GenerateSTLfileButtonPushed(app, event)
     cd(fileparts(mfilename('fullpath')));
    %=================================
    L = app.DomainLengthEditField.Value;
    W = app.DomainHeightEditField.Value;
    Nx = app.ElementNumberLEditField.Value;
    Ny = app.ElementNumberHEditField.Value;
    VolFrac = app.VolumefractionEditField.Value;
    rmin =app.R_minEditField.Value; 

    LocalPorosityRange = app.LocalporositySlider.Value;
    LocalPorosityMin = LocalPorosityRange(1);
    LocalPorosityMax = LocalPorosityRange(2);
    Fkey = app.BoundaryConditionDropDown.Value;
  
    % Validate input parameters
    if L <= 0 || W <= 0 || Nx <= 0 || Ny <= 0 || VolFrac <= 0 || VolFrac > 1
        uialert(app.UIFigure, 'Invalid input parameters. Please check your values!', 'Error');
        return;
    end

    paths = genpath(pwd);
    addpath(paths);
    disp('All subfolders added to path.');         
    myPrint(app, 'Building STL file:...');
    run_build_stl(L, W,  Nx, Ny, VolFrac, LocalPorosityMin ,LocalPorosityMax,app);
    % Display completion message
    uialert(app.UIFigure, ['STL built in the directory successfully!'],  'All Set', ...
    'Icon', 'success');
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1043 737];
            app.UIFigure.Name = 'MATLAB App';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout.ColumnSpacing = 2;
            app.GridLayout.RowSpacing = 2;

            % Create UIAxes
            app.UIAxes = uiaxes(app.GridLayout);
            title(app.UIAxes, 'Convergence history')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.FontSize = 12;
            app.UIAxes.Layout.Row = [1 12];
            app.UIAxes.Layout.Column = [13 35];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.GridLayout);
            title(app.UIAxes2, 'Anisotropy')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.ClippingStyle = 'rectangle';
            app.UIAxes2.TitleHorizontalAlignment = 'left';
            app.UIAxes2.FontSize = 12;
            app.UIAxes2.Layout.Row = [13 22];
            app.UIAxes2.Layout.Column = [13 19];

            % Create UIAxes4
            app.UIAxes4 = uiaxes(app.GridLayout);
            title(app.UIAxes4, '\alpha')
            zlabel(app.UIAxes4, 'Z')
            app.UIAxes4.ClippingStyle = 'rectangle';
            app.UIAxes4.FontSize = 12;
            app.UIAxes4.Layout.Row = [13 22];
            app.UIAxes4.Layout.Column = [29 35];

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.GridLayout);
            title(app.UIAxes3, '\rho')
            zlabel(app.UIAxes3, 'Z')
            app.UIAxes3.ClippingStyle = 'rectangle';
            app.UIAxes3.FontSize = 12;
            app.UIAxes3.Layout.Row = [13 22];
            app.UIAxes3.Layout.Column = [21 27];

            % Create UIAxes5
            app.UIAxes5 = uiaxes(app.GridLayout);
            title(app.UIAxes5, 'Overview material distribution')
            zlabel(app.UIAxes5, 'Z')
            app.UIAxes5.ClippingStyle = 'rectangle';
            app.UIAxes5.FontSize = 12;
            app.UIAxes5.Layout.Row = [24 32];
            app.UIAxes5.Layout.Column = [12 22];

            % Create Image
            app.Image = uiimage(app.GridLayout);
            app.Image.ScaleMethod = 'fill';
            app.Image.Layout.Row = [23 26];
            app.Image.Layout.Column = [24 35];
            app.Image.ImageSource = fullfile(pathToMLAPP, 'icon.png');

            % Create Panel
            app.Panel = uipanel(app.GridLayout);
            app.Panel.Layout.Row = [1 36];
            app.Panel.Layout.Column = [1 11];

            % Create LocalporosityMin
            app.LocalporosityMin = uieditfield(app.Panel, 'numeric');
            app.LocalporosityMin.Limits = [0 1];
            app.LocalporosityMin.ValueChangedFcn = createCallbackFcn(app, @LocalporosityMinValueChanged, true);
            app.LocalporosityMin.HorizontalAlignment = 'center';
            app.LocalporosityMin.FontName = 'Arial';
            app.LocalporosityMin.FontSize = 18;
            app.LocalporosityMin.Position = [24 303 78 32];
            app.LocalporosityMin.Value = 0.3;

            % Create LocalporosityMax
            app.LocalporosityMax = uieditfield(app.Panel, 'numeric');
            app.LocalporosityMax.Limits = [0 1];
            app.LocalporosityMax.ValueChangedFcn = createCallbackFcn(app, @LocalporosityMaxValueChanged, true);
            app.LocalporosityMax.HorizontalAlignment = 'center';
            app.LocalporosityMax.FontName = 'Arial';
            app.LocalporosityMax.FontSize = 18;
            app.LocalporosityMax.Position = [184 303 78 32];
            app.LocalporosityMax.Value = 0.8;

            % Create StartButton
            app.StartButton = uibutton(app.Panel, 'push');
            app.StartButton.ButtonPushedFcn = createCallbackFcn(app, @StartButtonPushed, true);
            app.StartButton.FontName = 'Arial';
            app.StartButton.FontSize = 18;
            app.StartButton.FontWeight = 'bold';
            app.StartButton.Position = [9 13 137 67];
            app.StartButton.Text = '1.Start';

            % Create PauseResumeButton
            app.PauseResumeButton = uibutton(app.Panel, 'push');
            app.PauseResumeButton.ButtonPushedFcn = createCallbackFcn(app, @PauseResumeButtonPushed, true);
            app.PauseResumeButton.FontName = 'Arial';
            app.PauseResumeButton.FontSize = 18;
            app.PauseResumeButton.Position = [173 21 105 50];
            app.PauseResumeButton.Text = 'Pause';

            % Create DomainlengthLabel
            app.DomainlengthLabel = uilabel(app.Panel);
            app.DomainlengthLabel.FontName = 'arial';
            app.DomainlengthLabel.FontSize = 16;
            app.DomainlengthLabel.Position = [9 664 117 32];
            app.DomainlengthLabel.Text = 'Domain Length';

            % Create DomainLengthEditField
            app.DomainLengthEditField = uieditfield(app.Panel, 'numeric');
            app.DomainLengthEditField.Limits = [0 Inf];
            app.DomainLengthEditField.HorizontalAlignment = 'center';
            app.DomainLengthEditField.FontName = 'Arial';
            app.DomainLengthEditField.FontSize = 18;
            app.DomainLengthEditField.Position = [157 664 132 32];
            app.DomainLengthEditField.Value = 100;

            % Create DomainHeightEditFieldLabel
            app.DomainHeightEditFieldLabel = uilabel(app.Panel);
            app.DomainHeightEditFieldLabel.FontName = 'arial';
            app.DomainHeightEditFieldLabel.FontSize = 16;
            app.DomainHeightEditFieldLabel.Position = [9 613 117 32];
            app.DomainHeightEditFieldLabel.Text = 'Domain Height';

            % Create DomainHeightEditField
            app.DomainHeightEditField = uieditfield(app.Panel, 'numeric');
            app.DomainHeightEditField.Limits = [0 Inf];
            app.DomainHeightEditField.HorizontalAlignment = 'center';
            app.DomainHeightEditField.FontName = 'Arial';
            app.DomainHeightEditField.FontSize = 18;
            app.DomainHeightEditField.Position = [158 613 132 32];
            app.DomainHeightEditField.Value = 60;

            % Create ElementnumberLLabel
            app.ElementnumberLLabel = uilabel(app.Panel);
            app.ElementnumberLLabel.HorizontalAlignment = 'center';
            app.ElementnumberLLabel.FontName = 'arial';
            app.ElementnumberLLabel.FontSize = 16;
            app.ElementnumberLLabel.Position = [4 551 148 52];
            app.ElementnumberLLabel.Text = 'Element Number  -L';

            % Create ElementNumberLEditField
            app.ElementNumberLEditField = uieditfield(app.Panel, 'numeric');
            app.ElementNumberLEditField.Limits = [1 Inf];
            app.ElementNumberLEditField.ValueDisplayFormat = '%.0f';
            app.ElementNumberLEditField.HorizontalAlignment = 'center';
            app.ElementNumberLEditField.FontName = 'Arial';
            app.ElementNumberLEditField.FontSize = 18;
            app.ElementNumberLEditField.Position = [157 561 132 32];
            app.ElementNumberLEditField.Value = 50;

            % Create ElementNumberHLabel
            app.ElementNumberHLabel = uilabel(app.Panel);
            app.ElementNumberHLabel.HorizontalAlignment = 'center';
            app.ElementNumberHLabel.FontName = 'arial';
            app.ElementNumberHLabel.FontSize = 16;
            app.ElementNumberHLabel.Position = [4 499 146 43];
            app.ElementNumberHLabel.Text = 'Element Number -H';

            % Create ElementNumberHEditField
            app.ElementNumberHEditField = uieditfield(app.Panel, 'numeric');
            app.ElementNumberHEditField.Limits = [1 Inf];
            app.ElementNumberHEditField.ValueDisplayFormat = '%.0f';
            app.ElementNumberHEditField.HorizontalAlignment = 'center';
            app.ElementNumberHEditField.FontName = 'Arial';
            app.ElementNumberHEditField.FontSize = 18;
            app.ElementNumberHEditField.Position = [157 509 132 32];
            app.ElementNumberHEditField.Value = 30;

            % Create BoundaryconditionLabel
            app.BoundaryconditionLabel = uilabel(app.Panel);
            app.BoundaryconditionLabel.FontName = 'arial';
            app.BoundaryconditionLabel.FontSize = 16;
            app.BoundaryconditionLabel.Position = [7 449 149 41];
            app.BoundaryconditionLabel.Text = 'Boundary Condition';

            % Create BoundaryConditionDropDown
            app.BoundaryConditionDropDown = uidropdown(app.Panel);
            app.BoundaryConditionDropDown.Items = {'Cantilever', 'MBB', 'Compression', 'Compression(Paper)', 'User-defined', 'User-Demo1', 'User-Demo2'};
            app.BoundaryConditionDropDown.ItemsData = [1 2 3 4 5 6 7];
            app.BoundaryConditionDropDown.FontName = 'Arial';
            app.BoundaryConditionDropDown.FontSize = 18;
            app.BoundaryConditionDropDown.Position = [158 453 132 32];
            app.BoundaryConditionDropDown.Value = 1;

            % Create LocalRelativeDensityrhoLabel
            app.LocalRelativeDensityrhoLabel = uilabel(app.Panel);
            app.LocalRelativeDensityrhoLabel.HorizontalAlignment = 'center';
            app.LocalRelativeDensityrhoLabel.FontName = 'Arial';
            app.LocalRelativeDensityrhoLabel.FontSize = 18;
            app.LocalRelativeDensityrhoLabel.Interpreter = 'latex';
            app.LocalRelativeDensityrhoLabel.Position = [30 389 226 32];
            app.LocalRelativeDensityrhoLabel.Text = 'Local Relative Density \rho';

            % Create LocalporositySlider
            app.LocalporositySlider = uislider(app.Panel, 'range');
            app.LocalporositySlider.Limits = [0 1];
            app.LocalporositySlider.ValueChangedFcn = createCallbackFcn(app, @LocalporositySliderValueChanged, true);
            app.LocalporositySlider.FontName = 'Arial';
            app.LocalporositySlider.FontSize = 18;
            app.LocalporositySlider.Position = [31 381 227 3];
            app.LocalporositySlider.Value = [0.3 0.8];

            % Create R_minEditFieldLabel
            app.R_minEditFieldLabel = uilabel(app.Panel);
            app.R_minEditFieldLabel.HorizontalAlignment = 'center';
            app.R_minEditFieldLabel.FontName = 'Arial';
            app.R_minEditFieldLabel.FontSize = 16;
            app.R_minEditFieldLabel.Interpreter = 'tex';
            app.R_minEditFieldLabel.Position = [24 252 105 32];
            app.R_minEditFieldLabel.Text = 'R_{min}';

            % Create R_minEditField
            app.R_minEditField = uieditfield(app.Panel, 'numeric');
            app.R_minEditField.HorizontalAlignment = 'center';
            app.R_minEditField.FontName = 'Arial';
            app.R_minEditField.FontSize = 18;
            app.R_minEditField.Position = [145 252 132 32];
            app.R_minEditField.Value = 4.5;

            % Create VolumeFractionEditFieldLabel
            app.VolumeFractionEditFieldLabel = uilabel(app.Panel);
            app.VolumeFractionEditFieldLabel.FontName = 'arial';
            app.VolumeFractionEditFieldLabel.FontSize = 16;
            app.VolumeFractionEditFieldLabel.Position = [16 200 125 32];
            app.VolumeFractionEditFieldLabel.Text = 'Volume  Fraction';

            % Create VolumefractionEditField
            app.VolumefractionEditField = uieditfield(app.Panel, 'numeric');
            app.VolumefractionEditField.Limits = [0 1];
            app.VolumefractionEditField.HorizontalAlignment = 'center';
            app.VolumefractionEditField.FontName = 'Arial';
            app.VolumefractionEditField.FontSize = 18;
            app.VolumefractionEditField.Position = [145 200 133 32];
            app.VolumefractionEditField.Value = 0.5;

            % Create Switch
            app.Switch = uiswitch(app.Panel, 'slider');
            app.Switch.Items = {'', ''};
            app.Switch.ItemsData = [1 0];
            app.Switch.Orientation = 'vertical';
            app.Switch.FontSize = 18;
            app.Switch.FontColor = [0 1 0];
            app.Switch.Position = [30 107 28 62];
            app.Switch.Value = 1;

            % Create AnisotropicdesignLabel
            app.AnisotropicdesignLabel = uilabel(app.Panel);
            app.AnisotropicdesignLabel.FontName = 'arial';
            app.AnisotropicdesignLabel.FontSize = 16;
            app.AnisotropicdesignLabel.Position = [67 99 136 22];
            app.AnisotropicdesignLabel.Text = 'Anisotropic design';

            % Create IsotropicdesignLabel
            app.IsotropicdesignLabel = uilabel(app.Panel);
            app.IsotropicdesignLabel.FontName = 'arial';
            app.IsotropicdesignLabel.FontSize = 16;
            app.IsotropicdesignLabel.Position = [68 159 117 22];
            app.IsotropicdesignLabel.Text = 'Isotropic design';

            % Create ConsoleOutput
            app.ConsoleOutput = uitextarea(app.GridLayout);
            app.ConsoleOutput.HorizontalAlignment = 'center';
            app.ConsoleOutput.Layout.Row = [34 36];
            app.ConsoleOutput.Layout.Column = [13 36];

            % Create BusyLampLabel
            app.BusyLampLabel = uilabel(app.GridLayout);
            app.BusyLampLabel.HorizontalAlignment = 'center';
            app.BusyLampLabel.FontSize = 18;
            app.BusyLampLabel.Layout.Row = [31 32];
            app.BusyLampLabel.Layout.Column = [34 35];
            app.BusyLampLabel.Text = 'Busy';

            % Create BusyLamp
            app.BusyLamp = uilamp(app.GridLayout);
            app.BusyLamp.Tag = 'Ctrl+ C for shut down';
            app.BusyLamp.Layout.Row = [31 32];
            app.BusyLamp.Layout.Column = 36;
            app.BusyLamp.Color = [0.702 0.702 0.702];

            % Create ComputegeometryButton
            app.ComputegeometryButton = uibutton(app.GridLayout, 'push');
            app.ComputegeometryButton.ButtonPushedFcn = createCallbackFcn(app, @ComputegeometryButtonPushed, true);
            app.ComputegeometryButton.FontName = 'Arial';
            app.ComputegeometryButton.FontSize = 16;
            app.ComputegeometryButton.FontWeight = 'bold';
            app.ComputegeometryButton.Layout.Row = [27 29];
            app.ComputegeometryButton.Layout.Column = [24 29];
            app.ComputegeometryButton.Text = {'2.Compute'; 'geometry'};

            % Create GenerateSTLfileButton
            app.GenerateSTLfileButton = uibutton(app.GridLayout, 'push');
            app.GenerateSTLfileButton.ButtonPushedFcn = createCallbackFcn(app, @GenerateSTLfileButtonPushed, true);
            app.GenerateSTLfileButton.FontName = 'Arial';
            app.GenerateSTLfileButton.FontSize = 16;
            app.GenerateSTLfileButton.FontWeight = 'bold';
            app.GenerateSTLfileButton.Layout.Row = [27 29];
            app.GenerateSTLfileButton.Layout.Column = [31 35];
            app.GenerateSTLfileButton.Text = {'3. Generate'; ' STL file'};

            % Create Tree
            app.Tree = uitree(app.GridLayout, 'checkbox');
            app.Tree.Layout.Row = [30 32];
            app.Tree.Layout.Column = [24 32];

            % Create AllstepsNode
            app.AllstepsNode = uitreenode(app.Tree);
            app.AllstepsNode.NodeData = 1;
            app.AllstepsNode.Text = 'All steps';

            % Create VariableDehomogenazationNode
            app.VariableDehomogenazationNode = uitreenode(app.AllstepsNode);
            app.VariableDehomogenazationNode.NodeData = [1 0];
            app.VariableDehomogenazationNode.Text = '1- Variable Dehomogenazation';

            % Create VoronoiseedingNode
            app.VoronoiseedingNode = uitreenode(app.AllstepsNode);
            app.VoronoiseedingNode.Text = '2- Voronoi seeding';

            % Create VoronoitessellationNode
            app.VoronoitessellationNode = uitreenode(app.AllstepsNode);
            app.VoronoitessellationNode.NodeData = 1;
            app.VoronoitessellationNode.Text = '3- Voronoi tessellation';

            % Create TrussinfocomputationNode
            app.TrussinfocomputationNode = uitreenode(app.AllstepsNode);
            app.TrussinfocomputationNode.Text = '4- Truss info computation';

            % Create BuildPolyshapeNode
            app.BuildPolyshapeNode = uitreenode(app.AllstepsNode);
            app.BuildPolyshapeNode.Text = '5- Build Polyshape';

            % Create BoundarySimplifyNode
            app.BoundarySimplifyNode = uitreenode(app.AllstepsNode);
            app.BoundarySimplifyNode.Text = '6- Boundary Simplify';

            % Assign Checked Nodes
            app.Tree.CheckedNodes = [app.VariableDehomogenazationNode, app.VoronoiseedingNode, app.VoronoitessellationNode, app.TrussinfocomputationNode, app.BuildPolyshapeNode, app.BoundarySimplifyNode, app.AllstepsNode];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = run_main_app

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end