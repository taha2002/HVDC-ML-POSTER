% =========================================================================
% MATLAB Script for Automating HVDC Fault Simulations in Parallel with Saving Files
% =========================================================================
clear; clc; close all;
%clc;

% Load necessary data (ensure that xFinalState.mat exists)
load('xFinalState.mat');

% Generate fault scenarios
faultScenarios = scenarioGenerator(); % Ensure this function is available

% print length of faultScenarios
fprintf('Generated %d fault scenarios:\n', numel(faultScenarios));

%error("force stop");

%% --- Simulation Setup ---
modelName = 'HVDCThyristorBased';  % Replace with your actual .slx file name
load_system(modelName);

Ts = 50e-6;                        % Discrete sample time (seconds)
SimulationStopTime = 3.5;          % Simulation stop time (3.5s: 1.5s stabilization + 2s study)
saveMat = false;                    % Set to true to save .mat results
saveCSV = true;                    % Set to true to save CSV for ML

%% --- Output Folders ---
outputBaseFolder = 'HVDC_Simulation_Results';
if ~exist(outputBaseFolder, 'dir')
    mkdir(outputBaseFolder);
end
csvFolder = fullfile(outputBaseFolder, 'CSV');
if ~exist(csvFolder, 'dir')
    mkdir(csvFolder);
end

% Create subfolders for each unique fault type
fixedFaultTypes = {'None', 'DC', 'AG', 'BG', 'CG', 'AB', 'AC', 'BC'};

csvFoldersByFaultType = containers.Map();
for k = 1:length(fixedFaultTypes)
    ft = fixedFaultTypes{k};
    folderPath = fullfile(csvFolder, ft);
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end
    csvFoldersByFaultType(ft) = folderPath;
end

%% --- Fault Scenario Definitions ---
% If faultScenarios variable is missing, define default scenarios.
if ~exist('faultScenarios','var') || isempty(faultScenarios)
    faultScenarios = {
        struct('id', 'NoFault',       'fault_type', 'None', 'fault_resistance', 0.01, 'dc_fault_location_pu', 0.5);
        struct('id', 'DC_Fault_25pct', 'fault_type', 'DC',   'fault_resistance', 0.01, 'dc_fault_location_pu', 0.25);
    };
end

%% --- Common Fault Timing and Constants ---
default_fault_time_start = 2.5;  % Default fault start time
default_fault_duration   = 0.1;  % Default fault duration (e.g., 6 cycles @ 60Hz)
L_total_km = 300;              % Total DC line length in km

fprintf('Starting HVDC fault simulation batch in parallel...\n');

%% --- Prepare SimulationInput Array for Parallel Runs ---
numScenarios = numel(faultScenarios);
%simInputs(numScenarios) = Simulink.SimulationInput(modelName);  % Preallocate array
simInputs(1:numScenarios) = Simulink.SimulationInput(modelName);

for i = 1:numScenarios
    scenario = faultScenarios{i};
    
    % Set basic simulation variables.
    simInputs(i) = simInputs(i).setVariable('Ts', Ts);
    simInputs(i) = simInputs(i).setModelParameter('SimulationMode', 'accelerator');
    % simInputs(i) = simInputs(i).setVariable('xFinalState', xFinalState);
    
    % Calculate and set DC line lengths.
    dc_line_X_km = scenario.dc_fault_location_pu * L_total_km;
    dc_line_Y_km = (1 - scenario.dc_fault_location_pu) * L_total_km;
    simInputs(i) = simInputs(i).setVariable('DcFaultDX', dc_line_X_km);
    simInputs(i) = simInputs(i).setVariable('DcFaultDY', dc_line_Y_km);
    
    % Use scenario timing if provided; otherwise, use defaults.
    if isfield(scenario, 'fault_time_start')
        fault_time_start = scenario.fault_time_start;
    else
        fault_time_start = default_fault_time_start;
    end
    if isfield(scenario, 'fault_duration')
        fault_duration = scenario.fault_duration;
    else
        fault_duration = default_fault_duration;
    end
    fault_switch_times = [fault_time_start, fault_time_start + fault_duration];
    no_fault_switch_times = [SimulationStopTime + 1, SimulationStopTime + 2];  % Ensures breaker remains open
    
    % Set all fault timing variables to a "no fault" state initially.
    simInputs(i) = simInputs(i).setVariable('DCFaultTime', no_fault_switch_times);
    simInputs(i) = simInputs(i).setVariable('AGFaultTime', no_fault_switch_times);
    simInputs(i) = simInputs(i).setVariable('BGFaultTime', no_fault_switch_times);
    simInputs(i) = simInputs(i).setVariable('CGFaultTime', no_fault_switch_times);
    simInputs(i) = simInputs(i).setVariable('ABFaultTime', no_fault_switch_times);
    simInputs(i) = simInputs(i).setVariable('ACFaultTime', no_fault_switch_times);
    simInputs(i) = simInputs(i).setVariable('BCFaultTime', no_fault_switch_times);
    
    % Set fault resistances (assumed same for all types).
    simInputs(i) = simInputs(i).setVariable('DCFaultRes', scenario.fault_resistance);
    simInputs(i) = simInputs(i).setVariable('AGFaultRes', scenario.fault_resistance);
    simInputs(i) = simInputs(i).setVariable('BGFaultRes', scenario.fault_resistance);
    simInputs(i) = simInputs(i).setVariable('CGFaultRes', scenario.fault_resistance);
    simInputs(i) = simInputs(i).setVariable('ABFaultRes', scenario.fault_resistance);
    simInputs(i) = simInputs(i).setVariable('ACFaultRes', scenario.fault_resistance);
    simInputs(i) = simInputs(i).setVariable('BCFaultRes', scenario.fault_resistance);
    
    % Activate the selected fault type by overriding the timing variable.
    switch scenario.fault_type
        case 'None'
            % Leave all fault times at no_fault_switch_times.
        case 'DC'
            simInputs(i) = simInputs(i).setVariable('DCFaultTime', fault_switch_times);
        case 'AG'
            simInputs(i) = simInputs(i).setVariable('AGFaultTime', fault_switch_times);
        case 'BG'
            simInputs(i) = simInputs(i).setVariable('BGFaultTime', fault_switch_times);
        case 'CG'
            simInputs(i) = simInputs(i).setVariable('CGFaultTime', fault_switch_times);
        case 'AB'
            simInputs(i) = simInputs(i).setVariable('ABFaultTime', fault_switch_times);
        case 'AC'
            simInputs(i) = simInputs(i).setVariable('ACFaultTime', fault_switch_times);
        case 'BC'
            simInputs(i) = simInputs(i).setVariable('BCFaultTime', fault_switch_times);
        otherwise
            warning('Scenario %d: Unknown fault_type: %s', i, scenario.fault_type);
    end
    
    % Optionally, store scenario ID for later use.
    simInputs(i) = simInputs(i).setVariable('scenarioID', scenario.id);
    
    % Set simulation stop time for this run.
    simInputs(i) = simInputs(i).setModelParameter('StopTime', num2str(SimulationStopTime));
end

%Simulink.BlockDiagram.buildRapidAcceleratorTarget(modelName);

%% --- Run Simulations in Parallel ---
% parsim executes the array of SimulationInput objects in parallel.
simOut = parsim(simInputs,'TransferBaseWorkspaceVariables','on','UseFastRestart','on');

%% --- Postprocessing and Saving Results ---
% Define the list of logged data names from your Simulink model.
loggedDataNames = {'FaultSignals', 'RectifierACSignals', 'RectifierControlSignals', ...
                   'RectifierProtectionStatus', 'RectifierValveSignals', 'InverterACSignals', ...
                   'InverterControlSignals', 'InverterProtectionStatus', 'InverterValve1Signals'};

for i = 1:numScenarios
    scenario = faultScenarios{i};
    fprintf('Processing Scenario %d: %s\n', i, scenario.id);
    simV2Data = simOut(i);
    
    % Create a structure to store results.
    Results = struct();
    Results.Scenario = scenario;
    Results.SimulationOutput = simData;
    Results.MergedData = Simulink.SimulationData.Dataset;  % Initialize empty dataset
    
    % Loop through each logged data group.
    for k = 1:length(loggedDataNames)
        dataName = loggedDataNames{k};
        % Check if the field exists in the simulation output.
        if isprop(simData, dataName) || isfield(simData, dataName)
            try
                currentDataset = simData.(dataName);
                if isa(currentDataset, 'Simulink.SimulationData.Dataset')
                    signalNames = currentDataset.getElementNames();
                    for j = 1:length(signalNames)
                        try
                            signalObj = currentDataset.get(signalNames{j});
                            if isa(signalObj, 'Simulink.SimulationData.Signal')
                                Results.MergedData = Results.MergedData.addElement(signalObj, signalObj.Name);
                            else
                                warning('Scenario %d: Element "%s" in %s is not a Signal object.', i, signalNames{j}, dataName);
                            end
                        catch ME_signal
                            warning('Scenario %d: Error accessing signal "%s" in %s: %s', i, signalNames{j}, dataName, ME_signal.message);
                        end
                    end
                    % fprintf('    -> Stored data for "%s".\n', dataName);
                else
                    warning('Scenario %d: Field "%s" is not a Dataset.', i, dataName);
                end
            catch ME_dataset
                warning('Scenario %d: Failed to extract "%s": %s', i, dataName, ME_dataset.message);
            end
        else
            warning('Scenario %d: Field "%s" not found in simulation output.', i, dataName);
        end
    end
    
    % Combine signals from the merged dataset into a single table.
    combinedTable = table();
    numSignals = Results.MergedData.numElements;
    for j = 1:numSignals
        signal = Results.MergedData.get(j);
        ts = signal.Values;
        if isempty(ts.Time)
            continue;
        end

        % Trim to t >= 2s only
        t = ts.Time;
        data = ts.Data;
        validIdx = t >= 2;
        t = t(validIdx);
        data = data(validIdx, :);

        if isempty(t)
            continue;
        end
        
        % If the table doesn't already contain the time vector, add it.
        if ~ismember('Time', combinedTable.Properties.VariableNames)
            combinedTable.Time = t;
        end
        
        % Add data columns
        if size(data, 2) == 1
            colName = matlab.lang.makeValidName(signal.Name);
            combinedTable.(colName) = data;
        else
            for col = 1:size(data, 2)
                newColName = sprintf('%s_%d', matlab.lang.makeValidName(signal.Name), col);
                combinedTable.(newColName) = data(:, col);
            end
        end
    end
    
    if ~isempty(simOut(i).ErrorMessage)
        fprintf('Error in Scenario %d (%s): %s\n', i, scenario.id, simOut(i).ErrorMessage);
    end

    % Save the results structure as a .mat file.
    if saveMat
        outputFileName = sprintf('Results_Run%03d_%s.mat', i, scenario.id);
        outputFilePath = fullfile(outputBaseFolder, outputFileName);
        save(outputFilePath, 'Results');
        fprintf('Saved MAT file: %s\n', outputFileName);
    end
    
    % Save the combined table as a CSV file.
    if saveCSV
        if isempty(combinedTable) || height(combinedTable) == 0
            warning('Scenario %d: combinedTable is empty. Skipping CSV save.', i);
            continue;
        end

        csvFileName = sprintf('ML_Run%03d_%s.csv', i, scenario.id);
        subfolderPath = csvFoldersByFaultType(scenario.fault_type);  % Get the folder for this fault type
        csvFilePath = fullfile(subfolderPath, csvFileName);

        maxRetries = 3;
        retryCount = 0;
        success = false;

        while ~success && retryCount < maxRetries
            try
                writetable(combinedTable, csvFilePath);
                fprintf('Saved CSV file: %s\n', csvFileName);
                success = true;
            catch ME
                retryCount = retryCount + 1;
                warning('Failed to save CSV (Attempt %d/%d) for scenario %d: %s', ...
                    retryCount, maxRetries, i, ME.message);
                pause(1);  % Wait 1 second before retrying
            end
        end

        if ~success
            warning('Scenario %d: Failed to save CSV after %d attempts.', i, maxRetries);
        end

    end
end

fprintf('Parallel batch simulation complete.\n');
