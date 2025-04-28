% =========================================================================
% MATLAB Script for Automating HVDC Fault Simulations
% =========================================================================
clear; clc; close all;

load('xFinalState.mat');

% Generate fault scenarios
%faultScenarios = scenarioGenerator();

% print length of faultScenarios
%fprintf('Generated %d fault scenarios:\n', numel(faultScenarios));

% --- Simulation Setup ---
modelName = 'HVDCThyristorBased'; % *** REPLACE WITH YOUR ACTUAL .slx FILE NAME ***
Ts = 500e-6;                     % Discrete sample time (seconds) - KEEP THIS VALUE
SimulationStopTime = 3.5;       % Total simulation time (seconds) - Adjust as needed
saveMat = false;   % Set to false to skip saving .mat results
saveCSV = true;   % Set to false to skip saving ML-ready CSV

%% --- Output Folders ---
outputBaseFolder = 'HVDC_Simulation_Results';
if ~exist(outputBaseFolder, 'dir')
    mkdir(outputBaseFolder);
end
csvFolder = fullfile(outputBaseFolder, 'CSV');
if ~exist(csvFolder, 'dir')
    mkdir(csvFolder);
end

% --- Fault Scenario Definitions ---
% Define each simulation run as a struct within a cell array
% Add more scenarios by adding more structs to this cell array
if ~exist('faultScenarios','var')
    faultScenarios = {
        struct('id', 'NoFault',       'fault_type', 'None', 'fault_resistance', 0.01, 'dc_fault_location_pu', 0.5);
        struct('id', 'DC_Fault_25pct', 'fault_type', 'DC',   'fault_resistance', 0.01, 'dc_fault_location_pu', 0.25);
    };
end

shortNameMap = struct( ...
    "FaultSignals_IDCFault_A_", "DCFaultCurrent", ...
    "FaultSignals_IA_GFault_A_", "FaultCurrent_AG_PhA", ...
    "FaultSignals_IB_GFault_A_", "FaultCurrent_BG_PhB", ...
    "FaultSignals_IC_GFault_A_", "FaultCurrent_CG_PhC", ...
    "FaultSignals_IA_BFault_A_", "FaultCurrent_AB_PhA", ...
    "FaultSignals_IA_CFault_A_", "FaultCurrent_AC_PhA", ...
    "FaultSignals_IB_CFault_A_", "FaultCurrent_BC_PhB", ...
    "RectifierACSignals_Vabc_pu__1", "Rectifier_Va_pu", ... % tt
    "RectifierACSignals_Vabc_pu__2", "Rectifier_Vb_pu", ...
    "RectifierACSignals_Vabc_pu__3", "Rectifier_Vc_pu", ...
    "RectifierACSignals_Iabc_pu_100MVA__1", "Rectifier_Ia_pu", ...
    "RectifierACSignals_Iabc_pu_100MVA__2", "Rectifier_Ib_pu", ...
    "RectifierACSignals_Iabc_pu_100MVA__3", "Rectifier_Ic_pu", ...
    "RectifierControlSignals_VdL_pu_", "Rectifier_VdL_pu", ... % tt
    "RectifierControlSignals_IdIdref_lim_pu__1", "Rectifier_Id_pu", ...
    "RectifierControlSignals_IdIdref_lim_pu__2", "Rectifier_IdrefLim_pu", ...
    "RectifierControlSignals_alpha_ord_deg_", "Rectifier_AlphaOrd_deg", ...
    "RectifierControlSignals_ControlMode_0_blocked1_current2_voltage", "Rectifier_ControlMode", ...
    "RectifierProtectionStatus_Low_ac_volt_R", "Rectifier_LowACVolt", ... % tt
    "RectifierProtectionStatus_Forced_alpha_R", "Rectifier_ForcedAlpha", ...
    "RectifierValveSignals_uSw1_V_", "Rectifier_Valve1_Voltage", ... % tt
    "RectifierValveSignals_iSw1_Sw3_A__1", "Rectifier_Valve1_Current", ...
    "RectifierValveSignals_iSw1_Sw3_A__2", "Rectifier_Valve3_Current", ...
    "RectifierValveSignals_alpha_ord_deg_", "Rectifier_Valve_AlphaOrd", ...
    "InverterACSignals_Vabc_pu__1", "Inverter_Va_pu", ... % tt
    "InverterACSignals_Vabc_pu__2", "Inverter_Vb_pu", ...
    "InverterACSignals_Vabc_pu__3", "Inverter_Vc_pu", ...
    "InverterACSignals_Iabc_pu_100MVA__1", "Inverter_Ia_pu", ...
    "InverterACSignals_Iabc_pu_100MVA__2", "Inverter_Ib_pu", ...
    "InverterACSignals_Iabc_pu_100MVA__3", "Inverter_Ic_pu", ...
    "InverterControlSignals_VdLVd_ref_pu__1", "Inverter_VdL_pu", ... % tt
    "InverterControlSignals_VdLVd_ref_pu__2", "Inverter_VdRef_pu", ...
    "InverterControlSignals_IdIdref_lim_pu__1", "Inverter_Id_pu", ...
    "InverterControlSignals_IdIdref_lim_pu__2", "Inverter_IdrefLim_pu", ...
    "InverterControlSignals_alpha_ord_deg_", "Inverter_AlphaOrd_deg", ...
    "InverterControlSignals_ControlMode_0_blocked1_current2_voltage3", "Inverter_ControlMode", ...
    "InverterControlSignals_gamma_meanGamma_ref_deg__1", "Inverter_GammaMean_deg", ...
    "InverterControlSignals_gamma_meanGamma_ref_deg__2", "Inverter_GammaRef_deg", ...
    "InverterProtectionStatus_Low_ac_volt_I", "Inverter_LowACVolt", ... % tt
    "InverterProtectionStatus_A_min_I", "Inverter_AMin", ...
    "InverterValve1Signals_uSw1_V_", "Inverter_Valve1_Voltage", ... % tt
    "InverterValve1Signals_iSw1Ucom1_pu__1", "Inverter_Valve1_Current", ...
    "InverterValve1Signals_iSw1Ucom1_pu__2", "Inverter_Ucom1_Current", ...
    "InverterValve1Signals_gammaMean_deg_", "Inverter_Valve_GammaMean_deg" ...
);


% --- Common Fault Timing ---
% Define when faults generally occur (can be overridden in scenario struct if needed)
default_fault_time_start = 0.5; % Start fault after system stabilizes
default_fault_duration   = 0.1; % e.g., 6 cycles @ 60Hz or 5 cycles @ 50Hz

% --- Constants ---
L_total_km = 300; % Total DC Line length

fprintf('Starting HVDC fault simulation batch...\n');

% =========================================================================
% Main Simulation Loop
% =========================================================================
for i = 1:length(faultScenarios)
    scenario = faultScenarios{i};
    fprintf('Running Scenario %d: %s\n', i, scenario.id);

    % --- Prepare Parameters for this Run ---
    fault_type           = scenario.fault_type;
    fault_resistance     = scenario.fault_resistance;
    dc_fault_location_pu = scenario.dc_fault_location_pu;

    % Use default timing unless specified in scenario
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

    % Calculate switching times
    fault_switch_times = [fault_time_start, fault_time_start + fault_duration];
    no_fault_switch_times = [SimulationStopTime + 1, SimulationStopTime + 2]; % Ensures breaker stays open

    % Calculate DC line lengths
    dc_line_X_km = dc_fault_location_pu * L_total_km;
    dc_line_Y_km = (1 - dc_fault_location_pu) * L_total_km;

    % --- Assign Variables to Base Workspace for Simulink ---
    % Assign Ts (critical for discrete blocks)
    assignin('base', 'Ts', Ts);

    % Assign DC line lengths (Simulink blocks read these variables)
    assignin('base', 'DcFaultDX', dc_line_X_km);
    assignin('base', 'DcFaultDY', dc_line_Y_km);

    % Initialize all fault times to 'no fault' state
    assignin('base', 'DCFaultTime', no_fault_switch_times); % Changed name
    assignin('base', 'AGFaultTime', no_fault_switch_times); % Changed name
    assignin('base', 'BGFaultTime', no_fault_switch_times); % Changed name
    assignin('base', 'CGFaultTime', no_fault_switch_times); % Changed name
    assignin('base', 'ABFaultTime', no_fault_switch_times); % Changed name
    assignin('base', 'ACFaultTime', no_fault_switch_times); % Changed name
    assignin('base', 'BCFaultTime', no_fault_switch_times); % Changed name

    % Assign fault resistances (Simulink blocks read these)
    % Assuming all breaker blocks use the same variable name structure
    assignin('base', 'DCFaultRes', fault_resistance); % Changed name
    assignin('base', 'AGFaultRes', fault_resistance); % Changed name
    assignin('base', 'BGFaultRes', fault_resistance); % Changed name
    assignin('base', 'CGFaultRes', fault_resistance); % Changed name
    assignin('base', 'ABFaultRes', fault_resistance); % Changed name
    assignin('base', 'ACFaultRes', fault_resistance); % Changed name
    assignin('base', 'BCFaultRes', fault_resistance); % Changed name

    % Activate the selected fault by setting its specific switch times
    switch fault_type
        case 'None'
            % All times already set to no_fault_switch_times
        case 'DC'
            assignin('base', 'DCFaultTime', fault_switch_times); % Changed name
        case 'AG'
            assignin('base', 'AGFaultTime', fault_switch_times); % Changed name
        case 'BG'
            assignin('base', 'BGFaultTime', fault_switch_times); % Changed name
        case 'CG'
            assignin('base', 'CGFaultTime', fault_switch_times); % Changed name
        case 'AB'
            assignin('base', 'ABFaultTime', fault_switch_times); % Changed name
        case 'AC'
            assignin('base', 'ACFaultTime', fault_switch_times); % Changed name
        case 'BC'
            assignin('base', 'BCFaultTime', fault_switch_times); % Changed name
        otherwise
            warning('Scenario %d: Unknown fault_type specified: %s. Skipping fault activation.', i, fault_type);
    end

    % --- Run Simulation ---
    try
        fprintf('   Simulating...');
        % Use sim command, capture output in simOut object
        simOut = sim(modelName, 'SimulationMode', 'accelerator', 'StopTime', num2str(SimulationStopTime),'SaveOutput', 'on', 'OutputSaveName', 'simDataLog'); % Log signals to simOut.simDataLog
        fprintf(' Done.\n');

        % --- Process and Save Results ---
        fprintf('   Saving results...');

        % Create a structure to hold the results for this run
        Results = struct();
        Results.Scenario = scenario; % Store the parameters used for this run
        Results.SimulationOutput = struct(); % Store the actual logged data
        Results.MergedData = Simulink.SimulationData.Dataset; % Initialize an empty Dataset for merged signals

        % --- IMPORTANT: Adapt these names based on your Scope logging ---
        % Check the simOut.simDataLog structure for the exact variable names
        % These names MUST match the 'Variable name' you set in Scope properties -> Logging
        loggedDataNames = {'FaultSignals', 'RectifierACSignals', 'RectifierControlSignals', ...
                           'RectifierProtectionStatus', 'RectifierValveSignals', 'InverterACSignals', ...
                           'InverterControlSignals', 'InverterProtectionStatus', 'InverterValve1Signals'};
        
        for k = 1:length(loggedDataNames)
            dataName = loggedDataNames{k}; % e.g., 'FaultCurrents'
            % Check if simOut has this property
            if isprop(simOut, dataName) || isfield(simOut, dataName) % Check both for robustness
                currentDataset = simOut.(dataName);

                % Verify it's the expected type (Dataset or Signal)
                if isa(currentDataset, 'Simulink.SimulationData.Dataset')
                    signalNames = currentDataset.getElementNames(); % Get names of signals within this Dataset

                    for j = 1:length(signalNames)
                        signalName = signalNames{j};
                        try
                            % Combine dataName and signalName
                            mergedName = sprintf('%s_%s', dataName, signalName);

                            % Get short name from the dictionary
                            if isfield(shortNameMap, mergedName)
                                shortName = shortNameMap.(mergedName);
                            else
                                % shortName = matlab.lang.makeValidName(mergedName); % Fallback to a valid MATLAB name
                                shortName = mergedName;
                            end

                            signalObj = currentDataset.get(signalName);

                            if isa(signalObj, 'Simulink.SimulationData.Signal')
                                % Add the individual signal to the merged dataset with the new unique name
                                Results.MergedData = Results.MergedData.addElement(signalObj, shortName);
                            else
                                warning('Scenario %d: Element "%s" within Dataset "%s" is not a Signal object.', i, signalName, dataName);
                            end

                        catch ME_getSignal
                                warning('Scenario %d: Could not get element "%s" from Dataset "%s". Error: %s', i, signalName, dataName, ME_getSignal.message);
                        end
                    end
                    % fprintf('    -> Stored data for "%s".\n', dataName);
                else
                    warning('Scenario %d: Output field "%s" is not a recognized Dataset.', i, dataName);
                end

            else
                warning('Scenario %d: Logged data property/field "%s" not found directly in simOut object. Check Scope logging/variable name.', i, dataName);
            end
        end
        
        % Combine all signals into a single table
        combinedTable = table();
        for j = 1:numElements(Results.MergedData)
            signal = Results.MergedData.get(j);
            ts = signal.Values;
            if isempty(ts.Time)
                continue
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
        
            if ~ismember('Time', combinedTable.Properties.VariableNames)
                combinedTable.Time = t;
            end
        
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
        
        if saveMat
            outputFileName = sprintf('Results_Run%03d_%s.mat', i, scenario.id);
            outputFilePath = fullfile(outputBaseFolder, outputFileName);
            save(outputFilePath, 'Results');
            fprintf(' Saved to %s\n', outputFileName);
        end
        
        if saveCSV
            csvFileName = sprintf('ML_Run%03d_%s.csv', i, scenario.id);
            writetable(combinedTable, fullfile(csvFolder, csvFileName));
            fprintf('    -> ML CSV saved: %s\n', csvFileName);
        end        

    catch ME % Catch simulation errors
        warning('Scenario %d (%s) failed: %s', i, scenario.id, ME.message);
        fprintf('   Skipping saving for this scenario.\n');
        % Display detailed error information
        disp(ME.getReport);
    end

    % Optional: Clear base workspace variables specific to the run if needed
    clearvars -except modelName Ts SimulationStopTime outputBaseFolder faultScenarios default_* L_total_km i scenario loggedDataNames csvFolder xFinalState saveMat saveCSV faultScenarios shortNameMap % Keep loop variables etc.

end % End of simulation loop

fprintf('Batch simulation complete.\n');