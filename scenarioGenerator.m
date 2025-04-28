function faultScenarios = scenarioGenerator()
% =========================================================================
% Generate Fault Scenarios for HVDC Fault Diagnosis ML Training
% =========================================================================

% --- Define Fault Types to Include ---
% None: No fault (healthy system)
% DC: DC line faults at different locations
% AG/BG/CG: Single phase to ground faults
% AB/BC/AC: Phase to phase faults
% ABC: Three phase fault
faultTypes = {'None', ...                    % Healthy system
    'DC10', 'DC25', 'DC50', 'DC75', 'DC90'...   % DC faults at different locations
    'AG', 'BG', 'CG',...               % Single phase to ground
    'AB', 'BC', 'AC',...               % Phase to phase
    %'ABC'                           % Three phase
    };

% --- Define Parameter Arrays ---
% Fault times centered at 2.5s (1.5s stabilization + 1s into study period)
faultTimes = [2.5];  % Fixed time for consistent analysis

% --- Define Parameter Arrays ---
% Number of repetitions for each scenario
numRepetitions = [20, 30];

% Fault durations (in seconds)
% - 0.05s (2.5 cycles @ 50Hz): Very fast clearing
% - 0.1s  (5 cycles @ 50Hz): Normal clearing
% - 0.15s (7.5 cycles @ 50Hz): Slow clearing
% - 0.2s  (10 cycles @ 50Hz): Very slow clearing
% faultDurations = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5];
faultDurations = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5];

% Fault resistances (in ohms)
% - 0.01: Very low impedance fault
% - 0.1: Low impedance fault
% - 1.0: Medium impedance fault
% - 10.0: High impedance fault
% - 50.0: Very high impedance fault
% faultResistances = [0.01, 0.1, 1.0, 10.0, 50.0];
faultResistances = [0.01, 0.1, 1, 10, 50, 100, 500];

% --- Initialize Scenario Cell Array ---
faultScenarios = {};   % Each scenario will be stored as a struct in this cell array
scenarioIndex  = 1;

% --- Generate Scenarios with Random Repetitions ---
for ft = 1:length(faultTypes)
    currentFaultType = faultTypes{ft};
    for t = faultTimes
        for rep = 1:randi(numRepetitions)
            for d = faultDurations
                for r = faultResistances
                    randomOffset = 0.5 * rand;
                    
                    % Create a unique scenario ID for identification
                    scenarioID = sprintf('%s_t%.2f_d%.2f_r%.3f', currentFaultType, t, d, r);
                    
                    % Create scenario struct with base parameters
                    scenario = struct();
                    scenario.id = scenarioID;
                    scenario.fault_type = currentFaultType;
                    scenario.fault_time_start = t + randomOffset;
                    scenario.fault_duration = d;
                    scenario.fault_resistance = r;
                    
                    % For DC fault types, extract the numeric value if available.
                    if contains(currentFaultType, 'DC')
                        scenario.fault_type = 'DC';
                        numericPart = regexp(currentFaultType, '\d+', 'match');
                        if ~isempty(numericPart)
                            % Convert percentage to per unit (e.g., 'DC25' becomes 0.25)
                            scenario.dc_fault_location_pu = str2double(numericPart{1}) / 100;
                        else
                            % If fault type is just 'DC', use default of 50%
                            scenario.dc_fault_location_pu = 0.5;
                        end
                    else
                        % For non-DC faults, assign a default fault location (if applicable)
                        scenario.dc_fault_location_pu = 0.5;
                    end
                    
                    % Store the scenario in the cell array
                    faultScenarios{scenarioIndex} = scenario;
                    scenarioIndex = scenarioIndex + 1;
                end
            end
        end
    end
end

% --- Display All Generated Scenarios ---
% fprintf('Generated %d fault scenarios:\n', numel(faultScenarios));
% for i = 1:numel(faultScenarios)
%     disp(faultScenarios{i});
% end
end