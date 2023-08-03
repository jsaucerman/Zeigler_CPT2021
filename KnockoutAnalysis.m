% Written by Anirudha Chandrabhatla 
% Last Updated 11/24/2020
% Version 3.0

% % If any inputs have changed, run the following section: 
% 
%  warning off;
% % Species information from fib617_references.xlsx, 'species' tab 
% nodeGenesTable = readtable('fib617_references.xlsx', 'Sheet', 'species');
% 
% % Reaction information from fib417_references.xlsx, 'reactions' tab 
% networkReactions = readtable('fib617_references.xlsx', 'Sheet', 'reactions');
% 
% % Drug information from AllPharmActive.xlsx
% targetsTable = readtable('AllPharmActive.xlsx');
% %targetsTable = readtable('SmallMolec_pharmacologically_active_918.xlsx');
% 
% % Drug ID to drug name matching from drugIDsAndNamesFull516.xlsx
% %drugIDToNameConversion = readtable('drugIDsAndNamesFull516.xlsx');
% drugIDToNameConversion = readtable('drugIDsAndNamesFull918.xlsx');
% 
% % Manually curated competitive/non-competitive classifications for all drugs that target a node in the network
% compNonCompClassification = readtable('All Drugs Comp NonComp Classification.xlsx', 'Sheet', 'All Drugs Condensed');
% 
% % Output from the Webscraper. Lists the drug name, target(s) and action(s) (Agonist / Antagonist)
% drugInformation = readtable('drugTargetsAndActions_918.xlsx');
% 
% %load('finalDrugOutputNetworkTargets.mat'); % Result from the webscraper NEED TO CLEAN UP WEBSCRAPER OUTPUT AND THEN CHANGE THIS TO BE READTABLE LIKE THE OTHER VARIABLES
% %drugInformation = finalDrugOutputNetworkTargets;
% 
% geneNames = nodeGenesTable.('geneName');  % Gene names from the model
% geneIDs = nodeGenesTable.('ID');          % Gene IDs from the model
% 
% % Drug IDs from the database of approved targets
% drugIDs_fromTargetList = targetsTable.('DrugIDs'); 
% % Gene names of the drug targets. Retrieved from the database of approved targets.
% drugGeneTargets = targetsTable.('GeneName'); 
% 
% % Drug names from the database of all drugs
% drugNames = drugIDToNameConversion.('Name');
% % Drug IDs from the database of all drugs 
% drugIDs_fromTargetList_forNameConversion = drugIDToNameConversion.('DrugBankID');
% 
% [drugsToSimulate, formattedReactions] = formatInputs(geneNames, geneIDs, drugGeneTargets, drugIDs_fromTargetList, drugIDs_fromTargetList_forNameConversion, drugNames, networkReactions, compNonCompClassification, drugInformation);

%% Inputs for simulations

drugsToSimulate = load('drugsToSimulate.mat');
drugsToSimulate = drugsToSimulate.drugsToSimulate; % Extract from struct

formattedReactions = load('formattedReactions.mat');
formattedReactions = formattedReactions.formattedReactions; % Extract from struct

clusteredRowLabels = load('clusteredRowLabels.mat');
clusteredRowLabels = clusteredRowLabels.clusteredRowLabels; % Extract from struct

% Set drug dose or doses (as a vector)
dose_antag = [0.85]; 
dose_ag = [0.85];

%% Creating Input Curves
[InputCsim,tInSim,inputNode,resNorm,resNormConvert] = InputCurve(0.6, 0.6);
%% Part 2: Generate graph of results

% Parameters and initial conditions
[params,y0] = fib617_params;
tspan = [0 2329]; 
options = [];
timePointsOfInterest = [168,192,336,1176];
% Pull parameters out in order to alter them in the loop without negating the previous alteration
[rpar,tau,ymax,speciesNames,KI]=params{:}; 
w = rpar(1,:);
wNew = rpar(1,:); % ONLY NEED WHEN CHANGING INPUT MECHANICAL WEIGHT
wNew(3) = 0.6; % Altering mechanical input
wNew(23) = 0.6; % Altering the MMP9 & TGFBmRNA => TGFB weight
wNew(24) = 0.6; % Altering the MMP2 & TGFBmRNA => TGFB weight
wNew(25) = 0; % Altering the TGFBmRNA => TGFB weight. 
n = rpar(2,:);
EC50 = rpar(3,:);
dose = rpar(4,:);
drugBinding = rpar(5,:); 
drugAgonism = rpar(6,:); 

%rpar = [wNew;n;EC50;alt;drugType]; % ONLY NEED WHEN CHANGING INPUT MECHANICAL WEIGHT
rpar = [wNew;n;EC50;dose;drugBinding;drugAgonism];

% Default simulation creating a vector showing steady state activation of all species
params = {rpar,tau,ymax,speciesNames,KI,InputCsim,inputNode,tInSim};
[t,y] = ode15s(@dynamicODE,tspan,y0,options,params);


y1 = real(interp1(t,y,tInSim));
yEnd0 = y1(timePointsOfInterest,:);

% Preturbed simulations
sens = zeros(length(timePointsOfInterest),length(speciesNames), height(drugsToSimulate));

drugTargetResponse = zeros(1, height(drugsToSimulate)); 
drugTargetResponseID = zeros(1, height(drugsToSimulate));
colImRNA_ANI = [];
colImRNA_ANI = [colImRNA_ANI, y1(:,101)];
prolif_ANI = [];
prolif_ANI = [prolif_ANI, y1(:,84)];

knockoutData_drug = table('Size', [108 2], 'VariableTypes', {'string', 'double'});

%% Drug knockouts
% Iterate through 'drugsToSimulate' and simulate the addition of each drug into the network.
for i = 32 % Input the row number (in 'drugsToSimulate.mat') of the drug you want to analyze.
    
    disp(drugsToSimulate{i, 1}); % Prints drug name to the command window 
    for r = 0:length(speciesNames)
        disp(num2str(r));
        ymax_Knockout = ymax;
        if r > 0
            ymax_Knockout(r) = 0;
        end
        for j = 1:length(dose_ag)
        
            doseNew = dose;
            drugBindingNew = drugBinding;
            drugAgonismNew = drugAgonism;

            % The drug is an agonist
            if strcmp(drugsToSimulate.IsAgonist{i}, 'Yes') == 1
                if isempty(find(drugsToSimulate.AgonistTarget{i} == ';')) % Drug has one agonist target
                    locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, drugsToSimulate.AgonistTargetGeneID{i});
                    if strcmpi('Competitive', drugsToSimulate.DrugAction{i}) % Competitive, agonist
                        drugBindingNew(locationOfReactions) = 1;
                        drugAgonismNew(locationOfReactions) = 1;
                        doseNew(locationOfReactions) = -1*dose_ag(j);
                    else % Non-Competitive, agonist
                        drugBindingNew(locationOfReactions) = -1;
                        drugAgonismNew(locationOfReactions) = 1;
                        doseNew(locationOfReactions) = dose_ag(j);
                    end
                else % Drug has multiple targets
                    geneIDsOfTargets = strsplit(drugsToSimulate.AgonistTargetGeneID{i}, ';');
                    for k = 1:length(geneIDsOfTargets)
                        locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, geneIDsOfTargets{k});
                        if strcmpi('Competitive', drugsToSimulate.DrugAction{i})
                            drugBindingNew(locationOfReactions) = 1;
                            drugAgonismNew(locationOfReactions) = 1;
                            doseNew(locationOfReactions) = -1*dose_ag(j);
                        else 
                            drugBindingNew(locationOfReactions) = -1;
                            drugAgonismNew(locationOfReactions) = 1;
                            doseNew(locationOfReactions) = dose_ag(j);
                        end
                    end
                end
            end

            % The drug is an antagonist
            if strcmp(drugsToSimulate.IsAntagonist{i}, 'Yes') == 1
                if isempty(find(drugsToSimulate.AntagonistTarget{i} == ';')) % Drug has one antagonist target
                    locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, drugsToSimulate.AntagonistTargetGeneID{i});
                    if strcmpi('Competitive', drugsToSimulate.DrugAction{i}) % Competitive, antagonist
                        drugBindingNew(locationOfReactions) = 1;
                        drugAgonismNew(locationOfReactions) = -1;
                        doseNew(locationOfReactions) = dose_antag(j);
                    else
                        drugBindingNew(locationOfReactions) = -1;
                        drugAgonismNew(locationOfReactions) = -1;
                        doseNew(locationOfReactions) = dose_antag(j);
                    end
                else
                    geneIDsOfTargets = strsplit(drugsToSimulate.AntagonistTargetGeneID{i}, ';');
                    for p = 1:length(geneIDsOfTargets)
                        locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, geneIDsOfTargets{p});
                        if strcmpi('Competitive', drugsToSimulate.DrugAction{i})
                            drugBindingNew(locationOfReactions) = 1;
                            drugAgonismNew(locationOfReactions) = -1;
                            doseNew(locationOfReactions) = dose_antag(j);
                        else 
                            drugBindingNew(locationOfReactions) = -1;
                            drugAgonismNew(locationOfReactions) = -1;
                            doseNew(locationOfReactions) = dose_antag(j);
                        end
                    end
                end
            end

            %rparNew = [wNew;n;EC50;altNew;drugTypeNew]; % ONLY NEED WHEN CHANGING INPUT MECHANICAL WEIGHT
            rparNew = [wNew;n;EC50;doseNew;drugBindingNew;drugAgonismNew];
            disp(num2str(wNew(3)));
            params = {rparNew,tau,ymax_Knockout,speciesNames,KI,InputCsim,inputNode,tInSim};
            tspan = [0 2329]; 
            options = []; 
            [t2,y2] = ode15s(@dynamicODE,tspan,y0,options,params); 
            y2I = real(interp1(t2,y2,tInSim));
            compiledYData(:, :, i) = y2I(:,:)';
            yEnd1 = y2I(timePointsOfInterest,:);

            colImRNA_ANI = [colImRNA_ANI, y2I(:,101)];
            prolif_ANI = [prolif_ANI, y2I(:,84)];

            sens(:, :, i) = (yEnd1 - yEnd0); %expressing sensitivity as the difference in activity
        end 
        % Calculate change in final area fraction
        CmRNA = sum(y1(:,[101,102]),2);
        peakCol = max(CmRNA);

        if r == 0
            ySim(:) = sum(real(y2I(:,[101,102]))');
            [c1] = MISimODE(ySim, tInSim, peakCol);
            c1End_control = c1(end);
            c1End_knockout = c1(end);
            label = 'Control';
        else
            ySim(:) = sum(real(y2I(:,[101,102]))');
            [c1] = MISimODE(ySim, tInSim, peakCol);
            c1End_knockout = c1(end);
            label = speciesNames{r};
        end

        dataVector = {label, c1End_knockout-c1End_control, max(ySim(:))};
        knockoutData_control(r+1, :) = dataVector;
    end
    
    [sortedNodes, indicies] = sort(knockoutData_control(:,1));
    
    figure
    bar(cell2mat(knockoutData_control(indicies,2)))
    set(gca, 'xtick', 1:108);
    set(gca, 'xticklabel', sortedNodes,'FontSize', 11, 'FontWeight', 'bold');
    xtickangle(90)
    ylabel('Change in collagen area fraction ((drug+KO) - drug)')
    
end

%% Control knockouts

drugsToSimulate = load('drugsToSimulate.mat');
drugsToSimulate = drugsToSimulate.drugsToSimulate; % Extract from struct

formattedReactions = load('formattedReactions.mat');
formattedReactions = formattedReactions.formattedReactions; % Extract from struct

clusteredRowLabels = load('clusteredRowLabels.mat');
clusteredRowLabels = clusteredRowLabels.clusteredRowLabels; % Extract from struct

% Set drug dose or doses (as a vector)
dose_antag = [0.85]; 
dose_ag = [0.85];

[InputCsim,tInSim,inputNode,resNorm,resNormConvert] = InputCurve(0.6, 0.6);

% Parameters and initial conditions
[params,y0] = fib617_params;
tspan = [0 2329]; 
options = [];
timePointsOfInterest = [168,192,336,1176];
% Pull parameters out in order to alter them in the loop without negating the previous alteration
[rpar,tau,ymax,speciesNames,KI]=params{:}; 
w = rpar(1,:);
wNew = rpar(1,:); % ONLY NEED WHEN CHANGING INPUT MECHANICAL WEIGHT
wNew(3) = 0.6; % Altering mechanical input
wNew(23) = 0.6; % Altering the MMP9 & TGFBmRNA => TGFB weight
wNew(24) = 0.6; % Altering the MMP2 & TGFBmRNA => TGFB weight
wNew(25) = 0; % Altering the TGFBmRNA => TGFB weight. 
n = rpar(2,:);
EC50 = rpar(3,:);
dose = rpar(4,:);
drugBinding = rpar(5,:); 
drugAgonism = rpar(6,:); 

%rpar = [wNew;n;EC50;alt;drugType]; % ONLY NEED WHEN CHANGING INPUT MECHANICAL WEIGHT
rpar = [wNew;n;EC50;dose;drugBinding;drugAgonism];

% Default simulation creating a vector showing steady state activation of all species
params = {rpar,tau,ymax,speciesNames,KI,InputCsim,inputNode,tInSim};
[t,y] = ode15s(@dynamicODE,tspan,y0,options,params);


y1 = real(interp1(t,y,tInSim));
yEnd0 = y1(timePointsOfInterest,:);


% Preturbed simulations
sens = zeros(length(timePointsOfInterest),length(speciesNames), height(drugsToSimulate));

drugTargetResponse = zeros(1, height(drugsToSimulate)); 
drugTargetResponseID = zeros(1, height(drugsToSimulate));
colImRNA_ANI = [];
colImRNA_ANI = [colImRNA_ANI, y1(:,101)];
prolif_ANI = [];
prolif_ANI = [prolif_ANI, y1(:,84)];

for r = 0:length(speciesNames)
    disp(num2str(r))
    ymax_Knockout = ymax;
    if r > 0
        ymax_Knockout(r) = 0;
    end
    
    params = {rpar,tau,ymax_Knockout,speciesNames,KI,InputCsim,inputNode,tInSim};
    tspan = [0 2329]; 
    options = []; 
    [t2,y2] = ode15s(@dynamicODE,tspan,y0,options,params); 
    y2I = real(interp1(t2,y2,tInSim));
    
    % Calculate change in final area fraction
    CmRNA = sum(y1(:,[101,102]),2);
    peakCol = max(CmRNA);
    
    if r == 0
        ySim(:) = sum(real(y2I(:,[101,102]))');
        [c1] = MISimODE(ySim, tInSim, peakCol);
        c1End_control = c1(end);
        c1End_knockout = c1(end);
        label = 'Control';
    else
        ySim(:) = sum(real(y2I(:,[101,102]))');
        [c1] = MISimODE(ySim, tInSim, peakCol);
        c1End_knockout = c1(end);
        label = speciesNames{r};
    end
    
    dataVector = {label, c1End_knockout-c1End_control, max(ySim(:))};
    knockoutData_control(r+1, :) = dataVector;
end

figure
[sortedNodes, indicies] = sort(knockoutData_control(:,1));

bar(cell2mat(knockoutData_control(indicies,2)))
set(gca, 'xtick', 1:108);
set(gca, 'xticklabel', sortedNodes,'FontSize', 11, 'FontWeight', 'bold');
xtickangle(90)
ylabel('Change in collagen area fraction ((drug+KO) - drug)')

disp('Done')

    