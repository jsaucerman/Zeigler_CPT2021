% Written by Anirudha Chandrabhatla 
% Last Updated 11/24/2020
% Version 4.0

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

daysOfInterest = [3, 4]; % Days 7 and 42
speciesToKnockout = [find(strcmpi(speciesNames, "TGFBmRNA")), find(strcmpi(speciesNames, "MMP2"))];
drugOfInterest = find(strcmpi(drugsToSimulate.Drug, "Arsenic Trioxide"));
additionalDrug = find(strcmpi(drugsToSimulate.Drug, "Marimastat"));

barGraphData = zeros(2, 4+(length(speciesToKnockout)+1)*4);

c1 = [];
c2 = [];
c3 = [];

barGraphData(1,1) = yEnd0(daysOfInterest(1), 101); % Control, day 7
barGraphData(1,length(barGraphData)/2 + 1) = yEnd0(daysOfInterest(2), 101); % Control, day 42

CmRNA = sum(y1(:,[101,102]),2);
peakCol = max(CmRNA);
[c0, days, time_rat, AF_rat, AFsd_rat] = MISimODE(CmRNA, tInSim, peakCol); % Control

% Control knockout
[rpar,tau,ymax,speciesNames,KI]=params{:};

for x = 1:length(speciesToKnockout)+1
    ymaxNew = ymax;
    if x == length(speciesToKnockout)+1
        doseNew = dose;
            drugBindingNew = drugBinding;
            drugAgonismNew = drugAgonism;

            % The drug is an agonist
            if strcmp(drugsToSimulate.IsAgonist{additionalDrug}, 'Yes') == 1
                if isempty(find(drugsToSimulate.AgonistTarget{additionalDrug} == ';')) % Drug has one agonist target
                    locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, drugsToSimulate.AgonistTargetGeneID{additionalDrug});
                    if strcmpi('Competitive', drugsToSimulate.DrugAction{additionalDrug}) % Competitive, agonist
                        drugBindingNew(locationOfReactions) = 1;
                        drugAgonismNew(locationOfReactions) = 1;
                        doseNew(locationOfReactions) = -1*dose_ag(1);
                    else % Non-Competitive, agonist
                        drugBindingNew(locationOfReactions) = -1;
                        drugAgonismNew(locationOfReactions) = 1;
                        doseNew(locationOfReactions) = dose_ag(1);
                    end
                else % Drug has multiple targets
                    geneIDsOfTargets = strsplit(drugsToSimulate.AgonistTargetGeneID{additionalDrug}, ';');
                    for k = 1:length(geneIDsOfTargets)
                        locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, geneIDsOfTargets{k});
                        if strcmpi('Competitive', drugsToSimulate.DrugAction{additionalDrug})
                            drugBindingNew(locationOfReactions) = 1;
                            drugAgonismNew(locationOfReactions) = 1;
                            doseNew(locationOfReactions) = -1*dose_ag(1);
                        else 
                            drugBindingNew(locationOfReactions) = -1;
                            drugAgonismNew(locationOfReactions) = 1;
                            doseNew(locationOfReactions) = dose_ag(1);
                        end
                    end
                end
            end

            % The drug is an antagonist
            if strcmp(drugsToSimulate.IsAntagonist{additionalDrug}, 'Yes') == 1
                if isempty(find(drugsToSimulate.AntagonistTarget{additionalDrug} == ';')) % Drug has one antagonist target
                    locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, drugsToSimulate.AntagonistTargetGeneID{additionalDrug});
                    if strcmpi('Competitive', drugsToSimulate.DrugAction{additionalDrug}) % Competitive, antagonist
                        drugBindingNew(locationOfReactions) = 1;
                        drugAgonismNew(locationOfReactions) = -1;
                        doseNew(locationOfReactions) = dose_antag;
                    else
                        drugBindingNew(locationOfReactions) = -1;
                        drugAgonismNew(locationOfReactions) = -1;
                        doseNew(locationOfReactions) = dose_antag;
                    end
                else
                    geneIDsOfTargets = strsplit(drugsToSimulate.AntagonistTargetGeneID{additionalDrug}, ';');
                    for p = 1:length(geneIDsOfTargets)
                        locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, geneIDsOfTargets{p});
                        if strcmpi('Competitive', drugsToSimulate.DrugAction{additionalDrug})
                            drugBindingNew(locationOfReactions) = 1;
                            drugAgonismNew(locationOfReactions) = -1;
                            doseNew(locationOfReactions) = dose_antag;
                        else 
                            drugBindingNew(locationOfReactions) = -1;
                            drugAgonismNew(locationOfReactions) = -1;
                            doseNew(locationOfReactions) = dose_antag;
                        end
                    end
                end
            end

            rparNew = [wNew;n;EC50;doseNew;drugBindingNew;drugAgonismNew];
            params = {rparNew,tau,ymaxNew,speciesNames,KI,InputCsim,inputNode,tInSim};
            tspan = [0 2329]; 
            options = []; 
            [t,y] = ode15s(@dynamicODE,tspan,y0,options,params); 
    else
        ymaxNew(speciesToKnockout(x)) = 0;
        params = {rpar,tau,ymaxNew,speciesNames,KI,InputCsim,inputNode,tInSim};
        [t,y] = ode15s(@dynamicODE,tspan,y0,options,params);
    end
    
    y1_ControlKnockout = real(interp1(t,y,tInSim));
    
    barGraphData(1,x+1) = y1_ControlKnockout(timePointsOfInterest(daysOfInterest(1)), 101); % Control + KO
    
    barGraphData(1,x+length(barGraphData)/2 + 1) = y1_ControlKnockout(timePointsOfInterest(daysOfInterest(2)), 101);
    
    ySim(:) = sum(real(y1_ControlKnockout(:,[101,102])),2); % Control + KO
    c1 = horzcat(c1, MISimODE(ySim, tInSim, peakCol));
end


% Drug + drug / knockout
for j = 1:length(speciesToKnockout)+1
    for i = drugOfInterest:drugOfInterest
    
        for cycle = 1:2
            ymaxNew = ymax;
            
            doseNew = dose;
            drugBindingNew = drugBinding;
            drugAgonismNew = drugAgonism;
            
            %KI_new = KI;
            if cycle == 2
                if j == length(speciesToKnockout)+1
                    % The drug is an agonist
                    if strcmp(drugsToSimulate.IsAgonist{additionalDrug}, 'Yes') == 1
                        if isempty(find(drugsToSimulate.AgonistTarget{additionalDrug} == ';')) % Drug has one agonist target
                            locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, drugsToSimulate.AgonistTargetGeneID{additionalDrug});
                            if strcmpi('Competitive', drugsToSimulate.DrugAction{additionalDrug}) % Competitive, agonist
                                drugBindingNew(locationOfReactions) = 1;
                                drugAgonismNew(locationOfReactions) = 1;
                                doseNew(locationOfReactions) = -1*dose_ag(1);
                            else % Non-Competitive, agonist
                                drugBindingNew(locationOfReactions) = -1;
                                drugAgonismNew(locationOfReactions) = 1;
                                doseNew(locationOfReactions) = dose_ag(1);
                            end
                        else % Drug has multiple targets
                            geneIDsOfTargets = strsplit(drugsToSimulate.AgonistTargetGeneID{additionalDrug}, ';');
                            for k = 1:length(geneIDsOfTargets)
                                locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, geneIDsOfTargets{k});
                                if strcmpi('Competitive', drugsToSimulate.DrugAction{additionalDrug})
                                    drugBindingNew(locationOfReactions) = 1;
                                    drugAgonismNew(locationOfReactions) = 1;
                                    doseNew(locationOfReactions) = -1*dose_ag(1);
                                else 
                                    drugBindingNew(locationOfReactions) = -1;
                                    drugAgonismNew(locationOfReactions) = 1;
                                    doseNew(locationOfReactions) = dose_ag(1);
                                end
                            end
                        end
                    end

                    % The drug is an antagonist
                    if strcmp(drugsToSimulate.IsAntagonist{additionalDrug}, 'Yes') == 1
                        if isempty(find(drugsToSimulate.AntagonistTarget{additionalDrug} == ';')) % Drug has one antagonist target
                            locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, drugsToSimulate.AntagonistTargetGeneID{additionalDrug});
                            if strcmpi('Competitive', drugsToSimulate.DrugAction{additionalDrug}) % Competitive, antagonist
                                drugBindingNew(locationOfReactions) = 1;
                                drugAgonismNew(locationOfReactions) = -1;
                                doseNew(locationOfReactions) = dose_antag;
                            else
                                drugBindingNew(locationOfReactions) = -1;
                                drugAgonismNew(locationOfReactions) = -1;
                                doseNew(locationOfReactions) = dose_antag;
                            end
                        else
                            geneIDsOfTargets = strsplit(drugsToSimulate.AntagonistTargetGeneID{additionalDrug}, ';');
                            for p = 1:length(geneIDsOfTargets)
                                locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, geneIDsOfTargets{p});
                                if strcmpi('Competitive', drugsToSimulate.DrugAction{additionalDrug})
                                    drugBindingNew(locationOfReactions) = 1;
                                    drugAgonismNew(locationOfReactions) = -1;
                                    doseNew(locationOfReactions) = dose_antag;
                                else 
                                    drugBindingNew(locationOfReactions) = -1;
                                    drugAgonismNew(locationOfReactions) = -1;
                                    doseNew(locationOfReactions) = dose_antag;
                                end
                            end
                        end
                    end
                else
                    ymaxNew(speciesToKnockout(j)) = 0;
                end
            end

            % The drug is an agonist
            if strcmp(drugsToSimulate.IsAgonist{i}, 'Yes') == 1
                if isempty(find(drugsToSimulate.AgonistTarget{i} == ';')) % Drug has one agonist target
                    locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, drugsToSimulate.AgonistTargetGeneID{i});
                    if strcmpi('Competitive', drugsToSimulate.DrugAction{i}) % Competitive, agonist
                        drugBindingNew(locationOfReactions) = 1;
                        drugAgonismNew(locationOfReactions) = 1;
                        doseNew(locationOfReactions) = -1*dose_ag(1);
                    else % Non-Competitive, agonist
                        drugBindingNew(locationOfReactions) = -1;
                        drugAgonismNew(locationOfReactions) = 1;
                        doseNew(locationOfReactions) = dose_ag(1);
                    end
                else % Drug has multiple targets
                    geneIDsOfTargets = strsplit(drugsToSimulate.AgonistTargetGeneID{i}, ';');
                    for k = 1:length(geneIDsOfTargets)
                        locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, geneIDsOfTargets{k});
                        if strcmpi('Competitive', drugsToSimulate.DrugAction{i})
                            drugBindingNew(locationOfReactions) = 1;
                            drugAgonismNew(locationOfReactions) = 1;
                            doseNew(locationOfReactions) = -1*dose_ag(1);
                        else 
                            drugBindingNew(locationOfReactions) = -1;
                            drugAgonismNew(locationOfReactions) = 1;
                            doseNew(locationOfReactions) = dose_ag(1);
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
                        doseNew(locationOfReactions) = dose_antag;
                    else
                        drugBindingNew(locationOfReactions) = -1;
                        drugAgonismNew(locationOfReactions) = -1;
                        doseNew(locationOfReactions) = dose_antag;
                    end
                else
                    geneIDsOfTargets = strsplit(drugsToSimulate.AntagonistTargetGeneID{i}, ';');
                    for p = 1:length(geneIDsOfTargets)
                        locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, geneIDsOfTargets{p});
                        if strcmpi('Competitive', drugsToSimulate.DrugAction{i})
                            drugBindingNew(locationOfReactions) = 1;
                            drugAgonismNew(locationOfReactions) = -1;
                            doseNew(locationOfReactions) = dose_antag;
                        else 
                            drugBindingNew(locationOfReactions) = -1;
                            drugAgonismNew(locationOfReactions) = -1;
                            doseNew(locationOfReactions) = dose_antag;
                        end
                    end
                end
            end

            rparNew = [wNew;n;EC50;doseNew;drugBindingNew;drugAgonismNew];
            disp(num2str(wNew(3)));
            params = {rparNew,tau,ymaxNew,speciesNames,KI,InputCsim,inputNode,tInSim};
            tspan = [0 2329]; 
            options = []; 
            [t2,y2] = ode15s(@dynamicODE,tspan,y0,options,params); 
            y2I = real(interp1(t2,y2,tInSim));
            yEnd1 = y2I(timePointsOfInterest,:);
            
            if cycle == 1
                barGraphData(1,length(speciesToKnockout)+3) = yEnd1(daysOfInterest(1), 101); % Drug
                barGraphData(1,3*(length(speciesToKnockout)+2)+1) = yEnd1(daysOfInterest(2), 101);
                ySim(:) = sum(real(y2I(:,[101, 102])),2); % Drug
                c2 = horzcat(c2, MISimODE(ySim, tInSim, peakCol));
            else
                barGraphData(1,j+length(speciesToKnockout)+3) = yEnd1(daysOfInterest(1), 101); % Drug + KO
                barGraphData(1,j+3*(length(speciesToKnockout)+2)+1) = yEnd1(daysOfInterest(2), 101);
                ySim(:) = sum(real(y2I(:,[101, 102])),2); % Drug + KO
                c3 = horzcat(c3, MISimODE(ySim, tInSim, peakCol));
            end
            
        end
    end
end


barGraphData(2,:) = [c0(132), c1(132,:), c2(132), c3(132,:), c0(460), c1(460,:), c2(460), c3(460,:)]; % Area fraction at day 7 and 42

barGraphData = [barGraphData(:,1:length(barGraphData)/2), zeros(2,1), barGraphData(:,length(barGraphData)/2 + 1:end)];
speciesToKnockout = [find(strcmpi(speciesNames, "TGFBmRNA")), find(strcmpi(speciesNames, "MMP2"))];

% Collagen
figure
bar(barGraphData(1,:));
set(gca, 'xtick', 1:length(barGraphData));
set(gca, 'xticklabel', {'Control', 'Control + TGFBmRNA KO', 'Control + MMP2 KO', 'Control + Marimastat', 'ATO', 'ATO + TGFBmRNA KO', 'ATO + MMP2 KO', 'ATO + Marimastat', '', 'Control', 'Control + TGFBmRNA KO', 'Control + MMP2 KO', 'Control + Marimastat', 'ATO', 'ATO + TGFBmRNA KO', 'ATO + MMP2 KO', 'ATO + Marimastat'},'FontSize', 11, 'FontWeight', 'bold');
xtickangle(45)
ylim([0 0.5])
box off
ylabel('Normalized Col I mRNA expression', 'FontSize', 11, 'FontWeight', 'bold')

% Area Fraction
figure
bar(barGraphData(2,:))
set(gca, 'xtick', 1:1:length(barGraphData));
set(gca, 'xticklabel', {'Control', 'Control + TGFBmRNA KO', 'Control + MMP2 KO', 'Control + Marimastat', 'ATO', 'ATO + TGFBmRNA KO', 'ATO + MMP2 KO', 'ATO + Marimastat', '', 'Control', 'Control + TGFBmRNA KO', 'Control + MMP2 KO', 'Control + Marimastat', 'ATO', 'ATO + TGFBmRNA KO', 'ATO + MMP2 KO', 'ATO + Marimastat'},'FontSize', 11, 'FontWeight', 'bold');
xtickangle(45)
ylim([0 100])
box off
ylabel('Collagen area fraction', 'FontSize', 11, 'FontWeight', 'bold')
