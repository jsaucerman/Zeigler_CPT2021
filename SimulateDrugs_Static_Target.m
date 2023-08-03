% Written by Anirudha Chandrabhatla 
% Last Updated 11/24/2020
% Version 9.0

% %% If any inputs have changed, run the following section: 
% 
%   warning off;
% % Species information from fib617_references.xlsx, 'species' tab 
% nodeGenesTable = readtable('fib617_references.xlsx', 'Sheet', 'species');
% 
% % Reaction information from fib417_references.xlsx, 'reactions' tab 
% networkReactions = readtable('fib617_references.xlsx', 'Sheet', 'reactions');
% 
% % Drug information from AllPharmActive.xlsx
% approvedTargetsTable = readtable('AllPharmActive.xlsx');
% 
% % Drug ID to drug name matching from drugIDsAndNamesFull516.xlsx
% drugIDToNameConversion = readtable('drugIDsAndNamesFull516.xlsx');
% 
% % Manually curated competitive/non-competitive classifications for all drugs that target a node in the network
% compNonCompClassification = readtable('All Drugs Comp NonComp Classification.xlsx', 'Sheet', 'All Drugs Condensed');
% 
% load('finalDrugOutputNetworkTargets.mat'); % Result from the webscraper
% drugInformation = finalDrugOutputNetworkTargets;
% 
% geneNames = nodeGenesTable.('geneName');  % Gene names from the model
% geneIDs = nodeGenesTable.('ID');          % Gene IDs from the model
% 
% % Drug IDs from the database of approved targets
% drugIDs_fromTargetList = approvedTargetsTable.('DrugIDs'); 
% % Gene names of the drug targets. Retrieved from the database of approved targets.
% drugGeneTargets = approvedTargetsTable.('GeneName'); 
% 
% % Drug names from the database of all drugs
% drugNames = drugIDToNameConversion.('Name');
% % Drug IDs from the database of all drugs 
% drugIDs_fromTargetList_forNameConversion = drugIDToNameConversion.('DrugBankID');
% 
% [drugsToSimulate, formattedReactions] = formatInputs(geneNames, geneIDs, drugGeneTargets, drugIDs_fromTargetList, drugIDs_fromTargetList_forNameConversion, drugNames, networkReactions, compNonCompClassification, drugInformation);

%% Inputs for simulations


drugsToSimulate = load('drugsToSimulate_doseResponse.mat');
drugsToSimulate = drugsToSimulate.drugsToSimulate;

formattedReactions = load('formattedReactions.mat');
formattedReactions = formattedReactions.formattedReactions; % Extract from struct

clusteredRowLabels = load('clusteredRowLabels.mat');
clusteredRowLabels = clusteredRowLabels.clusteredRowLabels; % Extract from struct

% Set drug dose or doses (as a vector)
doseResponse = [0.85];


%% Part 2: Generate graph of results

inputNodeW = {5,2,1,[5,10],[2,9],[7,1]};
inputLabels = {'Control','IL1','TGFB','AngII','IL1 + NP', 'TGFB + ET1', 'AngII + NE'};
deltaIn = 0.6;

% Preturbed simulations
sens = zeros(height(drugsToSimulate), length(inputLabels));


drugTargetResponse = [];
drugTargetResponseID = zeros(1, height(drugsToSimulate));

drugsWithNoChangeInCollagen = [];
dataVector = [];
entrestoAgonistResponse = [];

antagIndicies = drugsToSimulate.AntagonistTargetIndex;
agIndicies = drugsToSimulate.AgonistTargetIndex;
antagIndicies(strcmp("", antagIndicies)) = agIndicies(strcmp("", antagIndicies));
targetIndicies = double(antagIndicies);

for k = 1:length(inputLabels)
    % Parameters and initial conditions
    [params,y0] = fib617_params;
    tspan = [0 500]; 
    options = [];

    % Pull parameters out in order to alter them in the loop without negating the previous alteration
    [rpar,tau,ymax,speciesNames,KI]=params{:}; 
    w = rpar(1,:);
    wNew = rpar(1,:);
    wNew(1:11) = 0.25;
    if k > 1 % First run is a control simulation. 
        wNew(cell2mat(inputNodeW(k-1))) = deltaIn;
    end
    wNew(3) = 0.85; % Setting mechanical input
    wNew(23) = 0.6;
    wNew(24) = 0.6;
    wNew(25) = 0; % Altering the TGFBmRNA --> TGFB reaction weight
    n = rpar(2,:);
    EC50 = rpar(3,:);
    dose = rpar(4,:);
    drugType = rpar(5,:);
    drugBinding = rpar(5,:); 
    drugAgonism = rpar(6,:); 

    rpar = [wNew;n;EC50;dose;drugBinding;drugAgonism];

    % Default simulation creating a vector showing steady state activation of all species
    %params = {rpar,tau,ymaxNew,speciesNames,KI}; % Only need if doing a knockout
    params = {rpar,tau,ymax,speciesNames,KI};
    [t,y] = ode15s(@fib617_modelODE_drug_static,tspan,y0,options,params);
    ycol0(k) = y(end,101);
    
    dataVector = [];
    % Iterate through 'drugsToSimulate' and simulate the addition of each drug into the network.
    for i = 1:height(drugsToSimulate)
        if i == 32
            disp('here')
        end
        disp(drugsToSimulate{i, 1}); % Prints drug name to the command window  

        for j = 1:length(doseResponse)

            doseNew = dose;
            drugBindingNew = drugBinding;
            drugAgonismNew = drugAgonism;

            % The drug is an agonist
            if strcmp(drugsToSimulate.IsAgonist{i}, 'Yes') == 1
                if isempty(find(drugsToSimulate.AgonistTarget{i} == ';')) % Drug has one agonist target
                    locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, drugsToSimulate.AgonistTargetGeneID{i});
                    if strcmp('Competitive', drugsToSimulate.DrugAction{i}) % Competitive, agonist
                        drugBindingNew(locationOfReactions) = 1;
                        drugAgonismNew(locationOfReactions) = 1;
                        doseNew(locationOfReactions) = -1*doseResponse(j);
                    else % Non-Competitive, agonist
                        drugBindingNew(locationOfReactions) = -1;
                        drugAgonismNew(locationOfReactions) = 1;
                        doseNew(locationOfReactions) = doseResponse(j);
                    end
                else % Drug has multiple targets
                    geneIDsOfTargets = strsplit(drugsToSimulate.AgonistTargetGeneID{i}, ';');
                    for m = 1:length(geneIDsOfTargets)
                        locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, geneIDsOfTargets{m});
                        if strcmp('Competitive', drugsToSimulate.DrugAction{i})
                            drugBindingNew(locationOfReactions) = 1;
                            drugAgonismNew(locationOfReactions) = 1;
                            doseNew(locationOfReactions) = -1*doseResponse(j);
                        else 
                            drugBindingNew(locationOfReactions) = -1;
                            drugAgonismNew(locationOfReactions) = 1;
                            doseNew(locationOfReactions) = doseResponse(j);
                        end
                    end
                end
            end

            % The drug is an antagonist
            if strcmp(drugsToSimulate.IsAntagonist{i}, 'Yes') == 1
                if isempty(find(drugsToSimulate.AntagonistTarget{i} == ';')) % Drug has one antagonist target
                    locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, drugsToSimulate.AntagonistTargetGeneID{i});
                    if strcmp('Competitive', drugsToSimulate.DrugAction{i}) % Competitive, antagonist
                        drugBindingNew(locationOfReactions) = 1;
                        drugAgonismNew(locationOfReactions) = -1;
                        doseNew(locationOfReactions) = doseResponse(j);
                    else
                        drugBindingNew(locationOfReactions) = -1;
                        drugAgonismNew(locationOfReactions) = -1;
                        doseNew(locationOfReactions) = doseResponse(j);
                    end
                else
                    geneIDsOfTargets = strsplit(drugsToSimulate.AntagonistTargetGeneID{i}, ';');
                    for p = 1:length(geneIDsOfTargets)
                        locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, geneIDsOfTargets{p});
                        if strcmp('Competitive', drugsToSimulate.DrugAction{i})
                            drugBindingNew(locationOfReactions) = 1;
                            drugAgonismNew(locationOfReactions) = -1;
                            doseNew(locationOfReactions) = doseResponse(j);
                        else 
                            drugBindingNew(locationOfReactions) = -1;
                            drugAgonismNew(locationOfReactions) = -1;
                            doseNew(locationOfReactions) = doseResponse(j);
                        end
                    end
                end
            end
            disp(num2str(wNew(3)));
            rparNew = [wNew;n;EC50;doseNew;drugBindingNew;drugAgonismNew];
            %paramsNew = {rparNew,tau,ymaxNew,speciesNames,KI}; % Only need if doing a knockout
            paramsNew = {rparNew,tau,ymax,speciesNames,KI};
            [t2,y2] = ode15s(@fib617_modelODE_drug_static,tspan,y0,options,paramsNew); % Make sure y0 is correct here.
                    
            sens(i,k) = y2(end,targetIndicies(i)) - y(end,targetIndicies(i)); % expressing sensitivity as difference in activity

        end   
    end
end


figure;
targetData_clustered = [];
for i = 1:height(drugsToSimulate)
    targetData_clustered = vertcat(targetData_clustered, sens(clusteredRowLabels(i,:),:));
end
cmaprange = [0.5:0.005:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
caxis([-0.25 0.25]); imagesc(targetData_clustered, [-1, 1]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(inputLabels));
set(gca,'XTickLabel',inputLabels,'fontsize',11)
set(gca,'YTick',1:height(drugsToSimulate));
drugsToSimulate{21, 7} = {"MMP1;2;9;14"};
[padding,list_numDrugs, list_drugAction, list_drugType, list_drugTarget, rowLabels_drugs] = createRowLabels(drugsToSimulate);
rowLabels_drugs_clustered = [];
for i = 1:height(drugsToSimulate)
    rowLabels_drugs_clustered = vertcat(rowLabels_drugs_clustered, rowLabels_drugs(clusteredRowLabels(i,:),:));
end
set(gca,'YTickLabel',rowLabels_drugs_clustered(:,1),'fontsize',11)
set(gca, 'fontname', 'FixedWidth');
xtickangle(45)
title('Effect of drug on target', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Drugs', 'FontSize', 11, 'FontWeight', 'bold')

disp('here')

