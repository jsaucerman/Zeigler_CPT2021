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

drugsToSimulate = load('drugsToSimulate.mat');
drugsToSimulate = drugsToSimulate.drugsToSimulate; % Extract from struct

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
sens_col = zeros(height(drugsToSimulate), length(inputLabels));
sens_edafn = zeros(height(drugsToSimulate), length(inputLabels));
sens_prolif = zeros(height(drugsToSimulate), length(inputLabels)); 

sens_mmp1 = zeros(height(drugsToSimulate), length(inputLabels)); 
sens_mmp2 = zeros(height(drugsToSimulate), length(inputLabels)); 
sens_mmp9 = zeros(height(drugsToSimulate), length(inputLabels)); 


drugTargetResponse = [];
drugTargetResponseID = zeros(1, height(drugsToSimulate));

drugsWithNoChangeInCollagen = [];
dataVector = [];
entrestoAgonistResponse = [];

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
            [t2,y2] = ode15s(@fib617_modelODE_drug_static,tspan,y0,options,paramsNew);

            sens_col(i,k) = y2(end,101) - y(end,101); %expressing sensitivity as difference in activity
            sens_prolif(i,k) = y2(end,84) - y(end,84);
            sens_edafn(i,k) = y2(end,86) - y(end,86);
            
            sens_mmp1(i,k) = y2(end,96) - y(end,96);
            sens_mmp2(i,k) = y2(end,97) - y(end,97);
            sens_mmp9(i,k) = y2(end,98) - y(end,98);
        end   
    end
end

% Col I mRNA + Proliferation + EDAFN subplots
subplot(1,3,1)
colImRNAData_clustered = [];
for i = 1:height(drugsToSimulate)
    colImRNAData_clustered = vertcat(colImRNAData_clustered, sens_col(str2double(clusteredRowLabels{i,:}),:));
end
cmaprange = [0.5:0.005:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
caxis([-0.25 0.25]); imagesc(colImRNAData_clustered, [-0.25, 0.25]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(inputLabels));
set(gca,'XTickLabel',inputLabels,'fontsize',11)
set(gca,'YTick',1:height(drugsToSimulate));
drugsToSimulate{21, 7} = {"MMP1;2;9;14"};
[padding,list_numDrugs, list_drugAction, list_drugType, list_drugTarget, rowLabels_drugs] = createRowLabels(drugsToSimulate);
rowLabels_drugs_clustered = [];
for i = 1:height(drugsToSimulate)
    rowLabels_drugs_clustered = vertcat(rowLabels_drugs_clustered, rowLabels_drugs(str2double(clusteredRowLabels{i,:}),:));
end
set(gca,'YTickLabel',rowLabels_drugs_clustered(:,1),'fontsize',11)
set(gca, 'fontname', 'FixedWidth');
xtickangle(45)
title('Effect of Drug on Col I mRNA', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Drugs', 'FontSize', 11, 'FontWeight', 'bold')

subplot(1,3,2)
prolifData_clustered = [];
for i = 1:height(drugsToSimulate)
    prolifData_clustered = vertcat(prolifData_clustered, sens_prolif(str2double(clusteredRowLabels{i,:}),:));
end
caxis([-0.25 0.25]); imagesc(prolifData_clustered, [-0.25, 0.25]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(inputLabels));
set(gca,'XTickLabel',inputLabels,'fontsize',11)
set(gca,'YTickLabel','')
xtickangle(45)
title('Effect of Drug on Proliferation', 'FontSize', 14, 'FontWeight', 'bold')
xlabel('Days Post-Myocardial Infarction', 'FontSize', 11, 'FontWeight', 'bold')
set(gca, 'fontname', 'FixedWidth');

subplot(1,3,3)
edafnData_clustered = [];
for i = 1:height(drugsToSimulate)
    edafnData_clustered = vertcat(edafnData_clustered, sens_edafn(str2double(clusteredRowLabels{i,:}),:));
end
cmaprange = [0.5:0.005:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
caxis([-0.25 0.25]); imagesc(edafnData_clustered, [-0.25, 0.25]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(inputLabels));
set(gca,'XTickLabel',inputLabels,'fontsize',11)
title('Effect of Drug on EDAFN', 'FontSize', 14, 'FontWeight', 'bold')
c = colorbar; 
ylabel(c, {'Change in Activity', '(drugged - non drugged)'}, 'Rotation', 270, 'FontWeight', 'bold', 'FontSize', 14);
pos = c.Label.Position; 
c.Label.Position = [pos(1)+1, pos(2), pos(3)];
xlabel('Days Post-Myocardial Infarction', 'FontSize', 11, 'FontWeight', 'bold')
xtickangle(45)
set(gca, 'fontname', 'FixedWidth');
set(gca, 'YTickLabel', '');


% MMP1 + MMP2 + MMP9 subplots
subplot(1,3,1);
mmp1Data_clustered = [];
for i = 1:height(drugsToSimulate)
    mmp1Data_clustered = vertcat(mmp1Data_clustered, sens_mmp1(str2double(clusteredRowLabels{i,:}),:));
end
cmaprange = [0.5:0.005:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
caxis([-0.25 0.25]); imagesc(mmp1Data_clustered, [-0.25, 0.25]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(inputLabels));
set(gca,'XTickLabel',inputLabels,'fontsize',11)
title('Effect of Drug on MMP1', 'FontSize', 14, 'FontWeight', 'bold')
xlabel('Days Post-Myocardial Infarction', 'FontSize', 11, 'FontWeight', 'bold')
xtickangle(45)
set(gca, 'fontname', 'FixedWidth');
set(gca,'YTick',1:height(drugsToSimulate));
set(gca,'YTickLabel',rowLabels_drugs_clustered(:,1),'fontsize',11)
title('Effect of Drug on MMP1', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Drugs', 'FontSize', 11, 'FontWeight', 'bold')

subplot(1,3,2)
mmp2Data_clustered = [];
for i = 1:height(drugsToSimulate)
    mmp2Data_clustered = vertcat(mmp2Data_clustered, sens_mmp2(str2double(clusteredRowLabels{i,:}),:));
end
cmaprange = [0.5:0.005:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
caxis([-0.25 0.25]); imagesc(mmp2Data_clustered, [-0.25, 0.25]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(inputLabels));
set(gca,'XTickLabel',inputLabels,'fontsize',11)
title('Effect of Drug on MMP2', 'FontSize', 14, 'FontWeight', 'bold')
xlabel('Days Post-Myocardial Infarction', 'FontSize', 11, 'FontWeight', 'bold')
xtickangle(45)
set(gca, 'YTickLabel', '');
set(gca, 'fontname', 'FixedWidth');

subplot(1,3,3)
mmp9Data_clustered = [];
for i = 1:height(drugsToSimulate)
    mmp9Data_clustered = vertcat(mmp9Data_clustered, sens_mmp9(str2double(clusteredRowLabels{i,:}),:));
end
cmaprange = [0.5:0.005:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
caxis([-0.25 0.25]); imagesc(mmp9Data_clustered, [-0.25, 0.25]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(inputLabels));
set(gca,'XTickLabel',inputLabels,'fontsize',11)
title('Effect of Drug on MMP9', 'FontSize', 14, 'FontWeight', 'bold')
c = colorbar; 
ylabel(c, {'Change in Activity', '(drugged - non drugged)'}, 'Rotation', 270, 'FontWeight', 'bold', 'FontSize', 14);
pos = c.Label.Position; 
c.Label.Position = [pos(1)+1, pos(2), pos(3)];
xlabel('Days Post-Myocardial Infarction', 'FontSize', 11, 'FontWeight', 'bold')
xtickangle(45)
set(gca, 'YTickLabel', '');
set(gca, 'fontname', 'FixedWidth');
