% Written by Anirudha Chandrabhatla 
% Last Updated 11/24/2020
% Version 9.0

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
% %load('finalDrugOutputNetworkTargets.mat'); % Result from the webscraper
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
wNew = rpar(1,:);
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

% Iterate through 'drugsToSimulate' and simulate the addition of each drug into the network.
for i = 1:height(drugsToSimulate)
    ymaxNew = ymax;
    disp(drugsToSimulate{i, 1}); % Prints drug name to the command window 
    for j = 1:length(dose_ag)
        
        doseNew = dose;
        drugBindingNew = drugBinding;
        drugAgonismNew = drugAgonism;

        % Drug is an agonist
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

        % Drug is an antagonist
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
        params = {rparNew,tau,ymaxNew,speciesNames,KI,InputCsim,inputNode,tInSim};
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
end

% Col I mRNA + Proliferation Subplots
time = {'0 days','1 day','7 days','42 days'};
subplot(1,3,1)
colImRNAData = reshape(sens(:,101,:), [length(timePointsOfInterest) height(drugsToSimulate)])'; 
colImRNAData_clustered = [];
for i = 1:height(drugsToSimulate)
    colImRNAData_clustered = vertcat(colImRNAData_clustered, colImRNAData(str2double(clusteredRowLabels{i,:}),:));
end
cmaprange = [0.5:0.005:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
caxis([-0.25 0.25]); imagesc(colImRNAData_clustered, [-0.25, 0.25]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(timePointsOfInterest));
set(gca,'XTickLabel',time,'fontsize',11)
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
prolifData = reshape(sens(:,84,:), [length(timePointsOfInterest) height(drugsToSimulate)])'; 
prolifData_clustered = [];
for i = 1:height(drugsToSimulate)
    prolifData_clustered = vertcat(prolifData_clustered, prolifData(str2double(clusteredRowLabels{i,:}),:));
end
caxis([-0.25 0.25]); imagesc(prolifData_clustered, [-0.25, 0.25]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(timePointsOfInterest));
set(gca,'XTickLabel',time,'fontsize',11)
set(gca,'YTickLabel','')
xtickangle(45)
title('Effect of Drug on Proliferation', 'FontSize', 14, 'FontWeight', 'bold')
xlabel('Days Post-Myocardial Infarction', 'FontSize', 11, 'FontWeight', 'bold')
set(gca, 'fontname', 'FixedWidth');

subplot(1,3,3)
edafnData = reshape(sens(:,86,:), [length(timePointsOfInterest) height(drugsToSimulate)])'; 
edafnData_clustered = [];
for i = 1:height(drugsToSimulate)
    edafnData_clustered = vertcat(edafnData_clustered, edafnData(str2double(clusteredRowLabels{i,:}),:));
end
cmaprange = [0.5:0.005:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
caxis([-0.25 0.25]); imagesc(edafnData_clustered, [-0.25, 0.25]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(timePointsOfInterest));
set(gca,'XTickLabel',time,'fontsize',11)
title('Effect of Drug on EDAFN', 'FontSize', 14, 'FontWeight', 'bold')
c = colorbar; 
ylabel(c, {'Change in Activity', '(drugged - non drugged)'}, 'Rotation', 270, 'FontWeight', 'bold', 'FontSize', 14);
pos = c.Label.Position; 
c.Label.Position = [pos(1)+1, pos(2), pos(3)];
xlabel('Days Post-Myocardial Infarction', 'FontSize', 11, 'FontWeight', 'bold')
xtickangle(45)
set(gca, 'fontname', 'FixedWidth');
set(gca, 'YTickLabel', '');

% MMP1
subplot(1,3,1);
mmp1Data = reshape(sens(:,96,:), [length(timePointsOfInterest) height(drugsToSimulate)])'; 
mmp1Data_clustered = [];
for i = 1:height(drugsToSimulate)
    mmp1Data_clustered = vertcat(mmp1Data_clustered, mmp1Data(str2double(clusteredRowLabels{i,:}),:));
end
cmaprange = [0.5:0.005:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
caxis([-0.25 0.25]); imagesc(mmp1Data_clustered, [-0.25, 0.25]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(timePointsOfInterest));
set(gca,'XTickLabel',time,'fontsize',11)
title('Effect of Drug on MMP1', 'FontSize', 14, 'FontWeight', 'bold')
xlabel('Days Post-Myocardial Infarction', 'FontSize', 11, 'FontWeight', 'bold')
xtickangle(45)
set(gca, 'fontname', 'FixedWidth');
set(gca,'YTick',1:height(drugsToSimulate));
set(gca,'YTickLabel',rowLabels_drugs_clustered(:,1),'fontsize',11)
title('Effect of Drug on MMP1', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Drugs', 'FontSize', 11, 'FontWeight', 'bold')

subplot(1,3,2)
mmp2Data = reshape(sens(:,97,:), [length(timePointsOfInterest) height(drugsToSimulate)])'; 
mmp2Data_clustered = [];
for i = 1:height(drugsToSimulate)
    mmp2Data_clustered = vertcat(mmp2Data_clustered, mmp2Data(str2double(clusteredRowLabels{i,:}),:));
end
cmaprange = [0.5:0.005:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
caxis([-0.25 0.25]); imagesc(mmp2Data_clustered, [-0.25, 0.25]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(timePointsOfInterest));
set(gca,'XTickLabel',time,'fontsize',11)
title('Effect of Drug on MMP2', 'FontSize', 14, 'FontWeight', 'bold')
xlabel('Days Post-Myocardial Infarction', 'FontSize', 11, 'FontWeight', 'bold')
xtickangle(45)
set(gca, 'YTickLabel', '');
set(gca, 'fontname', 'FixedWidth');

subplot(1,3,3)
mmp9Data = reshape(sens(:,98,:), [length(timePointsOfInterest) height(drugsToSimulate)])'; 
mmp9Data_clustered = [];
for i = 1:height(drugsToSimulate)
    mmp9Data_clustered = vertcat(mmp9Data_clustered, mmp9Data(str2double(clusteredRowLabels{i,:}),:));
end
cmaprange = [0.5:0.005:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
caxis([-0.25 0.25]); imagesc(mmp9Data_clustered, [-0.25, 0.25]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(timePointsOfInterest));
set(gca,'XTickLabel',time,'fontsize',11)
title('Effect of Drug on MMP9', 'FontSize', 14, 'FontWeight', 'bold')
c = colorbar; 
ylabel(c, {'Change in Activity', '(drugged - non drugged)'}, 'Rotation', 270, 'FontWeight', 'bold', 'FontSize', 14);
pos = c.Label.Position; 
c.Label.Position = [pos(1)+1, pos(2), pos(3)];
xlabel('Days Post-Myocardial Infarction', 'FontSize', 11, 'FontWeight', 'bold')
xtickangle(45)
set(gca, 'YTickLabel', '');
set(gca, 'fontname', 'FixedWidth');

figure;
mmp14Data = reshape(sens(:,99,:), [length(timePointsOfInterest) height(drugsToSimulate)])'; 
mmp14Data_clustered = [];
for i = 1:height(drugsToSimulate)
    mmp14Data_clustered = vertcat(mmp14Data_clustered, mmp14Data(str2double(clusteredRowLabels{i,:}),:));
end
cmaprange = [0.5:0.005:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
caxis([-0.25 0.25]); imagesc(mmp14Data_clustered, [-0.25, 0.25]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(timePointsOfInterest));
set(gca,'XTickLabel',time,'fontsize',11)
title('Effect of Drug on MMP14', 'FontSize', 14, 'FontWeight', 'bold')
c = colorbar; 
ylabel(c, {'Change in Activity', '(drugged - non drugged)'}, 'Rotation', 270, 'FontWeight', 'bold', 'FontSize', 14);
pos = c.Label.Position; 
c.Label.Position = [pos(1)+1, pos(2), pos(3)];
xlabel('Days Post-Myocardial Infarction', 'FontSize', 11, 'FontWeight', 'bold')
xtickangle(45)
set(gca, 'fontname', 'FixedWidth');
set(gca,'YTick',1:height(drugsToSimulate));
set(gca,'YTickLabel',rowLabels_drugs_clustered(:,1),'fontsize',11)
set(gca, 'fontname', 'FixedWidth');
title('Effect of Drug on MMP14', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Drugs', 'FontSize', 11, 'FontWeight', 'bold')

figure;
aSMAData = reshape(sens(:,87,:), [length(timePointsOfInterest) height(drugsToSimulate)])'; 
aSMAData_clustered = [];
for i = 1:height(drugsToSimulate)
    aSMAData_clustered = vertcat(aSMAData_clustered, aSMAData(str2double(clusteredRowLabels{i,:}),:));
end
cmaprange = [0.5:0.005:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
caxis([-0.25 0.25]); imagesc(aSMAData_clustered, [-0.25, 0.25]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(timePointsOfInterest));
set(gca,'XTickLabel',time,'fontsize',11)
title('Effect of Drug on aSMA', 'FontSize', 14, 'FontWeight', 'bold')
c = colorbar; 
ylabel(c, {'Change in Activity', '(drugged - non drugged)'}, 'Rotation', 270, 'FontWeight', 'bold', 'FontSize', 14);
pos = c.Label.Position; 
c.Label.Position = [pos(1)+1, pos(2), pos(3)];
xlabel('Days Post-Myocardial Infarction', 'FontSize', 11, 'FontWeight', 'bold')
xtickangle(45)
set(gca, 'fontname', 'FixedWidth');
set(gca,'YTick',1:height(drugsToSimulate));
set(gca,'YTickLabel',rowLabels_drugs_clustered(:,1),'fontsize',11)
set(gca, 'fontname', 'FixedWidth');
ylabel('Drugs', 'FontSize', 11, 'FontWeight', 'bold')

figure;
periostinData = reshape(sens(:,100,:), [length(timePointsOfInterest) height(drugsToSimulate)])'; 
periostinData_clustered = [];
for i = 1:height(drugsToSimulate)
    periostinData_clustered = vertcat(periostinData_clustered, periostinData(str2double(clusteredRowLabels{i,:}),:));
end
cmaprange = [0.5:0.005:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
caxis([-0.25 0.25]); imagesc(periostinData_clustered, [-0.25, 0.25]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(timePointsOfInterest));
set(gca,'XTickLabel',time,'fontsize',11)
title('Effect of Drug on Periostin', 'FontSize', 14, 'FontWeight', 'bold')
c = colorbar; 
ylabel(c, {'Change in Activity', '(drugged - non drugged)'}, 'Rotation', 270, 'FontWeight', 'bold', 'FontSize', 14);
pos = c.Label.Position; 
c.Label.Position = [pos(1)+1, pos(2), pos(3)];
xlabel('Days Post-Myocardial Infarction', 'FontSize', 11, 'FontWeight', 'bold')
xtickangle(45)
set(gca, 'fontname', 'FixedWidth');
set(gca,'YTick',1:height(drugsToSimulate));
set(gca,'YTickLabel',rowLabels_drugs_clustered(:,1),'fontsize',11)
set(gca, 'fontname', 'FixedWidth');
ylabel('Drugs', 'FontSize', 11, 'FontWeight', 'bold')


drugResponse = colImRNAData;