function [padding, list_numDrugs, list_drugAction, list_drugType, list_drugTarget, rowLabels] = createRowLabels(drugsToSimulate)
% Formats row labels for drug simulation heatmaps. 

rowLabels = {};
list_numDrugs = {};
list_drugAction = {};
list_drugType = {};
list_drugTarget = {};

for i = 1:height(drugsToSimulate)
    rowOfInterest = drugsToSimulate(i,:);
    if isempty(rowOfInterest.SimilarDrugs{1})
        numOfDrugs = 1;
    else
        numOfDrugs = length(strsplit(rowOfInterest.SimilarDrugs{1}, ';')) + 1; % Adds 1 to account for the drug that is not listed in "SimilarDrugs"
    end
    
    if numOfDrugs == 1
        list_numDrugs = [list_numDrugs; rowOfInterest.Drug];
    else
        list_numDrugs = [list_numDrugs; strcat(num2str(numOfDrugs), {' '}, 'drugs')];
    end
   
    drugAction = strcat(rowOfInterest.DrugAction{1});
    list_drugAction = [list_drugAction; drugAction];
    
    if strcmpi(rowOfInterest.IsAgonist{1}, 'Yes')
        drugType = 'Agonist';
        drugTarget = rowOfInterest.AgonistTargetGeneID{1};
    else
        drugType = 'Antagonist';
        drugTarget = rowOfInterest.AntagonistTargetGeneID{1};
    end
    
    if strcmpi(rowOfInterest.Drug{1}, 'Entresto') && strcmpi(rowOfInterest.IsAgonist{1}, 'Yes') && strcmpi(rowOfInterest.IsAntagonist{1}, 'Yes')
        drugType = 'Both';
        drugTarget = 'AT1R;NPRA';
    end
    
    list_drugType = [list_drugType; drugType];    
    list_drugTarget = [list_drugTarget; drugTarget];

end

max_numDrugs = max(cellfun('length', list_numDrugs));
max_drugAction = max(cellfun('length', list_drugAction));
max_drugType = max(cellfun('length', list_drugType));
max_drugTarget = max(cellfun('length', list_drugTarget));

for j = 1:length(list_numDrugs)
    % Number of drugs
    amount = strcat('%', num2str(max_numDrugs), 's');
    list_numDrugs{j} = sprintf(amount, list_numDrugs{j});
    
    % Drug action
    amount = strcat('%', num2str(max_drugAction), 's');
    list_drugAction{j} = sprintf(amount, list_drugAction{j});
    
    % Drug type
    amount = strcat('%', num2str(max_drugType), 's');
    list_drugType{j} = sprintf(amount, list_drugType{j});
    
    % Drug target
    if strcmp(list_numDrugs{j}, 'Marimastat')
        list_drugTarget{j} = 'MMP1,2,9,14';
    end
    padding = strcat('%', num2str(max_drugTarget), 's');
    list_drugTarget{j} = sprintf(padding, list_drugTarget{j});
end

for row = 1:length(list_numDrugs)
   rowLabels = [rowLabels; strcat(strcat(list_numDrugs{row}, {'|'}), strcat(list_drugAction{row}, {'|'}), strcat(list_drugType{row}, {'|'}), list_drugTarget{row})];
end

