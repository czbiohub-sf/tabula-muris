%% cell_ontology_counts

clear all
close all
clc

%% load files

d = dir('*_annotation.csv');

for i = 1:length(d)
    
    workingFile = d(i).name;
    workingTable = readtable(workingFile, 'delimiter', ',');
    
    workingFile = strsplit(workingFile,'_');
    if length(workingFile) == 3
        tissueName = workingFile(1);
    elseif length(workingFile) == 4
        tissueName = strcat(workingFile(1:2));
    elseif length(workingFile) == 5
        tissueName = strcat(workingFile(1:3));
    end
    
    methodAux = workingFile(2:end);
    methodName = {};
    [a,b] = ismember('facs',workingFile);
    if a
        methodName = workingFile(b);
    else
        [a,b] = ismember('droplet',workingFile);
        if a
            methodName = workingFile(b);
        else methodName = {};
        end
    end
    
    if ~isempty(methodName)
        %if ~(strcmp(methodName,'droplet') & (strcmp(tissueName,'Heart')))
        [cellOntologies,ia,ic] = unique(workingTable.cell_ontology_class);
        [a,b]=hist(ic,unique(ic));
        
        if length(tissueName)>1 & strcmp(tissueName{2},'Non-Myeloid')
            tissueName{2} = 'NonMyeloid';
        end
        assignin('base',strjoin([tissueName, methodName],'_'), [cellOntologies num2cell(a')]);
%         end
    end
end
