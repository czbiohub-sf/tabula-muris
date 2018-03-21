%% averages_nreads_ngenes

clear all
close all
clc

%% load files

d = dir('*_nreads_ngenes.csv');
droplet_nreads_ngenes = cell(length(d),3);

acrossAllTissues = [];
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
    
    if length(tissueName)>1
        if strcmp(tissueName{2},'Non-Myeloid')
            tissueName{2} = 'NonMyeloid';
        end
        tissueName = {strjoin(tissueName,'_')};
    end
    
    droplet_nreads_ngenes(i,:) = [tissueName mean(workingTable.nGene)...
        mean(workingTable.nUMI)];
    acrossAllTissues = [acrossAllTissues;[workingTable.nGene workingTable.nUMI]];
    
end
droplet_nreads_ngenes = [{'acrossAllTissues' num2cell(mean(acrossAllTissues(:,1))) num2cell(mean(acrossAllTissues(:,2)))};droplet_nreads_ngenes];
droplet_nreads_ngenes = cell2table(droplet_nreads_ngenes,'VariableNames',{'TissueName' 'avg_nGene' 'avg_nUMI'});
writetable(droplet_nreads_ngenes,'dropletNreadsNgenes.csv')

