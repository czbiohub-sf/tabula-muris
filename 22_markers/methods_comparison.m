clear all
close all
clc

%% goal

% the idea is to compare 10x with FACS with microwell-seq
% for this to automate start the files names with the tissue


%% load files

tissuesAvailable = {'Bladder', 'Kidney', 'Liver', 'Lung','Marrow', 'Limb_Muscle', 'Spleen', 'Thymus'}';
[s,v] = listdlg('PromptString','Select a tissue:',...
    'SelectionMode','single',...
    'ListString',tissuesAvailable);

tissuename = tissuesAvailable{s};filename1 = sprintf('facs_%s_cell_ontology_class_markers.csv',tissuename);
tissueFACS = readtable(filename1);
filename11 = sprintf('facs_%s_cell_ontology_class_classes.csv',tissuename);
tissueFACSnames = readtable(filename11, 'delimiter', ',');
filename2 = sprintf('droplet_%s_cell_ontology_class_markers.csv',tissuename);
tissue10x = readtable(filename2);
filename22 = sprintf('droplet_%s_cell_ontology_class_classes.csv',tissuename);
tissue10xnames = readtable(filename22, 'delimiter', ',');

if strcmp(tissuename,'Limb_Muscle')
    tissuenameHan = 'Muscle';
else tissuenameHan = tissuename;
end
filename3 = sprintf('%s_Han.csv',tissuenameHan);
tissueMicrowellSeq = readtable(filename3);

%% get tissue cell signatures for Tabula muris FACS
cellNamesTMfacs = unique(table2cell(tissueFACSnames(:,2)));
cellNamesTMfacs = sortrows(cellNamesTMfacs);
for i = 1:length(cellNamesTMfacs)
    [a,b] = ismember(table2cell(tissueFACSnames(:,2)),cellNamesTMfacs(i));
    cellClustersTMfacs(i,1:length(find(a))) = str2double(table2cell(tissueFACSnames(a,1)))';
end
cellClustersTMfacs = cellClustersTMfacs - ones;

tissueTMfacs = table2cell(tissueFACS);
for i = 1:size(cellClustersTMfacs,1)
    tissueTMfacsCellSigsAUX = [];
    tissueTMfacsCellSigsFCaux = [];
    for j = 1:size(cellClustersTMfacs,2)
        [m,n] = find(str2double(tissueTMfacs(:,7))==cellClustersTMfacs(i,j));
        
        tissueTMfacsCellSigsAUX = [tissueTMfacsCellSigsAUX;tissueTMfacs(m,8)];
        tissueTMfacsCellSigsFCaux = [tissueTMfacsCellSigsFCaux;tissueTMfacs(m,3)];
    end
    [tissueTMfacsCellSigsAUX,ia,ic] = unique(tissueTMfacsCellSigsAUX);
    tissueTMfacsCellSigs(1:length(tissueTMfacsCellSigsAUX)+1,i) = [cellNamesTMfacs(i);tissueTMfacsCellSigsAUX];
    tissueTMfacsCellSigsFC(1:length(tissueTMfacsCellSigsAUX)+1,i) = [-cellClustersTMfacs(i,j);tissueTMfacsCellSigsFCaux(ia)];
end


%% get tissue cell signatures for Tabula muris 10x
cellNamesTM10x = unique(table2cell(tissue10xnames(:,2)));
cellNamesTM10x = sortrows(cellNamesTM10x);
for i = 1:length(cellNamesTM10x)
    [a,b] = ismember(table2cell(tissue10xnames(:,2)),cellNamesTM10x(i));
        cellClustersTM10x(i,1:length(find(a))) = str2double(table2cell(tissue10xnames(a,1)))';
end
cellClustersTM10x = cellClustersTM10x - ones;

tissueTM10x = table2cell(tissue10x);
for i = 1:size(cellClustersTM10x,1)
    tissueTM10xCellSigsAUX = [];
    tissueTM10xCellSigsFCaux = [];
    for j = 1:size(cellClustersTM10x,2)
        [m,n] = find(str2double(tissueTM10x(:,7))==cellClustersTM10x(i,j));
        
        tissueTM10xCellSigsAUX = [tissueTM10xCellSigsAUX;tissueTM10x(m,8)];
        tissueTM10xCellSigsFCaux = [tissueTM10xCellSigsFCaux;tissueTM10x(m,3)];
    end
    [tissueTM10xCellSigsAUX,ia,ic] = unique(tissueTM10xCellSigsAUX);
    tissueTM10xCellSigs(1:length(tissueTM10xCellSigsAUX)+1,i) = [cellNamesTM10x(i);tissueTM10xCellSigsAUX];
    tissueTM10xCellSigsFC(1:length(tissueTM10xCellSigsAUX)+1,i) = [-cellClustersTM10x(i,j);tissueTM10xCellSigsFCaux(ia)];
end


%% get tissue cells signatures for microwell-seq

microwellSeq = table2cell(tissueMicrowellSeq);
microwellSeqcelltypes = {};
aux = 1;
microwellSeqCellSigs = cell(length(microwellSeq),length(unique(microwellSeq(:,8))));
for i = 1:length(microwellSeq)
    if isequal(microwellSeq(i,8),{[]})
        microwellSeq(i,8) = cellaux;
        microwellSeq{i,9} = aux;
    else
        cellaux = microwellSeq(i,8);
        aux = aux+1;
        microwellSeq{i,9} = aux;
        microwellSeqcelltypes = [microwellSeqcelltypes;cellaux];
        if i > 1
            microwellSeqCellSigs(1:(i-baux+1),length(microwellSeqcelltypes)-1) = [microwellSeq(baux,8);microwellSeq(baux:i-1,6)];
            microwellSeqcellSigsFC(1:(i-baux+1),length(microwellSeqcelltypes)-1) = [-cell2mat(microwellSeq(baux,5));cell2mat(microwellSeq(baux:i-1,2))];
        end
        baux = i;
    end
end
microwellSeqCellSigs(1:(i-baux+1),end) = [microwellSeq(baux,8);microwellSeq(baux:i-1,6)];
microwellSeqcellSigsFC(1:(i-baux+1),end+1) = [-cell2mat(microwellSeq(baux,5));cell2mat(microwellSeq(baux:i-1,2))];

[microwellSeqcelltypesSorted,index] = sortrows(microwellSeqcelltypes);
microwellSeqCellSigs = microwellSeqCellSigs(:,index);
microwellSeqcellSigsFC = microwellSeqcellSigsFC(:,index);


%% overlap tissue cell types signature

% pick a tissue to compare: bladder, (brain), kidney, liver, lung, marrow, muscle, pancreas, spleen, thymus

if strcmp(tissuename,'Bladder')
    tissueOverlaps = cell(1,6);
    tissueOverlapsSizes = [];
%     overlapMatrix = [[1;2;3],[1;2;3],[[1 0 0];[8 11 12];[5 6 14]]];% when using clusters
    overlapMatrix = [[1;2] [1;2] [[8 11 12];[14 0 0]]];% when using cell ontology
end

if strcmp(tissuename , 'Kidney')
    tissueOverlaps = cell(1,6);
    tissueOverlapsSizes = [];
    overlapMatrix = [[1;2;3;4],[1;5;3;6],[[8 0 0 0 0 0 0];[12 13 17 18 19 20 21];[5 10 11 0 0 0 0];[14 15 0 0 0 0 0]]];
end

if strcmp(tissuename,'Liver')
    tissueOverlaps = cell(1,6);
    tissueOverlapsSizes = [];
    overlapMatrix = [[3;4;1],[2;3;4],...
        [[5 0 0 0];[11 12 16 17];[1 2 19 20]]];
%     overlapMatrix = [[3;4;1],[2;3;4],...
%         [[5 0 0 0 0 0 0 0 0 0 0];[11 12 16 17 0 0 0 0 0 0 0];[1 2 3 4 8 9 10 14 15 19 20]]];
end

if strcmp(tissuename, 'Lung')
    tissueOverlaps = cell(1,6);
    tissueOverlapsSizes = [];
    overlapMatrix = [[1;3;4;8;10;11;12],[1;3;5;8;10;11;13],[[6 22 0 0 0 0 0 0 0 0 0];[15 32 0 0 0 0 0 0 0 0 0];[8 0 0 0 0 0 0 0 0 0 0];...
        [18 19 20 0 0 0 0 0 0 0 0];[10 11 12 13 14 17 21 23 24 26 28];[25 0 0 0 0 0 0 0 0 0 0];[29 30 31 0 0 0 0 0 0 0 0]]];
end

if strcmp(tissuename,'Marrow')
    tissueOverlaps = cell(1,6);
    tissueOverlapsSizes = [];
    overlapMatrix = [[8;14;17],[8;11;12],[[5 0];[6 7];[9 10]]];
end

if strcmp(tissuename,'Limb_Muscle')
    tissueOverlaps = cell(1,6);
    tissueOverlapsSizes = [];
    overlapMatrix = [[1;2;3;4;5;6],[1;3;4;5;6;7],[[1 2 0];[17 0 0];[4 0 0];[8 9 0];[16 0 0];[10 11 12]]];
end

if strcmp(tissuename , 'Spleen')
    tissueOverlaps = cell(1,6);
    tissueOverlapsSizes = [];
    overlapMatrix = [[1;2;3],[1;2;4],[[6];[11];[5]]];
end

if strcmp(tissuename , 'Thymus')
    tissueOverlaps = cell(1,6);
    tissueOverlapsSizes = [];
    overlapMatrix = [[1;2],[1;2],[[3 0];[4 2]]];
%     overlapMatrix = [[1;2;3],[1;2;3],[[3 0 0 0];[4 2 0 0];[5 6 7 8]]];
end

% find the overlaps
overlapGenesTissue = {};
for i = 1:size(overlapMatrix,1)
        
        % FACS
        aux1 = tissueTMfacsCellSigs(~cellfun('isempty',tissueTMfacsCellSigs(:,overlapMatrix(i,1))),overlapMatrix(i,1));
        aux = aux1(2:end);
        % 10x
        baux1 = tissueTM10xCellSigs(~cellfun('isempty',tissueTM10xCellSigs(:,overlapMatrix(i,2))),overlapMatrix(i,2));
        baux = baux1(2:end);
        % microwellSeq
        caux = [];
        for j = 3:size(overlapMatrix,2)
            if overlapMatrix(i,j)
                cauxaux = microwellSeqCellSigs(~cellfun('isempty',microwellSeqCellSigs(:,overlapMatrix(i,j))),overlapMatrix(i,j));
                cauxaux = cauxaux(2:end);
                caux = [caux;cauxaux];
            end
        end
        caux = unique(caux);
        
        % check overlaps for FACS
        [a1,b1] = ismember(aux,baux);
        overlap_FACS_10x = aux(a1);
        [a2,b2] = ismember(overlap_FACS_10x,caux);
        overlap_FACS_10x_microwell = overlap_FACS_10x(a2);
        overlap_FACS_10x_nomicrowell = overlap_FACS_10x(~a2);
        
        overlap_FACSno10x = aux(~a1);
        [a3,b3] = ismember(overlap_FACSno10x,caux);
        overlap_FACSno10x_microwell = overlap_FACSno10x(a3);
        overlap_FACSno10xnomicrowell = overlap_FACSno10x(~a3);
        
        % check overlaps for 10x
        [c1,d1] = ismember(baux,aux);
        overlap_10x_FACS = baux(c1);
        [c2,d2] = ismember(overlap_10x_FACS,caux);
        overlap_10x_FACS_microwell = overlap_10x_FACS(c2); % should be the same as overlap_FACS_10x_microwell
        overlap_10x_FACS_nomicrowell = overlap_10x_FACS(~c2); % should be the same as overlap_FACS_10xnomicrowell
        
        overlap_10xnoFACS = baux(~c1);
        [c3,d3] = ismember(overlap_10xnoFACS,caux);
        overlap_10xnoFACS_microwell = overlap_10xnoFACS(c3);
        overlap_10xnoFACSnomicrowell = overlap_10xnoFACS(~c3);
        
        % check overlaps for microwellSeq
        [e1,f1] = ismember(caux,aux);
        overlap_microwell_FACS = caux(e1);
        [e2,f2] = ismember(overlap_microwell_FACS,baux);
        overlap_microwell_FACS_10x = overlap_microwell_FACS(e2); % should be the same as overlap_FACS_10x_microwell
        overlap_microwell_FACSno10x = overlap_microwell_FACS(~e2); % should be the same as overlap_FACSno10x_microwell
        
        overlap_microwellnoFACS = caux(~e1);
        [e3,f3] = ismember(overlap_microwellnoFACS,baux);
        overlap_microwellnoFACS_10x = overlap_microwellnoFACS(e3); % should be the same as overlap_10xnoFACS_microwell
        overlap_microwellnoFACSno10x = overlap_microwellnoFACS(~e3);
            
        tissueOverlaps{1,i} = {overlap_FACS_10x_microwell; overlap_FACS_10x_nomicrowell; overlap_FACSno10x_microwell; overlap_FACSno10xnomicrowell; ...
            overlap_10xnoFACS_microwell; overlap_10xnoFACSnomicrowell; overlap_microwellnoFACSno10x};
            
        tissueOverlapsSizes = [tissueOverlapsSizes [length(overlap_FACS_10x_microwell);length(overlap_FACS_10x_nomicrowell);length(overlap_FACSno10x_microwell);length(overlap_FACSno10xnomicrowell);...
            length(overlap_10xnoFACS_microwell);length(overlap_10xnoFACSnomicrowell);length(overlap_microwellnoFACSno10x)]];
        
        % plot Venn diagrams
        vennNumbers = [length(overlap_FACS_10x_microwell);length(overlap_FACS_10x_nomicrowell);length(overlap_FACSno10x_microwell);length(overlap_FACSno10xnomicrowell);...
            length(overlap_10xnoFACS_microwell);length(overlap_10xnoFACSnomicrowell);length(overlap_microwellnoFACSno10x)];
        vennNumbers = vennNumbers([4,2,6,5,7,3,1]);
        eRR = vennX(vennNumbers,.1);
        title([sprintf('Tissue: %s',tissuename) sprintf(' -- CellType: %s',cellNamesTMfacs{overlapMatrix(i,1)})])
        saveas(gcf,[sprintf('Tissue_%s',tissuename) sprintf('_CellType_%s',cellNamesTMfacs{overlapMatrix(i,1)})],'pdf')

%         % plot FC distribution
%         figure()
%         % FACS + 10x + microwellSeq: overlap_FACS_10x_microwell
%         subplot(3,3,3)
%         [a,x] = ismember(aux1,overlap_FACS_10x_microwell);
%         FACSfc = cell2mat(tissueTMfacsCellSigsFC(a,overlapMatrix(i,1)));
%         hist(FACSfc,25)
%         hold on
%         [b,x] = ismember(baux1,overlap_FACS_10x_microwell);
%         tenxfc = cell2mat(tissueTM10xCellSigsFC(b,overlapMatrix(i,1)));
%         histogram(tenxfc,25)
%         microwellSeqfc = [];
%         for j = 3:size(overlapMatrix,2)
%             if overlapMatrix(i,j)
%                 caux1 = microwellSeqCellSigs(~cellfun('isempty',microwellSeqCellSigs(:,overlapMatrix(i,j))),overlapMatrix(i,j));
%                 [c,x] = ismember(caux1,overlap_FACS_10x_microwell);
%                 microwellSeqfc = [microwellSeqfc;microwellSeqcellSigsFC(c,overlapMatrix(i,j))];
%             end
%         end
%         histogram(microwellSeqfc,25)
%         legend;
%         hold off
%         
%         % FACS + 10x: overlap_FACS_10x_nomicrowell
%         subplot(3,3,2)
%         [a,x] = ismember(aux1,overlap_FACS_10x_nomicrowell);
%         FACSfc = cell2mat(tissueTMfacsCellSigsFC(a,overlapMatrix(i,1)));
%         hist(FACSfc,25)
%         hold on
%         [b,x] = ismember(baux1,overlap_FACS_10x_nomicrowell);
%         tenxfc = cell2mat(tissueTM10xCellSigsFC(b,overlapMatrix(i,1)));
%         histogram(tenxfc,25)
%         legend;
%         hold off
%         
%         % FACS + microwellSeq: overlap_FACSno10x_microwell
%         subplot(3,3,5)
%         [a,x] = ismember(aux1,overlap_FACSno10x_microwell);
%         FACSfc = cell2mat(tissueTMfacsCellSigsFC(a,overlapMatrix(i,1)));
%         hist(FACSfc,25)
%         hold on
%         microwellSeqfc = [];
%         for j = 3:size(overlapMatrix,2)
%             if overlapMatrix(i,j)
%                 caux1 = microwellSeqCellSigs(~cellfun('isempty',microwellSeqCellSigs(:,overlapMatrix(i,j))),overlapMatrix(i,j));
%                 [c,x] = ismember(caux1,overlap_FACSno10x_microwell);
%                 microwellSeqfc = [microwellSeqfc;microwellSeqcellSigsFC(c,overlapMatrix(i,j))];
%             end
%         end
%         histogram(microwellSeqfc,25)
%         legend;
%         hold off
%         
%         % FACs only: overlap_FACSno10xnomicrowell
%         subplot(3,3,1)
%         [a,x] = ismember(aux1,overlap_FACSno10xnomicrowell);
%         FACSfc = cell2mat(tissueTMfacsCellSigsFC(a,overlapMatrix(i,1)));
%         hist(FACSfc,25)
%         legend;
%         hold off
%         
%         % 10x + microwellSeq: overlap_10xnoFACS_microwell
%         subplot(3,3,8)
%         [b,x] = ismember(baux1,overlap_10xnoFACS_microwell);
%         tenxfc = cell2mat(tissueTM10xCellSigsFC(b,overlapMatrix(i,1)));
%         histogram(tenxfc,25)
%         hold on
%         microwellSeqfc = [];
%         for j = 3:size(overlapMatrix,2)
%             if overlapMatrix(i,j)
%                 caux1 = microwellSeqCellSigs(~cellfun('isempty',microwellSeqCellSigs(:,overlapMatrix(i,j))),overlapMatrix(i,j));
%                 [c,x] = ismember(caux1,overlap_10xnoFACS_microwell);
%                 microwellSeqfc = [microwellSeqfc;microwellSeqcellSigsFC(c,overlapMatrix(i,j))];
%             end
%         end
%         histogram(microwellSeqfc,25)
%         legend;
%         hold off
%         
%         % 10x only: overlap_10xnoFACSnomicrowell
%         subplot(3,3,4)
%         [b,x] = ismember(baux1,overlap_10xnoFACSnomicrowell);
%         tenxfc = cell2mat(tissueTM10xCellSigsFC(b,overlapMatrix(i,1)));
%         histogram(tenxfc,25)
%         legend;
%         hold off
%         % microwellSeq only: overlap_microwellnoFACSno10x
%         subplot(3,3,7)
%         microwellSeqfc = [];
%         for j = 3:size(overlapMatrix,2)
%             if overlapMatrix(i,j)
%                 caux1 = microwellSeqCellSigs(~cellfun('isempty',microwellSeqCellSigs(:,overlapMatrix(i,j))),overlapMatrix(i,j));
%                 [c,x] = ismember(caux1,overlap_microwellnoFACSno10x);
%                 microwellSeqfc = [microwellSeqfc;microwellSeqcellSigsFC(c,overlapMatrix(i,j))];
%                 hold on
%             end
%         end
%         histogram(microwellSeqfc,25)
%         legend;
%         hold off
        

        % write gene overlaps
        overlapGenes = cell(max([length(overlap_FACS_10x_microwell);length(overlap_FACS_10x_nomicrowell);length(overlap_FACSno10x_microwell);length(overlap_FACSno10xnomicrowell);...
            length(overlap_10xnoFACS_microwell);length(overlap_10xnoFACSnomicrowell);length(overlap_microwellnoFACSno10x)]),7);
        
        overlapGenesAux = {overlap_FACSno10xnomicrowell overlap_10xnoFACSnomicrowell overlap_microwellnoFACSno10x...
            overlap_FACS_10x_nomicrowell overlap_FACSno10x_microwell overlap_10xnoFACS_microwell...
            overlap_FACS_10x_microwell};
        
        for k = 1:length(overlapGenesAux)
            overlapGenes(1:length(overlapGenesAux{k}),k) = overlapGenesAux{k};
        end
        
        overlapGenes = [{'FACS','Droplet10x','microwellSeq','FACS_10x','FACS_microwellSeq','microwellSeq_10x','FACS_10x_microwellSeq'};overlapGenes];
        [oga,ogb] = size(overlapGenes);
        aux = cell(oga,2);
        aux(1:2,2) = ['Cell Ontology';cellNamesTMfacs(overlapMatrix(i,1))];
        if ~isempty(overlapGenesTissue)
            [ogta,ogtb] = size(overlapGenesTissue);
            if oga > ogta
                overlapGenesTissue = [overlapGenesTissue;cell(oga-ogta,ogtb)];
                overlapGenesTissue = [overlapGenesTissue aux overlapGenes];
            else
%                 overlapGenesTissue = [overlapGenesTissue cell(ogta,ogb+2)];
                overlapGenesTissue(1:oga,end+1:end+2+ogb) = [aux overlapGenes];
            end
        else overlapGenesTissue = [aux overlapGenes];
        end
            
        

end


% write to table overlapping genes in the Venn diagrams
overlapGenesTissue = cell2table(overlapGenesTissue);%,'VariableNames',{'FACS','Droplet10x','microwellSeq','FACS_10x','FACS_microwellSeq','microwellSeq_10x','FACS_10x_microwellSeq'});
%         fileNameOverlapGenes = ;
writetable(overlapGenesTissue,sprintf('%s_genes.csv',tissuename))

    
% write to table cell ontology matches
overlapMatrixNames = cell(size(overlapMatrix));
overlapMatrixNames(:,1) = cellNamesTMfacs(overlapMatrix(:,1));
overlapMatrixNames(:,2) = cellNamesTM10x(overlapMatrix(:,2));
for i=1:size(overlapMatrix,1)
    for j = 3:size(overlapMatrix,2)
        if overlapMatrix(i,j)
            overlapMatrixNames(i,j) = microwellSeqcelltypesSorted(overlapMatrix(i,j));
        end
    end
end

microwellSeqVars = {'microwellSeq1','microwellSeq2','microwellSeq3','microwellSeq4',...
    'microwellSeq5','microwellSeq6','microwellSeq7','microwellSeq8','microwellSeq9',...
    'microwellSeq10','microwellSeq11','microwellSeq12','microwellSeq13','microwellSeq14'};
overlapMatrixNames = cell2table(overlapMatrixNames,'VariableNames',{'FACS','Droplet10x',microwellSeqVars{1:size(overlapMatrix,2)-2}});
writetable(overlapMatrixNames,sprintf('%s_cell_ontologies.csv',tissuename))





