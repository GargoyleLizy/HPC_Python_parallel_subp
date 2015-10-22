function [BHFBASol] = GetBHFBASol(model, list, targetRxn,substrateRxn, population, x,scores,BPCY,isGeneList,output_dir)
%GetOptGeneSol save the solution from optGene and optGeneR in same format as OptKnock
%
% [optGeneSol] = GetOptGeneSol(model, targetRxn, substrateRxn, generxnList, population, x, scores, isGeneList)
%
%INPUTS
% model
% targetRxn
% substrateRxn
% generxnList
% population
% x                 the best solution
% scores
% isGeneList
%
%OUTPUT
% optGeneSol

% Modification !!!
% writeDirect where the files should be saved
% writeDirect = 'C:\Users\Ewen\Documents\PhD\BAFBA2\Results\';
writeDirect = output_dir;

BHFBASol = struct();
BHFBASol.substrateRxn = substrateRxn;
BHFBASol.targetRxn = targetRxn;

if sum(x) == 0
    save (strcat(writeDirect, 'BHFBASol--target-', char(targetRxn),...
        '--sub-',char(substrateRxn),'--KOs-0-no_solution_better_than_WT'...
        ), 'BHFBASol')
    return;
end

%if isGeneList
%    BHFBASol.rxnList = list(logical(x));
%    BHFBASol.numDel = length(BHFBASol.rxnList);
%   [tmp,tmp2,BHFBASol.rxnList] = deleteModelGenes(model,BHFBASol.rxnList); %finds just the reactions that are KOed b/c of gene removal
%else
    BHFBASol.rxnList = list(logical(x));
    BHFBASol.numDel = length(BHFBASol.rxnList);

%end


BHFBASol.BPCY = -min(scores);
%ba_obj = min(scores);
BHFBASol.population = population;
BHFBASol.scores = scores;


[growthRate,minProd,maxProd] = testOptKnockSol(model,BHFBASol.targetRxn,BHFBASol.rxnList);
if (-BHFBASol.BPCY - maxProd) / maxProd < .001 % acculacy must be within .1%
    slnCheck = 'valid_sln';
else slnCheck = 'unsound_sln';
end
if (maxProd - minProd) / maxProd < .001 % acculacy must be within .1%
    slnType = 'unique_point';
else slnType = 'non_unique';
end
BPCYOK = maxProd * (.98^BHFBASol.numDel)* growthRate;
BHFBASol.sln=slnCheck;
BHFBASol.Genes = findGenesFromRxns(model,BHFBASol.rxnList);
BHFBASol.GR = growthRate;
% storage

    save (strcat(writeDirect, 'BHFBASol--rxns--target-', char(BHFBASol.targetRxn),...
        '--sub-',char(BHFBASol.substrateRxn),'--KOs-',num2str(BHFBASol.numDel),...
        '--BPCY-',num2str(BHFBASol.BPCY),...
        '--',slnCheck,'--',slnType,'--GR-',num2str(growthRate),...
        '--10CC.mat'...
        ), 'BHFBASol')