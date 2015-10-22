function [DBFBASol] = GetDBFBASol(model, list, targetRxn,substrateRxn, population, x,scores,BPCY,isGeneList,output_dir)
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

% Modification
% writeDirect where the files should be saved
% riteDirect = 'C:\Users\Ewen\Documents\PhD\BAFBA2\Results\';
writeDirect = output_dir;

DBFBASol = struct();
DBFBASol.substrateRxn = substrateRxn;
DBFBASol.targetRxn = targetRxn;

if sum(x) == 0
    save (strcat(writeDirect, 'DBFBASol--target-', char(targetRxn),...
        '--sub-',char(substrateRxn),'--KOs-0-no_solution_better_than_WT'...
        ), 'DBFBASol')
    return;
end


    DBFBASol.rxnList = list(logical(x));
    DBFBASol.numDel = length(DBFBASol.rxnList);


DBFBASol.BPCY = -min(scores);
DBFBASol.population = population;
DBFBASol.scores = scores;


[growthRate,minProd,maxProd] = testOptKnockSol(model,DBFBASol.targetRxn,DBFBASol.rxnList);
if (-DBFBASol.BPCY - maxProd) / maxProd < .001 % acculacy must be within .1%
    slnCheck = 'valid_sln';
else slnCheck = 'unsound_sln';
end
if (maxProd - minProd) / maxProd < .001 % acculacy must be within .1%
    slnType = 'unique_point';
else slnType = 'non_unique';
end
BPCYOK = maxProd * (.98^DBFBASol.numDel)* growthRate;
DBFBASol.sln=slnCheck;
DBFBASol.Genes = findGenesFromRxns(model,DBFBASol.rxnList);
DBFBASol.GR = growthRate;
% storage

    save (strcat(writeDirect, 'DBFBASol--rxns--target-', char(DBFBASol.targetRxn),...
        '--sub-',char(DBFBASol.substrateRxn),'--KOs-',num2str(DBFBASol.numDel),...
        '--BPCY-',num2str(DBFBASol.BPCY),...
        '--',slnCheck,'--',slnType,'--GR-',num2str(growthRate),...
        '--10CC.mat'...
        ), 'DBFBASol')
