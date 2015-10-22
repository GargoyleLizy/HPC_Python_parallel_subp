function [BAFBASol] = GetBAFBASol(model, list, targetRxn,substrateRxn, population, x,scores,BPCY,isGeneList,output_dir)
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

% Modification : the output_dir would be used to write over writeDirect
% to make the algorithm platform independent.
writeDirect = output_dir;
% writeDirect where the files should be saved
%writeDirect = 'C:\Users\Ewen\Documents\PhD\BAFBA2\Results\';

BAFBASol = struct();
BAFBASol.substrateRxn = substrateRxn;
BAFBASol.targetRxn = targetRxn;

if sum(x) == 0
    save (strcat(writeDirect, 'BAFBASol--target-', char(targetRxn),...
        '--sub-',char(substrateRxn),'--KOs-0-no_solution_better_than_WT'...
        ), 'BAFBASol')
    return;
end

%if isGeneList
%    BAFBASol.rxnList = list(logical(x));
%    BAFBASol.numDel = length(BAFBASol.rxnList);
%   [tmp,tmp2,BAFBASol.rxnList] = deleteModelGenes(model,BAFBASol.rxnList); %finds just the reactions that are KOed b/c of gene removal
%else
    BAFBASol.rxnList = list(logical(x));
    BAFBASol.numDel = length(BAFBASol.rxnList);

%end


BAFBASol.BPCY = -min(scores);
%ba_obj = min(scores);
BAFBASol.population = population;
BAFBASol.scores = scores;


[growthRate,minProd,maxProd] = testOptKnockSol(model,BAFBASol.targetRxn,BAFBASol.rxnList);
if (-BAFBASol.BPCY - maxProd) / maxProd < .001 % acculacy must be within .1%
    slnCheck = 'valid_sln';
else slnCheck = 'unsound_sln';
end
if (maxProd - minProd) / maxProd < .001 % acculacy must be within .1%
    slnType = 'unique_point';
else slnType = 'non_unique';
end
BPCYOK = maxProd * (.98^BAFBASol.numDel)* growthRate;
BAFBASol.sln=slnCheck;
BAFBASol.Genes = findGenesFromRxns(model,BAFBASol.rxnList);
BAFBASol.GR = growthRate;
[isValidRxn,removeInd] = ismember(BAFBASol.rxnList,model.rxns);
removeInd = removeInd(isValidRxn);
BAFBASol.rxnNames =model.rxnNames(removeInd);

%[ExRxns,MaxTheoOut]= theoretMaxProd(model, 'pr_mol',targetRxn);
%[isValidRxn,removeInd] = ismember(targetRxn,ExRxns);
%maxtheo = MaxTheoOut(removeInd);
%diff_max = maxtheo - BAFBASol.BPCY;
%disp(sprintf('Difference between Theoretical Max Production & BPCY %02.0f',diff_max));    

% storage

    save (strcat(writeDirect, 'BAFBASol--rxns--target-', char(BAFBASol.targetRxn),...
        '--sub-',char(BAFBASol.substrateRxn),'--KOs-',num2str(BAFBASol.numDel),...
        '--BPCY-',num2str(BAFBASol.BPCY),...
        '--',slnCheck,'--',slnType,'--GR-',num2str(growthRate),...
        '--10CC.mat'...
        ), 'BAFBASol')
