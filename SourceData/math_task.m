function [] = math_task(proj_dir,output_dir,model,task_fun,rxn_list,targetRxn,substrateRxn,MaxKOs,imax)
% initialize the essential directory that needed for alg
% For future uses, users should rewrite these directories to adjust their 
% developing environment
projDir=proj_dir
cobraDir=strcat(projDir,'toolboxes/cobratoolbox/cobratoolbox-master/');
datasetDir=strcat(projDir,'SourceData/Dataset');
BAFBADir=strcat(projDir,'SourceData/BAFBA');
BHFBADir=strcat(projDir,'SourceData/BHFBA');
DBFBADir=strcat(projDir,'SourceData/DBFBA');

% ----- Initialize the essential Toolboxes ------
disp(pwd);
% cd into CobraToolbox directory and initialize the CobraToolbox
cd(cobraDir);
disp(pwd);
initCobraToolbox


% cd to the dataset directory
cd(datasetDir)
disp(pwd)

% ----Reading SBML model -------
% Read the SBML models,just copy from the CobraToolbox source code
% I donot know what is the meaning of 1000
% Just know that 1000 is the maximum flux. But is 1000 right for this task?
sbml_model = readCbModel(model,1000,'SBML');
% Ecoli=readCbModel('E coli iJR904',1000,'SBML');

% TODO --Need to find a way to read list from xls .
if(strcmp(rxn_list,'all'))
	rxnList=sbml_model.rxns;
else
	% rxnList should be fectched from a xls file
	% !!Warning, This xls read code has not tested yet. Will test it once I get the right xls file.
	[num,rxnList,raw]=xlsread(rxn_list);	
end 

switch task_fun
    case 'BAFBA'
        % Finally, start applying the algorithms
        cd(BAFBADir)
        disp(pwd)       
        BAFBA(sbml_model,rxnList,targetRxn,substrateRxn,MaxKOs,output_dir)
    case 'BHFBA'
        cd(BHFBADir)
	disp(pwd)
	BHFBA(sbml_model,rxnList,targetRxn,substrateRxn,MaxKOs,imax,output_dir)
    case 'DBFBA'
       	cd(DBFBADir)
	disp(pwd)
	DBFBA(sbml_model,rxnList,targetRxn,substrateRxn,MaxKOs,imax,output_dir)
    otherwise 
	disp('Task function name not recognized. Only BAFBA/BHFBA/DBFBA are allowed.') 
end
disp('Task finished? Great')
