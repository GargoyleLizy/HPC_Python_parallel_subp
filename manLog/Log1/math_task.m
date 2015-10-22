function [] = math_task(output_dir,model,task_fun)
% initialize the essential directory that needed for alg
% For future uses, users should rewrite these directories to adjust their 
% developing environment
cobraDir='/home/zhenyal/rProject/toolboxes/cobratoolbox/cobratoolbox-master/';
datasetDir='/home/zhenyal/rProject/SourceData/Dataset';
BAFBADir='/home/zhenyal/rProject/SourceData/BAFBA';

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
EcoliList=Ecoli.rxns;


switch task_fun
    case 'BAFBA'
        % Finally, start applying the algorithms
        cd(BAFBADir)
        disp(pwd)
        % This is just a test case to see if this really can works.
        BAFBA(sbml_model,EcoliList,'UMPK','UNK3',5,output_dir)
    case 'BHFBA'
        % place holder
    case 'DBFBA'
        % place holder
end
disp('Task finished? Great')
