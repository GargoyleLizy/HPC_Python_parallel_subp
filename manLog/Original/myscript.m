function [] = myscript()
% initialize the essential directory that needed for alg
% For future uses, users should rewrite these directories to adjust their 
% developing environment
cobraDir='/home/zhenyal/rProject/toolboxes/cobratoolbox/cobratoolbox-master/'
datasetDir='/home/zhenyal/rProject/SourceData/Dataset'
BAFBADir='/home/zhenyal/rProject/SourceData/BAFBA'

disp(pwd)

% cd into CobraToolbox directory and initialize the CobraToolbox
cd(cobraDir)
disp(pwd) 
initCobraToolbox
% could test the initialization of CobraToolbox with OptKnock
% But whats the point? maybe later

% cd to the dataset directory
cd(datasetDir)
disp(pwd)
% Read the SBML models,just copy from the CobraToolbox source code
% I donot know what is the meaning of 1000
% Just know that 1000 is the maximum flux. But is 1000 right for this task?
Ecoli=readCbModel('E coli iJR904',1000,'SBML')
EcoliList=Ecoli.rxns

% Finally, start applying the algorithms
cd(BAFBADir)
disp(pwd)
% This is just a test case to see if this really can works.
BAFBA(Ecoli,EcoliList,'UMPK','UNK3',5)
disp('BAFBA finished? Then you should rewrite it.')
