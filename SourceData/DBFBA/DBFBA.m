function [DBFBASol] = DBFBA(model,list,targetRxn,substrateRxn,MaxKOs,imax,output_dir)

% Modfication: The interactive part are rewritten to allow parallel processing

%model_temp = input('Please Enter Model Name = ','s');
%model = readCbModel(model_temp);
%reply  =  input('Do  you  want to Enter Possible Knockout List?  Y/N  :  ',  's');
%
%if (reply == 'Y' || reply == 'y')
%
%list  = input('Please Enter Possible Knockout List = ','s');
%[num,list,raw]= xlsread(list);
%else   
%    list = model.rxns;
%end
%
%targetRxn = input('Please Enter Target Reaction = ','s');
%substrateRxn = input('Please Enter Substrate Reaction = ','s');
%MaxKOs = input('Please Enter Maximum Gene to Knockout = ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imax =100;   % number of iterations (e.g. 1000-5000)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim               = length(list); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% The Bees Algorithm%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n                 = round(dim + (dim/4));      % Number of scout bees (e.g. 40-100)
pop = n;
m    =  15;       % Number of selected locations  (e.g. 10-50)
%e    =   3;       % Number of elite selected points (e.g. 10-50)
e = MaxKOs;
%nsp  =    8;   % Number of Bees around each selected locations ( except the ellit location )
nsp = round(m/2);
%nep  =    6;      % Number of Bees around each elite locations
nep = (e * 2);
ngh  =    30;      % Patch radius for neighbourhood search
sc   =    1;	  % Shrinking constant; defined as percentage (%) and range is between 0-100
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global HTABLE
HTABLE = java.util.Hashtable;
global MaxKnockOuts

%%%%%%%%%%%%%%%%%
%Reactions/Genes%
%%%%%%%%%%%%%%%%%
rxnok = 1;
geneok = 1;
for i = 1:dim
    if(~ ismember(list{i}, model.rxns)),  rxnok = 0; end
    if(~ ismember(list{i}, model.genes)),geneok = 0; end
end
if geneok
    display('assuming list is genes');
elseif rxnok
    display('assuming list is reactions');
else
    display('list appears to be neither genes nor reactions:  aborting');
    return;
end


    MaxKnockOuts = MaxKOs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long


aa                 = double(zeros(1,pop));

nghx=ngh;              %%define patch size
tstart = clock;
	%%%%%%%%Initial Random generation %%%%%%
for(i=1:n)
               
	perm= randperm(dim);
    for k =1:MaxKnockOuts
     bPos(perm(k),i)=1;
    end		
end
    for(i=1:n)
               rxn_vector = bPos(:,i);
               rxnList = list(logical(rxn_vector));
    
        nummutations = sum(rxn_vector);
        if nummutations > MaxKnockOuts
           continue;
        end
    
        [isValidRxn,removeInd] = ismember(rxnList,model.rxns);
        removeInd = removeInd(isValidRxn);
        modelKO = model;
        modelKO.ub(removeInd) = 0;
        modelKO.lb(removeInd) = 0;
                
        FBAsolution = optimizeCbModel(modelKO);
        growthrate = FBAsolution.f;       
        
        if growthrate < .010
            continue;
        end 
        
       modelFixGR = modelKO;
       gamma = 0.9*growthrate;
       new_gr = FBAsolution.f - gamma;
       modelFixGR = changeRxnBounds(modelKO,modelKO.rxns(modelKO.c==1),new_gr,'l');
       modelFixGR = changeObjective(modelFixGR,targetRxn);
       solMin = optimizeCbModel(modelFixGR,'max');
       minProd = solMin.f;

       fit(i) = -minProd * (.98^nummutations)* growthrate;
          
end
            %%%%%%%%%%%%End of initial random generation 

			
                
                %%%%%%%Run until maximum number of iteration met%%%%%%%%%

   
for iter = 1 : imax
disp(sprintf('Iteration Number: %02.0f',iter));    
				
				%%%% Sorting fitnesses & positions
                
				[fit sortedFit bPosSort]= funcSort(fit, bPos, dim, n);

                          
            
			%%%%%%%%%% Choosing best m
            
				for(j= 1 : m)
					
                    for(d= 1 : dim)
						candidx(d,j) = double(bPosSort(d,j));
                    end 
                end
                           
                
               %%%%%% Search inrando the neighbourhood
               
pop = candidx(:,1:m);
popold    = zeros(size(pop));     % toggle population
val       = zeros(1,m);          % create and reset the "cost array"
bestmem   = zeros(1,dim);           % best population member ever
bestmemit = zeros(1,dim);           % best population member in iteration
nfeval    = 0;                    % number of function evaluations
CR = 0.5;
F = 1;
refresh = 5;
refresh = floor(refresh);
%------Evaluate the best member after initialization----------------------

ibest   = 1;                      % start with first population member
for i=1:m 
               rxn_vector = pop(:,i);
               rxnList = list(logical(rxn_vector));
    
        nummutations = sum(rxn_vector);
        if nummutations > MaxKnockOuts
           continue;
        end
    
        [isValidRxn,removeInd] = ismember(rxnList,model.rxns);
        removeInd = removeInd(isValidRxn);
        modelKO = model;
        modelKO.ub(removeInd) = 0;
        modelKO.lb(removeInd) = 0;
                
        FBAsolution = optimizeCbModel(modelKO);
        growthrate2 = FBAsolution.f;  
        
        if growthrate2 < .010
            continue;
        end 
        
       modelFixGR = modelKO;
       gamma = 0.9*growthrate;
       new_gr = FBAsolution.f - gamma;
       modelFixGR = changeRxnBounds(modelKO,modelKO.rxns(modelKO.c==1),new_gr,'l');
       modelFixGR = changeObjective(modelFixGR,targetRxn);
       solMin = optimizeCbModel(modelFixGR,'max');
       minProd = solMin.f;
      val(i) = -minProd * (.98^nummutations)* growthrate;

bestval = val(1);                 % best objective function value so far
nfeval  = nfeval + 1;

  if (val(i) < bestval)           % if member is better
     ibest   = i;                 % save its location
     bestval = val(i);
  end   
end
bestmemit = pop(: , ibest);         % best member of current iteration
bestvalit = bestval;              % best value of current iteration

bestmem = bestmemit;              % best member ever

%------DE-Minimization---------------------------------------------
%------popold is the population which has to compete. It is--------
%------static through one iteration. pop is the newly--------------
%------emerging population.----------------------------------------

pm1 = zeros(dim, m);              % initialize population matrix 1
pm2 = zeros(dim, m);              % initialize population matrix 2
pm3 = zeros(dim, m);              % initialize population matrix 3
pm4 = zeros(dim, m);              % initialize population matrix 4
pm5 = zeros(dim, m);              % initialize population matrix 5
bm  = zeros(dim, m);              % initialize bestmember  matrix
ui  = zeros(dim, m);              % intermediate population of perturbed vectors
mui = zeros(dim, m);              % mask for intermediate population
mpo = zeros(dim, m);              % mask for old population
rot = (0:1:m-1);               % rotating index array (size NP)
rotd= (0:1:dim-1);                % rotating index array (size D)
rt  = zeros(m);                % another rotating index array
rtd = zeros(dim);                 % rotating index array for exponential crossover
a1  = zeros(m);                % index array
a2  = zeros(m);                % index array
a3  = zeros(m);                % index array
a4  = zeros(m);                % index array
a5  = zeros(m);                % index array
ind = zeros(4);



  popold = pop;                   % save the old population

  ind = randperm(4);              % index pointer array

  a1  = randperm(m);             % shuffle locations of vectors
  rt = rem(rot+ind(1),m);        % rotate indices by ind(1) positions
  a2  = a1(rt+1);                 % rotate vector locations
  rt = rem(rot+ind(2),m);
  a3  = a2(rt+1);                
  rt = rem(rot+ind(3),m);
  a4  = a3(rt+1);               
  rt = rem(rot+ind(4),m);
  a5  = a4(rt+1);                

  pm1 = popold(:,a1);             % shuffled population 1
  pm2 = popold(:, a2);             % shuffled population 2
  pm3 = popold(:, a3);             % shuffled population 3
  pm4 = popold(:, a4);             % shuffled population 4
  pm5 = popold(:, a5);             % shuffled population 5

  for i=1:m                      % population filled with the best member
    bm(:,i) = bestmemit;          % of the last iteration
  end

  mui = rand(dim, m) < CR;          % all random numbers < CR are 1, 0 otherwise

    ui = pm3 + F*(pm1 - pm2);       % differential variation
    ui = popold.*mpo + (ui.*mui).^2;     % crossover


%-----Select which vectors are allowed to enter the new population------------
  for i=1:m
      rxn_vector = ui(:,i);
      rxnList = list(logical(rxn_vector));
    
        nummutations = sum(rxn_vector);
        if nummutations > MaxKnockOuts
           continue;
        end
    
        [isValidRxn,removeInd] = ismember(rxnList,model.rxns);
        removeInd = removeInd(isValidRxn);
        modelKO = model;
        modelKO.ub(removeInd) = 0;
        modelKO.lb(removeInd) = 0;
                
        FBAsolution = optimizeCbModel(modelKO);
        growthrate3 = FBAsolution.f;
        
        if growthrate3 < .010
            continue;
        end 
        
       modelFixGR = modelKO;
       gamma = 0.9*growthrate;
       new_gr = FBAsolution.f - gamma;
       modelFixGR = changeRxnBounds(modelKO,modelKO.rxns(modelKO.c==1),new_gr,'l');
       modelFixGR = changeObjective(modelFixGR,targetRxn);
       solMin = optimizeCbModel(modelFixGR,'max');
       minProd = solMin.f;
       tempval(i) = -minProd * (.98^nummutations)* growthrate;

  end
  temp = tempval(1,:);
  for i=1:temp
    nfeval  = nfeval + 1;
    if (tempval(i) <= val(i))  % if competitor is better than value in "cost array"
       pop(:, i) = ui(:,i);  % replace old vector with new one (for new iteration)
       val(i)   = tempval(i);  % save value in "cost array"
    
       %----we update bestval only in case of success to save time-----------
       if (min(tempval) < bestval)     % if competitor better than the best one ever
          bestval = min(tempval);      % new best value
          bestmem = ui(:,i);      % new best parameter vector ever
       end
    end
  end %---end for imember=1:NP

  bestmemit = bestmem;       % freeze the best member of this iteration for the coming 
                             % iteration. This is needed for some of the strategies.

%----Output section----------------------------------------------------------

  if (refresh > 0)
    if (rem(iter,refresh) == 0)
  %     fprintf(1,'Iteration: %d,  Best: %+5.2d,  F: %f,  CR: %f,  NP: %d\n',iter,bestval,F,CR,m);
  %     for n=1:m
  %       fprintf(1,'best(%d) = %f\n',m,bestmem(m));
  %     end
    end
  end

  if (val(1,m) < sortedFit(1,m))
        bPosSort(:,m)   =  ui(:,m);
     sortedFit(1,m)  =  val(1,m);
  end

  
  %%%end of  Neighbourhood Search

				%%%Shrink all the patches using the shrinking constant (sc) variable
				nghx=nghx*((100-sc)/100);

				%%% Send rest of the bees for random search...
                [a1 a2] = size(sortedFit);
                ranSearchBees= a2 -m;    %%%%% Number of bees for random search 
				%ranSearchBees= n - m;

                
                for(k= 1 : ranSearchBees)	
				
                    perm= randperm(dim);
                    for e =1:MaxKnockOuts
                        bPosSortcand(perm(e),k)=1;
                    end
                end	
                
                for(k= 1 : ranSearchBees)
                           rxn_vectorSortcand = bPosSortcand(:,k);
                           rxnListSortcand = list(logical(rxn_vectorSortcand));
    
                            nummutationsSortcand = sum(rxn_vectorSortcand);
                            if nummutationsSortcand > MaxKnockOuts
                                continue;
                            end
    
                    [isValidRxnSortcand,removeIndSortcand] = ismember(rxnListSortcand,model.rxns);
                    removeIndSortcand = removeIndSortcand(isValidRxnSortcand);
                    modelKOSortcand = model;
                    modelKOSortcand.ub(removeIndSortcand) = 0;
                    modelKOSortcand.lb(removeIndSortcand) = 0;
                
        FBAsolutionSortcand = optimizeCbModel(modelKOSortcand);
        growthrateSortcand = FBAsolutionSortcand.f;
                
        if growthrateSortcand < .010
            continue;
        end 
        
       modelFixGRSortcand = modelKOSortcand;
       gamma = 0.9*growthrate;
       new_grSortcand = FBAsolutionSortcand.f - gamma;
       modelFixGRSortcand = changeRxnBounds(modelKOSortcand,modelKOSortcand.rxns(modelKOSortcand.c==1),new_grSortcand,'l');
       modelFixGRSortcand = changeObjective(modelFixGRSortcand,targetRxn);
       solMinSortcand = optimizeCbModel(modelFixGRSortcand,'max');
       minProdSortcand = solMinSortcand.f;
       sortedFitcand(k) = -minProdSortcand * (.98^nummutations)* growthrate;
       
   
                          
                           if (sortedFitcand(1,k) < sortedFit(1,k+m))
                              bPosSort(:,m+k)   =  bPosSortcand(:,k);
                              sortedFit(1,k+m)  =  sortedFitcand(1,k);
                              
                           end
                           
                end
				
			
                
[fit sortedFit bPosSort]= funcSort(sortedFit, bPosSort, dim, n);


template(iter,:)        = double(-sortedFit(1,:));                      
      
sortedFit(1,1);
   
       
end


            
minProdFixGR = sortedFit(1,1);
Scores = sortedFit;
BPCY = sortedFit(1,1);
t = bPosSort(:,1);
r=list(logical(t));
Best = t;
DBFBASol = GetDBFBASol(model, list, targetRxn,substrateRxn, bPos,Best, Scores,BPCY,rxnok,output_dir);

sort_temp = sort(template);
Final_Result = sort_temp(1:iter,1);
plot(Final_Result)
title('Differential Bees Flux Balance Analysis')
xlabel('Generation Number')
ylabel('Objective Function')

% Modification: should save pictures to output_dir
%saveas(gcf, 'C:\Users\Ewen\Documents\PhD\BAFBA2\Results\Objective_Function', 'jpg')
[biomassValues,targetValues,lineHandle] = productionEnvelope(model,r,'k',targetRxn,model.rxns(model.c==1));
%saveas(gcf, 'C:\Users\Ewen\Documents\PhD\BAFBA2\Results\Production_Envelope', 'jpg')
End_Time = etime(clock,tstart)
end



