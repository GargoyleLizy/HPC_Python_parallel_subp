function [BAFBASol] = BAFBA(model,list,targetRxn,substrateRxn,MaxKOs,output_dir)

% Modification for Parallel.
% Added extra variable "output_dir", which is the path to the results Dir.
% 

%model_temp = input('Please Enter Model Name = ','s');
%model = readCbModel(model_temp);
%reply  =  input('Do  you  want to Enter Possible Knockout List?  Y/N  :  ',  's');

%if (reply == 'Y' || reply == 'y')

%list  = input('Please Enter Possible Knockout List = ','s');
%[num,list,raw]= xlsread(list);
%else
%    list = model.rxns;
%end
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
e    =   3;       % Number of elite selected points (e.g. 10-50)

%nsp  =    8;   % Number of Bees around each selected locations ( except the ellit location )
nsp = round(m/2);
%nep  =    6;      % Number of Bees around each elite locations
nep = (e * 2);
ngh  =    1;      % Patch radius for neighbourhood search
sc   =    2;	  % Shrinking constant; defined as percentage (%) and range is between 0-100
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


e = MaxKnockOuts;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long


aa                 = double(zeros(1,pop));
val               = double(zeros(1,n));


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
       val(i) = -minProd * (.98^nummutations)* growthrate;
       
end
            %%%%%%%%%%%%End of initial random generation 

			
                
                %%%%%%%Run until maximum number of iteration met%%%%%%%%%

   
for iter = 1 : imax
disp(sprintf('Iteration Number: %02.0f',iter));    
				
				%%%% Sorting fitnesses & positions
                
				[fit sortedFit bPosSort]= funcSort(val, bPos, dim, n);

                          
            
			%%%%%%%%%% Choosing best m
            
				for(i= 1 : m)
					
                    for(d= 1 : dim)
						candidx(d,i) = double(bPosSort(d,i));
                        

                    end 
                end
              
		     %%%%%%%%%Recruitment stage  
				for(i= 1 : m)
				 
                    
					 if(i<=e)
						aa(i)=nep; %%% Number of bees around each elite sites
					 else
						aa(i)=nsp; %%% Number of bees around other selected sites
                     end
                end
                
                
               %%%%%% Search inrando the neighbourhood
                              
                for(k= 1 : m)           %%%%%k site
				
					for(j= 1: aa(k))     %%%%% j recruited bee
					
						
                        for(d= 1 : dim)  %%%%%d dimension
						 
                              
                                      bNghPos(d,j)=(randfunc(candidx(d,k)-nghx , candidx(d,k)+nghx));
                                     
                                      if  (bNghPos(d,j) > 0.5 )
                                        
                                          bNghPos(d,j)=1;
                                      else
                                          bNghPos(d,j) = 0;
                                                                              
                                      end
                                  
                                                            
                                 double(bNghPos(d,j));

                         end
                         
                             
                        

						 %bNghFit(1,j)=(func(bNghPos,j));
                         rxn_vectorNgh = bNghPos(:,j);
                         rxnListNgh = list(logical(rxn_vectorNgh));
    
                        nummutationsNgh = sum(rxn_vectorNgh);
                        if nummutationsNgh > MaxKnockOuts
                            continue;
                        end
    
                    [isValidRxnNgh,removeIndNgh] = ismember(rxnListNgh,model.rxns);
                    removeIndNgh = removeIndNgh(isValidRxnNgh);
                    modelKONgh = model;
                    modelKONgh.ub(removeIndNgh) = 0;
                    modelKONgh.lb(removeIndNgh) = 0;
                
                    FBAsolutionNgh = optimizeCbModel(modelKONgh);
                    growthrateNgh = FBAsolutionNgh.f;
                
                    if growthrateNgh < .010
                        continue;
                    end 
        
                modelFixGRNgh = modelKONgh;
                gamma = 0.9*growthrate;
                new_grNgh = FBAsolutionNgh.f - gamma;
                modelFixGRNgh = changeRxnBounds(modelKONgh,modelKONgh.rxns(modelKONgh.c==1),new_grNgh,'l');
                modelFixGRNgh = changeObjective(modelFixGRNgh,targetRxn);
                solMinNgh = optimizeCbModel(modelFixGRNgh,'max');
                minProdNgh = solMinNgh.f;
                bNghFit(j) = -minProdNgh * (.98^nummutations)* growthrate;
       
 					
						 if(bNghFit(j)< sortedFit(k))
							
                             bPosSort(:,k)=double(bNghPos(:,j));
                             
                             candidx(:,k)=double(bNghPos(:,j));
                             
                             sortedFit(1,k)=double(bNghFit(1,j));
                        
% 							                        
                         end
                    end
                    
                end
                       
              
					%%%end of  Neighbourhood Search

				%%%Shrink all the patches using the shrinking constant (sc) variable
				nghx=nghx*((100-sc)/100);

				%%% Send rest of the bees for random search...
                %[a1 a2] = size(sortedFit);
                %ranSearchBees= a2 -m;    %%%%% Number of bees for random search 
				ranSearchBees= n -m;

                
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

%template(iter,:)        = double(-sortedFit(1,:));                       
                                                             
v(iter,:)=-sortedFit(1,1);      

       
end
Final_Result = v(1:iter,1);
plot(Final_Result)
title('A Hybrid of Bees Algorithm and Flux Balance Analysis')
xlabel('Generation Number')
ylabel('Objective Function')

% Addition: Need a FileName to specify the result of this task
max_str=num2str(MaxKOs);
task_name=strcat('BAFBA-',targetRxn,'-',substrateRxn,'-',max_str);
% Could just save pictures to output dir instead of direct to it.
% But then we bound ourselves to specific platforms.
% !Debug: I forgot that direct the directory to output directory would cause the 
% 	Matlab unable to understand "GetBAFBASol" function later. Hence, just comment this part,
%	may be improve this part later.
% cd(output_dir);

% Modification : need to rewrite the absolute directory path to a more
% flexiable place.
%saveas(gcf, 'C:\Users\Ewen\Documents\PhD\BAFBA2\Results\Objective_Function', 'jpg')
objective_pic_name = strcat(task_name,'-Objective_Function');
saveas(gcf,objective_pic_name,'jpg');


minProdFixGR = sortedFit(1,1);
Scores = sortedFit;
BPCY = sortedFit(1,1);
t = bPosSort(:,1);
r=list(logical(t));
Best = t;
% Modification : Since the GetBAFBASol function contains the save part, 
% The output_dir 
BAFBASol = GetBAFBASol(model, list, targetRxn,substrateRxn, bPos,Best, Scores,BPCY,rxnok,output_dir);

[biomassValues,targetValues,lineHandle] = productionEnvelope(model,r,'k',targetRxn,model.rxns(model.c==1));

% Modification : Same as previous picture
%saveas(gcf, 'C:\Users\Ewen\Documents\PhD\BAFBA2\Results\Production_Envelope', 'jpg')
production_pic_name = strcat(task_name,'-Production_Envelope');
saveas(gcf,production_pic_name,'jpg');

End_Time = etime(clock,tstart)
end



