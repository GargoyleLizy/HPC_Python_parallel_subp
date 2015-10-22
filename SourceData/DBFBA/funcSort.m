function [inP1 oP4 oP5] = funcSort(fit, bPos,  dim, n)
           %%%%%%%%%%dim = dimension of function;
          
inP1      = fit;

inP2      = bPos;

[a1 a2]   = size(inP1);
[b1 b2]   = size(inP2);
         
           inP1 = fit;
           [sortP1 IN] = sort(fit,'ascend');
           
           oP1= zeros((a1+ 1),a2);
          for j=1:a2
              
          oP1(1,j)     =IN(1,j);
          oP1(2,j)     =sortP1(:,j);
          end
          
          oP1;
         
          oP2= zeros(b1+ 1 ,b2);
         
          k1=0;
                   
          for k=1:b2
             k1=k1+1; 
              
          oP2(1,k)                   = k1;
          
          oP2((2:(1 + b1)),k)     =  bPos(:,k);
          end
          oP2;
            oP3=zeros(b1+ 1 ,b2);    
              
              for j=1:a2
                                        
                 oP3(:,j)=oP2(:,oP1(1,j));
            
                
          end
          oP2=oP3;
          oP4=oP1((a1+1),:);
          oP5=oP2(2:(b1+1),:);
                  
            
          
          
          
         
          
           
