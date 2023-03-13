%
%   Graph of Shock  
%
%
%


case_data_rich = 3;  % Case A =1, Case B = 2,  Case C =3 

if case_data_rich == 1
    load  shock_sample_A.txt;
    
     AAA=shock_sample_A; 
     
     figure(1)
     
end     

if case_data_rich == 2
     load  shock_sample_B.txt;  
     
      AAA=shock_sample_B; 
     
      figure(2)
end     


if case_data_rich == 3
    load  shock_sample_C.txt;
    
     AAA=shock_sample_C;
     
      figure(3)
end     


n = 7;
AA=zeros(59,3*n);

titlename={'Preference Shock ' ;
           'Investment Shock ' ; 
           'Equity Premium Shock ' ;
            'Labor Supply Shock ' ;            
           'Productivity Shock ' ;                          
           'Government Spending Shock ' ;
           'Monetary Policy Shock ' };

for i=1:n
    AA(:,1+(i-1)*3:3+(i-1)*3)=AAA((i-1)*59+1:(i-1)*59+59,1:3);
end


    %figure(1)
    
for j=1:7  
   subplot(3,3,j)
   hold on
 
  for i=0:0.01:1
       plot(1981.25:0.25:1995.75,AA(:,3+3*(j-1))-i*(AA(:,3+3*(j-1))-AA(:,2+3*(j-1))),'Color',[0.753 0.753 0.753],'LineWidth',2)    
         %plot(1981.25:0.25:1995.75,AA(:,3),'r--')
        %plot(1981.25:0.25:1995.75,AA(:,4),'r--')
   end
         plot(1981.25:0.25:1995.75,AA(:,1+3*(j-1)),'b','LineWidth',2)
         

  title( titlename(j),'FontSize',14 )  
            %title('éYèo','FontSize',14)
  xlim([min(1981.25),max(1995.75)])
  set(gca,'fontsize',12);
  hold off
      
end
 

