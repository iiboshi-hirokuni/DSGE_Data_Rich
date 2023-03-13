
%
%   Graph of Impulse Response Function 
%
%
%


case_data_rich = 3;  % Case A =1, Case B = 2,  Case C =3 

if case_data_rich == 1
     load irf_A_MP.txt
     n = 7
     AAA=irf_A_MP;
     
     figure(1)
     
end     

if case_data_rich == 2
     load irf_B_MP.txt
     n = 14
     AAA=irf_B_MP;
     
      figure(2)
end     


if case_data_rich == 3
     load irf_C_MP.txt
     n = 55;
     AAA=irf_C_MP;
     
      figure(3)
end     
     

AA=zeros(20,3*n);

title_name={ 'Output' ;    'Inflation' ;   'Real Wage' ;   'Investment' ;   'Consumption' ;           
           'Nominal Rate' ;    'Labor'  ;   'Output 2' ;  'Inflation 2' ;   'Real Wage 2' ;  
		   'Investment 2' ;  'Consumption 2' ;  'Nominal Rate 2' ;                      
           'Labor 2'  ;   'Output 3' ;  'Inflation 3' ;   'Real Wage 3' ;  
		   'Investment 3' ;      'Consumption 3' ;   'Nominal Rate 3' ;      
           'Labor 3'  ;
		   'Consumption 4' ; 'Consumption 5' ;'Consumption 6'; 'Consumption 7' ; 'Consumption 8' ; 
     'output 4' ; 'output 5' ; 'output 6' ;
     'output 7' ; 'output 8'  ; 'Labor 4' ; 'Labor 5' ; 'Labor 6' ; 'Interest Rate4'; 
      'Interest Rate 5'; 'Interest Rate 6' ; 'Interest Rate 7' ; 'Interest Rate 8' ;
      'Inflation 4' ; 'Inflation 5' ; 'Inflation 6' ; 'Inflation 7' ; 'Inflation 8' ; 
       'Inflation 9' ; 'Inflation 10'; 'Stock Price 1' ; 'Stock Price 2' ; 'Stock Price 3' ; 
	   'Commodity Price 1' ; 'Commodity Price 2' ; 'Money Stock 1' ; 'Money Stock 2' ;
	   'Exchange Rate' ; 'Real Exchange Rate'; };

for i=1:n
    AA(:,1+(i-1)*3:3+(i-1)*3)=AAA((i-1)*20+1:(i-1)*20+20,1:3);
end



   % figure(1)

    
for j=1:7
           
   subplot(3,3,j)
%for j=1+7*(k-1):k*7  
%   subplot(3,3,j-7*(k-1))
   hold on
 
  for i=0:0.01:1
      plot(1:20,AA(:,3+3*(j-1))-i*(AA(:,3+3*(j-1))-AA(:,2+3*(j-1))),'Color',[0.753 0.753 0.753],'LineWidth',2)
  end
      plot(1:20,AA(:,1+3*(j-1)),'b','LineWidth',2)
      
      hold off
      title( title_name(j),'FontSize',14  )  
end






