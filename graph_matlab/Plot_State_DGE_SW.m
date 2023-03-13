%
%   Graph of State Variable 
%
%
%


case_data_rich = 3;  % Case A =1, Case B = 2,  Case C =3 

T = 59; % number of period

if case_data_rich == 1
     load state_sample_A.txt
     n = 7
     AAA=state_sample_A;
     
     figure(1)
     
end     

if case_data_rich == 2
     load state_sample_B.txt
     n = 14
     AAA=state_sample_B
     
      figure(2)
end     


if case_data_rich == 3
     load state_sample_C.txt
     n = 55;
    AAA=state_sample_C;
     
      figure(3)
end

AA=zeros(T,4*n);

for i=1:n
    AA(:,1+(i-1)*4:4+(i-1)*4)=AAA((i-1)*T+1:(i-1)*T+T,1:4);
end

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


    
for j=1:7
    
  
    
   subplot(3,3,j)
   hold on
   
   
%    for i=0:0.01:1
%        plot(1981.25:0.25:1995.75,AA(:,4+4*(j-1))-i*(AA(:,4+4*(j-1))-AA(:,3+4*(j-1))),'Color',[0.753 0.753 0.753],'LineWidth',2)    
%          %plot(1981.25:0.25:1995.75,AA(:,3),'r--')
%         %plot(1981.25:0.25:1995.75,AA(:,4),'r--')
%    end

       h=area(1981.20:0.25:1995.75,[ AA(:,3+4*(j-1))  (AA(:,4+4*(j-1))-AA(:,3+4*(j-1))) ] );
          set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
          set(h(2),'FaceColor',[0.753 0.753 0.753])  % Gray     
          set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value
          
         plot(1981.25:0.25:1995.75,AA(:,1+4*(j-1)),'b','LineWidth',2)
         plot(1981.25:0.25:1995.75,AA(:,2+4*(j-1)),'r--','LineWidth',2)

  title( title_name(j),'FontSize',14 )  
            %title('éYèo','FontSize',14)
  xlim([min(1981.25),max(1995.75)])
  set(gca,'fontsize',12);
  hold off
 
end

