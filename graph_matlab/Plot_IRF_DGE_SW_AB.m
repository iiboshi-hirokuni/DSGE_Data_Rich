%
%
%   Figure 3 
%
%


%
%   Graph of Impulse Response Function 
%
%
%

case_A = 1 %  on-> 1,  off-> 0
case_B = 1 %  on-> 1,  off-> 0

case_data_rich = 3;  % Case A =1, Case B = 2,  Case C =3 


n = 7;  % number of shocks 
T = 20;  % number of horizon


load irf_A_MP.txt
      AAA=irf_A_MP;
      AA=zeros(T,3*n);

%       figure(1)
      
 if case_data_rich == 2    
     load irf_B_MP.txt
     BBB=irf_B_MP
     BB=zeros(T,3*n);
     figure(2)
 end
     
 if case_data_rich == 3    
     load irf_C_MP.txt
     BBB=irf_C_MP
     BB=zeros(T,3*n);
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
        
    BB(:,1+(i-1)*3:3+(i-1)*3)=BBB((i-1)*T+1:(i-1)*T+T,1:3);
end



   % figure(1)

    
for j=1:7
           
   subplot(3,3,j)
%for j=1+7*(k-1):k*7  
%   subplot(3,3,j-7*(k-1))
   hold on
   
   
  if case_A == 1
%       for i=0:0.01:1
%          plot(1:20,AA(:,3+3*(j-1))-i*(AA(:,3+3*(j-1))-AA(:,2+3*(j-1))),'c','LineWidth',2)    
%       end

         h=area(1:20,[ AA(:,2+3*(j-1))  (AA(:,3+3*(j-1))-AA(:,2+3*(j-1))) ] );
          set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
          set(h(2),'FaceColor',[0.5 1 1])      
          set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value


   
         plot(1:20,AA(:,1+3*(j-1)),'b','LineWidth',2)
  end
   
   
  if case_B == 1
%       for i=0:0.01:1
%          plot(1:20,BB(:,3+3*(j-1))-i*(BB(:,3+3*(j-1))-BB(:,2+3*(j-1))),'Color',[1 0.4 0.6],'LineWidth',2)    
%       end

        h=area(1:20,[ BB(:,2+3*(j-1))  (BB(:,3+3*(j-1))-BB(:,2+3*(j-1))) ] );
          set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
          set(h(2),'FaceColor',[1 0.4 0.6])   % Pink     
          set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value

   
         plot(1:20,BB(:,1+3*(j-1)),'r','LineWidth',2)
  end
  
  
      
      hold off
      title( title_name(j),'FontSize',14  )  
end



  figure(4)

    
for j=1:7
           
   subplot(3,3,j)
%for j=1+7*(k-1):k*7  
%   subplot(3,3,j-7*(k-1))
   hold on
   
     
  if case_B == 1
%       for i=0:0.01:1
%          plot(1:20,BB(:,3+3*(j-1))-i*(BB(:,3+3*(j-1))-BB(:,2+3*(j-1))),'Color',[1 0.4 0.6],'LineWidth',2)    
%       end

          h=area(1:20,[ BB(:,2+3*(j-1))  (BB(:,3+3*(j-1))-BB(:,2+3*(j-1))) ] );
          set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
          set(h(2),'FaceColor',[1 0.4 0.6])      
          set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value

   
         plot(1:20,BB(:,1+3*(j-1)),'r','LineWidth',2)
  end
  
  
  
  if case_A == 1
%       for i=0:0.01:1
%          plot(1:20,AA(:,3+3*(j-1))-i*(AA(:,3+3*(j-1))-AA(:,2+3*(j-1))),'c','LineWidth',2)    
%       end
      
          h=area(1:20,[ AA(:,2+3*(j-1))  (AA(:,3+3*(j-1))-AA(:,2+3*(j-1))) ] );
          set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
          set(h(2),'FaceColor',[0.5 1 1])      
          set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value
      
   
         plot(1:20,AA(:,1+3*(j-1)),'b','LineWidth',2)
  end
   
  
      
      hold off
      title( title_name(j),'FontSize',14  )  
end



