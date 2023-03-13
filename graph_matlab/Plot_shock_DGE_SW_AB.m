%
%
%   Figure 2 
%
%


case_A = 1 %  on-> 1,  off-> 0
case_B = 1 %  on-> 1,  off-> 0

n = 7;  % number of shocks 
T = 59;  % number of periods

load shock_sample_A.txt;
AAA=shock_sample_A;
AA=zeros(T,3*n);

case_data_rich == 3

if case_data_rich == 2
    load shock_sample_B.txt;
    BBB=shock_sample_B;
    BB=zeros(T,3*n);
elseif case_data_rich == 3
    load  shock_sample_C.txt;
    BBB=shock_sample_C;
    BB=zeros(T,3*n);
end     




 figure(3)

titlename={'Preference Shock ' ;
           'Investment Shock ' ; 
           'Equity Premium Shock ' ;
            'Labor Supply Shock ' ;            
           'Productivity Shock ' ;                          
           'Government Spending Shock ' ;
           'Monetary Policy Shock ' };
       
for i=1:n
    AA(:,1+(i-1)*3:3+(i-1)*3)=AAA((i-1)*T+1:(i-1)*T+T,1:3);

    BB(:,1+(i-1)*3:3+(i-1)*3)=BBB((i-1)*T+1:(i-1)*T+T,1:3);
end


    
for j=1:n  
   subplot(3,3,j)
   hold on
 
   if case_A == 1
%       for i=0:0.01:1
%          plot(1981.20:0.25:1995.75,AA(:,3+3*(j-1))-i*(AA(:,3+3*(j-1))-AA(:,2+3*(j-1))),'c','LineWidth',2)    
%       end

          h=area(1981.20:0.25:1995.75,[ AA(:,2+3*(j-1))  (AA(:,3+3*(j-1))-AA(:,2+3*(j-1))) ] );
          set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
          set(h(2),'FaceColor',[0.5 1 1])      
          set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value
   
         plot(1981.20:0.25:1995.75,AA(:,1+3*(j-1)),'b','LineWidth',2)
   end
         
         
  if case_B == 1          
%      for i=0:0.01:1
%         plot(1981.20:0.25:1995.75,BB(:,3+3*(j-1))-i*(BB(:,3+3*(j-1))-BB(:,2+3*(j-1))),'Color',[1 0.4 0.6],'LineWidth',2)    
%      end

       h=area(1981.20:0.25:1995.75,[ BB(:,2+3*(j-1))  (BB(:,3+3*(j-1))-BB(:,2+3*(j-1))) ] );
          set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
          set(h(2),'FaceColor',[1 0.4 0.6])   % Pink     
          set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value
         
        plot(1981.20:0.25:1995.75,BB(:,1+3*(j-1)),'r','LineWidth',2)
  end
  
        
        title( titlename(j),'FontSize',10 )  
            %title('éYèo','FontSize',14)
  xlim([min(1981.25),max(1996.0)])
  set(gca,'fontsize',12);
  hold off
      
end
 
      figure(4)
     
for j=1:n  
   subplot(3,3,j)
   hold on
   
   
          
  if case_B == 1          
%      for i=0:0.01:1
%         plot(1981.20:0.25:1995.75,BB(:,3+3*(j-1))-i*(BB(:,3+3*(j-1))-BB(:,2+3*(j-1))),'Color',[1 0.4 0.6],'LineWidth',2)    
%      end

       h=area(1981.20:0.25:1995.75,[ BB(:,2+3*(j-1))  (BB(:,3+3*(j-1))-BB(:,2+3*(j-1))) ] );
          set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
          set(h(2),'FaceColor',[1 0.4 0.6])   % Pink     
          set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value
         
        plot(1981.20:0.25:1995.75,BB(:,1+3*(j-1)),'r','LineWidth',2)
  end
   
   
   
 
   if case_A == 1
%       for i=0:0.01:1
%          plot(1981.20:0.25:1995.75,AA(:,3+3*(j-1))-i*(AA(:,3+3*(j-1))-AA(:,2+3*(j-1))),'c','LineWidth',2)    
%       end

        h=area(1981.20:0.25:1995.75,[ AA(:,2+3*(j-1))  (AA(:,3+3*(j-1))-AA(:,2+3*(j-1))) ] );
          set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
          set(h(2),'FaceColor',[0.5 1 1])      
          set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value        
   
         plot(1981.20:0.25:1995.75,AA(:,1+3*(j-1)),'b','LineWidth',2)
   end
         
  
  
        
        title( titlename(j),'FontSize',10 )  
            %title('éYèo','FontSize',14)
  xlim([min(1981.25),max(1996.0)])
  set(gca,'fontsize',12);
  hold off
      
end
