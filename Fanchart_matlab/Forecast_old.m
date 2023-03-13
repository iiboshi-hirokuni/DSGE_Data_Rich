%=======================================================
%
%     Forecasting by a Bayesain VAR model 
%
%=======================================================

load ./data/Phip.mat
load ./data/Sigmap.mat
load ./data/XXact.mat
load ./data/YYact.mat

nsim = 3000;
nburn = 2000;
p = 2;
n = 4;


% load mcmc_sample.mat

[T,k]=size(XXact);  % number of sample 

h = 60;           % horizon of forecasting

YY_f = zeros(nsim-nburn,h,n);
XX_f = zeros(nsim-nburn,h,n*p+1);

%*********************************************
%
%  MCMCによるVARの予測値
%
%*********************************************


for j = 1:1:nsim-nburn

    XX_f(j,1,:) = XXact(T,:);   % 1期目のX_fの値の設定　
    
    for i =1:1: h
    
        Phi = reshape(Phip(j,:,:),n*p+1,n);    
   
        Sigma = reshape(Sigmap(j,:,:),n,n); 
        
        XX = reshape(XX_f(j,i,:),1,n*p+1) ;
   
        YY_f(j,i,:) = XX*Phi(:,:) + randn(1,n)*chol(Sigma) ;
 
        XX_f(j,i+1,1:n) = YY_f(j,i,:);
        if p > 1
            XX_f(j,i+1,n+1:n*p) = YY_f(j,i,:);
        end
        XX_f(j,i+1, n*p+1) =1;
 
  end
end

%*********************************************
%
%  ファンチャートの作成
%
%*********************************************

a  =0.90;   % (1-2*rate) percentage of Credibile interval of Parameters  
  rate = (1-a)/2;
    
m = 10;   %  ファンチャートの階級の数
b = a/2/m; % ファンチャートの各階級の幅 

     % ファンチャートの色の設定      
      col1 = 0.75; %    青味の調整　 
      col2 = 0.65; %    青味の調整　
      col3 = 0.2; %  col3< col1  黒・白色の調整　大きいと黒くなる

  
sort_forecast = zeros( nsim-nburn, h, n );


for k = 1:1:n 
  for i=1:1:h
     sort_forecast(:,i,k) = sort(YY_f(:,i,k),1);
  end
end  
  
for i =1:n
  %  Graph of Forecast 
  YY_mean = mean(YY_f(:,:,i),1)';
  
  
  YY_band = zeros(h,2*m);
  sort_forecast_1 = squeeze(sort_forecast(:,:,i));
  YY_median = sort_forecast_1(round((nsim-nburn)*0.5),:)'; 
  
  for j = 1:m 
      
      YY_band(:,j) = sort_forecast_1(round((nsim-nburn)*(rate+b*(j-1))),:)'; 
      YY_band(:,2*m-j+1) = sort_forecast_1(round((nsim-nburn)*(1-(rate+b*(j-1)))),:)';
      
%       rate+b*(j-1)
%       1-(rate+b*(j-1))
  end
  
  YY_f_save = YY_band(:,1);
  for j = 1:m-1
    YY_f_save = [ YY_f_save YY_band(:,j+1)-YY_band(:,j) ];
  end  
    
    YY_f_save = [ YY_f_save YY_median-YY_band(:,m) YY_band(:,m+1)-YY_median ];
    
   for j = m:-1:2
     YY_f_save = [ YY_f_save YY_band(:,2*m+2-j)-YY_band(:,2*m+1-j) ];
   end
 
  if i == 1  
        figure('Name','Forecast of GDP','Color','w');
  elseif i==2
          figure('Name','Forecast of Inflation','Color','w');
  elseif i==3
          figure('Name','Forecast of Interest Rate','Color','w');              
  elseif i==4
          figure('Name','Forecast of Real Money Blance','Color','w');
  end
          
     tt = size([YYact(:,1);YY_mean],1);  
     ti=linspace(1964-p,1964-p+tt*0.25,tt)';
  
 
 hold on 
  
 hh = area(ti, [YYact(:,i) zeros(size(YYact(:,1),1),2*m);...
               YY_f_save ] ) ;
      
 
           
      for j = 1:m    
         set(hh(j),'FaceColor',[col1*(1-(j-1)*1/m)+(1-col1) col2*(1-(j-1)*1/m)+(1-col2) col3*(1-(j-1)*1/m)+(1-col3) ])
    
         set(hh(2*(m+1)-j),'FaceColor',[col1*(1-(j)*1/m)+(1-col1) col2*(1-(j)*1/m)+(1-col2) col3*(1-(j)*1/m)+(1-col3) ])

      end
         set(hh(m+1),'FaceColor',[1-col1 1-col2 1-col3 ])
         set(hh(1),'FaceColor',[1 1 1 ])
      
      %  罫線の設定 
      set(hh,'LineStyle','none') % Set all to same value
%       set(hh,'LineStyle','none','LineWidth',0.1) % Set all to same value
  
      plot(ti, [YYact(:,i); YY_median ],'LineStyle','-','Color','b', 'LineWidth',2.5); 
 
      hold off
end
  
 