/* filename: loaddata.g
** description: reads time series data from ASCII data set into
** Gauss variable SERIES
** created:  12/10/01
** Modifications: 10/12/4 change for BGGGK data(datasel=9(else) ), change data_us.txt for BGGGK data
**
*/

if datasel == 1;

   /* USE SIMULATED DATA
   */

   lsimdata  = datapath $+ "\\" $+ lmodel $+ "sim" $+ dataselstr;
   open fhsimdata = ^lsimdata for read;
   series_YT = readr(fhsimdata,Tsim[datasel]);
   closeall fhsimdata;   
   
   ti = seqa(1,1,Tsim[datasel]);

elseif datasel == 2;  

    titlestr =
     "output_gap " $|"Inflation" $|"Wage" $|"Investment" $|"Consumption" $|"Interest_Rate" $|"Labor"$|
     "output_gap_2 " $|"Inflation_2" $|"Wage_2" $|"Investment_2" $|"Consumption_2" $|"Interest_Rate_2" $|"Labor_2"$|
     "output_gap_3 " $|"Inflation_3" $|"Wage_3" $|"Investment_3" $|"Consumption_3" $|"Interest_Rate_3" $|"Labor_3"$|
     "Consumption_4" $| "Consumption_5" $|"Consumption_6"$| "Consumption_7" $| "Consumption_8" $| 
     "output_4" $| "output_5" $| "output_6" $|
     "output_7" $| "output_7"  $| "Labor_4" $| "labor_5" $| "labor_6" $| "Interest Rate_4"$| 
      "Interest_Rate_5"$| "Interest Rate_6" $| "Interest Rate_7" $| "Interest Rate_8" $|
      "Inflation_4" $| "Inflation_5" $| "Inflation_6" $| "Inflation_7" $| "Inflation_8" $| 
       "Inflation_9" $| "Inflation_10" $| "Stock_Price_1" $| "Stock_Price_2" $| "Stock_Price_3" $| 
	   "Commodity_Price_1" $| "Commodity_Price_2" $| "Money_Stock_1" $| "Money_Stock_2" $|
	   "Exchange_Rate" $| "Real_Exchange_Rate";

       titlestr_save ={};
	   
 load path   = ^datapath series_yt[77,55+1]= "data_rich_80-98.csv";

   ti          = seqa(1981.25,0.25,77);

      series_yt = real(series_yt[6:65,2:55+1]);
    
   
   if  mhrun == "3" ;    
       load path = ^priopath case_c_zz_1[55+1,12+1] = "case_a_zz.csv";
   elseif mhrun == "4";
       load path = ^priopath case_c_zz_1[55+1,12+1] = "case_b_zz.csv";
   elseif (mhrun == "5") or (mhrun == "6" );
       load path = ^priopath case_c_zz_1[55+1,12+1] = "case_c_zz.csv";
   endif; 
 
         case_c_zz_1 = real(case_c_zz_1[2:55+1,1+1:12+1]);

         series_yt_1 = series_yt;  

         index_type = 1;  case_c_zz={};  series_yt={}; 

         do until index_type > rows(case_c_zz_1);

            if case_c_zz_1[index_type,4] ne 0;
               
               case_c_zz = case_c_zz|case_c_zz_1[index_type,.]; 

               series_yt = series_yt~series_yt_1[.,index_type]; 

               titlestr_save = titlestr_save $|titlestr[index_type];	   
 
            endif;

           index_type = index_type + 1;

          endo;              
   
 if  mhrun == "6" ;    
       
         index_type = 1;  factor_yt={};  table_accel={};

         do until index_type > rows(case_c_zz);

            if case_c_zz[index_type,3] == 1;
               
               factor_yt = factor_yt~series_yt_1[.,index_type]; 
			   
			   table_accel = table_accel | index_type ;

            endif;

           index_type = index_type + 1;

          endo;              
		  
endif;	  
   
      
/****************************
 Drawing Graph of Data  
***************************/
                        
   graphset;
    begwind;
    //margin(0,0,0.5,0.5);
    window(3,3,0);           @@
    //fonts("microb");
    _ptitlht = 0.4;
    _paxht = 0.4;
    _pnumht = 0.3;
    //_pcolor = {4 4 4};
    //_pltype = {6 1 1};
    //_protate = 0; 
    //_plctrl = 0;
    //_plwidth = 5.5;
    //_psymsiz = 5;
    //_plegstr = "Mean\00090% HPD(U)\00090% HPD(L)";

       yy = series_yt ;
	   
       y_ind = 1;
       do until y_ind > 7;
       /* iterate over endogenous variables */
          
          title(titlestr_save[y_ind]);

          if (datasel ==2) or (datasel ==3);
              xtics(1981, 1995,2,4);
              x=seqa(0/4,1/4, rows(yy))+1981;
         
          endif;
          
          xy(x, yy[.,y_ind] );  

          nextwind;
                 
          y_ind = y_ind+1;

       endo;
       
    
       endwind;   
 
        
      
endif;

"data setting table" ;
print case_c_zz;  
"press any key!!";
wait;


