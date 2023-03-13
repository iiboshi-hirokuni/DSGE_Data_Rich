/* filename: loaddata.g
** description: reads time series data from ASCII data set into
** Gauss variable SERIES
** created:  12/10/01
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

elseif datasel == 2;   /*  Japan Data  */

    load path   = ^datapath series_yt[77,56]= "data_rich_80-98.csv";

      ti          = seqa(1981.25,0.25,77);

      series_yt = real(series_yt[6:65,2:56]);
    
       load path = ^priopath case_c_zz_1[55+1,11+1] = "case_sw_zz.csv";
 
         case_c_zz_1 = real(case_c_zz_1[2:55+1,1+1:11+1]);

         series_yt_1 = series_yt;  

         index_type = 1;  case_c_zz={};  series_yt={}; 

         do until index_type > rows(case_c_zz_1);

            if case_c_zz_1[index_type,4] ne 0;
               
               case_c_zz = case_c_zz|case_c_zz_1[index_type,.]; 

               series_yt = series_yt~series_yt_1[.,index_type]; 
 
            endif;

           index_type = index_type + 1;

          endo;              

          case_c_zz; "press any key" wait;
      


elseif datasel == 3;   /*  Sudo and Ueda data  */

    load path   = ^datapath series_yt[97,8]="sugo_ueda_data.txt";

      @from Mar 1980  to Dec 1995 @
                                                         @y	pi  w  inv  c	R  L @
                                                          @ n=76,  1980/3  -  1998/12 @
                                                          @Full: 141,   not 0% period:  117   @
      ti          = seqa(1980,0.25,76);
      series_yt = series_yt[2:61,2:8];

load path   = ^datapath series_yt1[77,8] = "data_rich_japan_B.txt";   // Data Rich B
    
       series_yt1 = series_yt1[2:61,2:8];


     load path   = ^datapath series_yt2[77,8] = "data_rich_japan_C_core.txt";  // Data Rich C
    
       series_yt2 = series_yt2[2:61,2:8];


     load path   = ^datapath series_yt3[77,39] = "data_rich_japan_C_no_core.txt";  // Data Rich C
    
       series_yt3 = series_yt3[2:61,2:39];


else; 

   /* USE ACTUAL DATA
     */

   load path   = ^datapath series_yt[231,8]= "data_us.csv";
   ti          = seqa(1947,0.5,230);
      series_yt = series_yt[2:231,2:8];
         
   if subT == 2;
      series_yt = series_yt[1:78,.];
      ti        = ti[1:78,.];
   elseif subT == 3;
      series_yt = series_yt[79:152,.];
      ti        = ti[79:152,.];
   elseif subT == 4;
      series_yt = series_yt[92:152,.];
      ti        = ti[92:152,.];
   endif;
         
endif;


/****************************
 Drawing Graph of Data  
***************************/

nobs  = rows(series_yT);  /* number of observations */
YY    = series_yT;   

   
   yy_m = meanc(yy);  
   yy = yy - yy_m';    /* demean */    


titlestr = "Output  " $|
           "Inflation " $|         
           "Wage " $|          
           "Investment " $|
           "Consumption " $|
           "Nominal Interest Rate " $|
           "Labor "      ;    @@

graphset;
    begwind;
    margin(0,0,0.2,0.2);
    window(3,3,0);           @@
    fonts("microb");
    _ptitlht = 0.4;
    _paxht = 0.4;
    _pnumht = 0.3;
    @_pcolor = {4 4 4};@
    _pltype = {6 1 1};
    _protate = 0; 
    _plctrl = 0;
    _plwidth = 5.5;
    _psymsiz = 5;
    _plegstr = "Mean\00090% HPD(U)\00090% HPD(L)";

       
       y_ind = 1;
       do until y_ind > 7;
       /* iterate over endogenous variables */
          
          title(titlestr[y_ind]);

          if (datasel ==2) or (datasel ==3);
              xtics(1981, 1998,2,4);
              x=seqa(0/4,1/4, rows(yy))+1980;
          else;
              xtics(1947, 2005,2,4);
              x=seqa(0/4,1/4, rows(yy))+1947;
          endif;
          
          xy(x, yy[.,y_ind] );

          nextwind;
                 
          y_ind = y_ind+1;
       endo;
       
    
    endwind;   



