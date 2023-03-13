/* filename:    CEEpm.g
** description: The program maximizes the posterior density of the 
**              DSGE model 
** created:     10/05/2010
*/

new;
 library optmum, pgraph, user; 

 
 dlibrary -a PszTgSen,PsdTgSen;
cls;

format /mat /on /mb1 /ros 16,8;


/* Global defaults
*/
  __output = 1;
  _fcmptol = 1E-10;

#include c:\dge-data-rich\pathspec.g
#include c:\dge-data-rich\ceespec.g
#include c:\dge-data-rich\ceemod.src;
#include c:\dge-data-rich\ceelh.src;

/* 
optimization method
*/

opt = 2;   // 1:csminwel  2:optmum

/******************************************************** 
** Import data
** series (nobs,7)
*/
#include c:\dge-data-rich\loaddata.g

nobs  = rows(series_yT);  /* number of observations */
YY    = series_yT;  // /100 ; 

   
   yy_m = meanc(yy);  
   yy = yy - yy_m';    /* demean */    



/******************************************************** 
** Import the information for the prior density
*/
priorfile = lmodel $+ lprior $+ dataselstr $+ subTstr $+ "prt.out";
load path = ^priopath prior[npara,3] = ^priorfile;

pmean  = prior[.,1];
pstdd  = prior[.,2];
pshape = prior[.,3];

/*
pshape = pshape.*para_maskinv[.,mspec];
*/
/******************************************************** 
** Starting Values for Maximization
*/


h       = 0.5;
sigma_c = 1;
sigma_L = 1;
phi     = 4;
phi0    = 2;
psi     = 0.15;
gam_p   = 0.3;
gam_w   = 0.7;
xi_p    = 0.88; 
xi_w    = 0.62;
rho_m   = 0.8;
mu_pi   = 1.7;
mu_y    = 0.12;
rho_z  =  0.6;
rho_c   =  0.5;
rho_g   =  0.9;
rho_L   =  0.9; 
rho_i   =  0.75;
e_c     =  2;
e_inv   = 1.8;
e_q     = 6;
e_L     = 1;
e_w     = 10;
e_z     = 1;
e_p     = 10;
e_g     = 1;
e_m     = 0.2;     

para    =  h | sigma_c | sigma_L | phi | phi0 | psi| gam_p |
            gam_w | xi_p |  xi_w | rho_m |
            mu_pi | mu_y | rho_z | rho_c | rho_g | rho_L | rho_i |
            e_c | e_inv | e_q |  e_L |  e_w | e_z | e_p | e_g | e_m ;



npara = rows(para);

para_names={h  sigma_c  sigma_L  phi  phi0  psi gam_p 
            gam_w  xi_p   xi_w  rho_m 
            mu_pi  mu_y  rho_z  rho_c  rho_g  rho_L  rho_i 
            e_c  e_inv  e_q    e_L   e_w  e_z  e_p e_g  e_m };
 

/* Load transformation scheme for parameters
*/
transfile = "trans.out";
load path=^priopath _trspec[npara,4] = ^transfile;


/******************************************************** 
**      Calculate posterior density at starting values
**
*/

cls;
"Prior Density at Starting Value:";  priodens(para,pmean,pstdd,pshape); 

{ lnpY,retcode,obsmean,obsvar, shock } = evaldsge(para,mspec,YY);

/* lnpy = -objfcn(invtrans(para)); */
"Posterior at Starting Value:    " lnpy;
"Press key to continue";
w = keyw;
cls;

/********************************************************
**       Maximize the posterior density
**
*/
cc = 0;
x0 = invtrans(para) ;  


H0= eye(npara)*1E-1;
nit = 1000;
crit= 1E-4;  //-5

/* Use personal optimization procedure
*/

if  opt == 1; 

 {fh,xh,g,H,itct,fcount,retcode} = csminwel(x0,H0,crit,nit,&objfcn);
 
elseif opt == 2;
 
  _opalgr  = 2;       // algorithm; 2:BFGS,  
  _opmiter = 1000;     // maximum number of iteration 
 
  {xh, fout, g, cout}=optmum(&objfcn, x0);

endif;
 

/********************************************************
** Save parameter estimates
*/

paraest = trans(xh);

oparaest         = postpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "pm";
create fhparaest = ^oparaest with PARAEST, 1, 8;
wr               = writer(fhparaest,paraest);
closeall fhparaest;   


/********************************************************
**       Estimation Output
**
*/
{lnpy, retcode, obsmean, obsvar,shock} = evaldsge( paraest,mspec,YY);

lnprio  = priodens( paraest, pmean, pstdd, pshape);

/* Compute gradient in terms of model parameters
*/
gradmod = g./ diag(gradp(&invtrans,paraest));

cls;

output file= c:\dge_sw\results\post_mode_para.txt reset;

output on;

"Model              " mspec;
"Dataset            " datasel;
"Subsample          " SubT;
"===================================";
"Posterior Mode     " (lnpy+lnprio);
"Likelihood at Mode " lnpy;
"Prior at Mode      " lnprio;
" ";
"Parameter | Estimate |  Start-V  | Prior Mean | Prior Stdd | Gradient";
outmat = para_names'~paraest~trans(x0)~pmean~pstdd~gradmod;
let mask[1,6] = 0 1 1 1 1 1;
let fmt[6,3] =
   "-*.*s " 8 10
   "*.*lf " 8 4
   "*.*lf " 8 4
   "*.*lf " 8 4
   "*.*lf " 8 4
   "*.*lf " 8 4;

d = printfm(outmat,(0~1~1~1~1~1),fmt);

if lnpy == -1E6 ;
   end;
endif;

output off;

/********************************************************
**       Compute Residuals
*/

/* Moments */

resid1 = YY - obsmean;
resid1 = resid1[2:rows(resid1),.];

resid1_m1 = meanc(resid1)';
resid1_m2 = meanc( (resid1 - resid1_m1)^2 )';
resid1_m3 = meanc( (resid1 - resid1_m1)^3 )';
resid1_m4 = meanc( (resid1 - resid1_m1)^4 )';

resid2 = (YY - obsmean)./sqrt(obsvar);
resid2 = resid2[2:rows(resid2),.];

resid2_m1 = meanc(resid2)';
resid2_m2 = meanc( (resid2 - resid2_m1)^2 )';
resid2_m3 = meanc( (resid2 - resid2_m1)^3 )';
resid2_m4 = meanc( (resid2 - resid2_m1)^4 )';

ti = ti[2:rows(ti)];

nresid = rows(resid1);

/* Autocorrelation */

resid1_ar = zeros(1,7);
resid2_ar = zeros(1,7);

ind = 1;
do until ind > 7;
   resid1_ar[1,ind]=inv(resid1[1:nresid-1,ind]'*resid1[1:nresid-1,ind])
                    *resid1[1:nresid-1,ind]'*resid1[2:nresid,ind];
   resid2_ar[1,ind]=inv(resid2[1:nresid-1,ind]'*resid2[1:nresid-1,ind])
                    *resid2[1:nresid-1,ind]'*resid2[2:nresid,ind];
   ind = ind+1;
endo;

output on;

"";
"";
"Raw Residuals";
"-------------";
"Mean:" resid1_m1;
"Stdd:" sqrt(resid1_m2);
"Skew:" resid1_m3./resid1_m2^(3/2);
"Kurt:" resid1_m4./(3*resid1_m2^2);
"Corr:" resid1_ar;
"";
"Standardized Residuals";
"----------------------";
"Mean:" resid2_m1;
"Stdd:" sqrt(resid2_m2);
"Skew:" resid2_m3./resid2_m2^(3/2);
"Kurt:" resid2_m4./(3*resid2_m2^2);
"Corr:" resid2_ar;
"";

output off;

/****************************
 Drawing Graph of Data  
***************************/
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
              xtics(1980, 2000,2,4);
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


    end;

/****************************************************/
/*                 PROCEDURES                       */
/****************************************************/

proc (1) = objfcn(para);
local lnpY, lnprio, obsmean, obsvar, shock, modelpara, retcode;

/* likelihood 
*/


modelpara = trans(para);



{lnpY,retcode,obsmean,obsvar,shock} = evaldsge(modelpara,mspec,YY); 



/* Evaluate the Prior distribution
*/

lnprio = priodens(modelpara, pmean, pstdd, pshape);


/*
locate 25,1;
"Likelihood" lnpY  "Return Code" retcode;
*/

retp(real(-lnpY-lnprio));  /* We minize the inverse of the likelihood fcn */
endp;


