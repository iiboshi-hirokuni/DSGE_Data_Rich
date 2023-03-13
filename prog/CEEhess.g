/* filename:    CEE-hess.g
** description: The program computes the hessian at the posterior mode
** created:     09/20/2010
*/

new;
closeall;
library user, pgraph;
dlibrary -a PszTgSen,PsdTgSen;

format /mb1 /ros 16,8;
cls;


/* Global defaults 
*/
  __output = 1;
  _fcmptol = 1E-10;

#include c:\dge-data-rich\pathspec.g      @@
#include c:\dge-data-rich\ceespec.g      @@

/******************************************************** 
** Import data
** series (nobs,3)
*/
#include c:\dge_sw\loaddata.g       @@

nobs  = rows(series_yT);  /* number of observations */
YY    = series_yT;   

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
** Load Posterior Mode Estimates
*/


lpara = postpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "pm";
open fhpara = ^lpara for read;
para = readr(fhpara,npara);
closeall fhpara;


para_names={h  sigma_c  sigma_L  phi  phi0  psi gam_p 
            gam_w  xi_p   xi_w  rho_m 
            mu_pi  mu_y  rho_z  rho_c  rho_g  rho_L  rho_i 
            e_c  e_inv  e_q   e_L   e_w e_z  e_p  e_g  e_m };




"Parameter   | Estimate ";
outmat = para_names'~para;
let mask[1,2] = 0 1;
let fmt[2,3] =
   "-*.*s " 10 4
   "*.*lf " 10 4;
d = printfm(outmat,(0 ~ 1),fmt);

"";
"Prior*Likelihood at Mode";
fcn(para);
"Press Key to Continue";
k = keyw;
cls;

/* Compute Hessian, element by element, fine tune with dxscale
*/
comphess:

ndx = 6;
dx =  exp(-seqa(6,2,ndx));
hessian  = zeros( npara, npara );
gradx    = zeros(ndx,1);
grady    = zeros(ndx,1);
gradxy   = zeros(ndx,1);
hessdiag = zeros(ndx,1);
dxscale  = ones(npara,1);


/* Compute Diagonal elements first
*/
seli = 1;
do until seli > npara;
     
     locate 1,1;
     "Hessian Element    (" seli seli ")";
     i=1;
     do until i > ndx;
      paradx = para;
      parady = para;
      paradx[seli] = paradx[seli] + dx[i]*dxscale[seli];
      parady[seli] = parady[seli] - dx[i]*dxscale[seli];
      paradxdy = paradx;
      paradxdy[seli] = paradxdy[seli] - dx[i]*dxscale[seli];
      fx  = fcn(para);
      fdx = fcn(paradx);
      fdy = fcn(parady);
      fdxdy = fcn(paradxdy);
      gradx[i] = -( fx - fdx )/ (dx[i]*dxscale[seli]);
      grady[i] = ( fx - fdy )/ (dx[i]*dxscale[seli]);
      gradxy[i] = -(fx -fdxdy)/ sqrt( (dx[i]*dxscale[seli])^2 + (dx[i]*dxscale[seli])^2 );
      hessdiag[i] = -( 2*fx - fdx - fdy)/(dx[i]*dxscale[seli])^2; 
      hessdiag[i] = -( fx - fdx - fdy + fdxdy )/(dx[i]*dx[i]*dxscale[seli]*dxscale[seli]);
      i = i+1;
     endo;
     "Values";
     -hessdiag';
     hessian[seli,seli] = -0.5*(hessdiag[3]+hessdiag[4]);
     locate 6,1;
     "Value Used:" hessian[seli,seli];
     /* k = keyw; */
   seli = seli+1;
endo;

/* Now compute off-diagonal elements
** Make sure that correlations are between -1 and 1
** errorij contains the index of elements that are invalid
*/
errorij = 0~0~0;

seli = 1;
do until seli == npara;
   selj = seli+1;
   do until selj > npara;
     
     locate 1,1;
     "Hessian Element    (" seli selj ")";
     i=1;
     do until i > ndx;
      paradx = para;
      parady = para;
      paradx[seli] = paradx[seli] + dx[i]*dxscale[seli];
      parady[selj] = parady[selj] - dx[i]*dxscale[selj];
      paradxdy = paradx;
      paradxdy[selj] = paradxdy[selj] - dx[i]*dxscale[selj];
      fx  = fcn(para);
      fdx = fcn(paradx);
      fdy = fcn(parady);
      fdxdy = fcn(paradxdy);
      gradx[i] = -( fx - fdx )/ (dx[i]*dxscale[seli]);
      grady[i] = ( fx - fdy )/ (dx[i]*dxscale[selj]);
      gradxy[i] = -(fx -fdxdy)/ sqrt( (dx[i]*dxscale[selj])^2 + (dx[i]*dxscale[seli])^2 );
      hessdiag[i] = -( 2*fx - fdx - fdy)/(dx[i]*dxscale[seli])^2; 
      hessdiag[i] = -( fx - fdx - fdy + fdxdy )/(dx[i]*dx[i]*dxscale[seli]*dxscale[selj]);
      i = i+1;
     endo;
     "Values";
     -hessdiag';

     hessian[seli,selj] = -0.5*(hessdiag[3]+hessdiag[4]);
     
     if ( hessian[seli,seli] == 0 ) or (hessian[selj,selj] == 0);
        corrij = 0;
     else;
        corrij = hessian[seli,selj]/sqrt(hessian[seli,seli]*hessian[selj,selj]);
     endif;

     if (corrij < -0.98) OR (corrij > 0.98);
        hessian[seli,selj]=0.9*sqrt(hessian[seli,seli]*hessian[selj,selj]);
        errorij = errorij |(seli~selj~corrij);
     elseif (corrij > -0.005) AND (corrij < 0.005);
        hessian[seli,selj]=0;
     endif;   
     hessian[selj,seli] = hessian[seli,selj];

     locate 6,1;
     "Value Used: " hessian[seli,selj];
     "Correlation:" corrij;
     "Number of Errors:" rows(errorij)-1;
     selj=selj+1;
   endo;
   seli = seli+1;
endo;

cls;
"Errors"  errorij;

/* Initialize Output files
*/
ohess = postpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "hs";
create fhhess=^ohess with hess, npara, 8;
wr    = writer(fhhess,hessian[1:npara,1:npara]);
closeall fhhess;
   

/*******************************************************************************
*/

/* Load Hessian, compute penalty
*/
evalhess:

lhess = postpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "hs";
open fhhess = ^lhess for read;
HHm   = readr( fhhess,npara );
closeall fhhess;

@rankHHM = sumc(para_maskinv[.,mspec]);@  
rankHHM = sumc(ones(npara,1) );

/* Create Inverse by Singular Value Decomposition
*/
{u , s, v} = svd1(HHM);
invHHMdet = 1;

i = 1;
do until i > npara;
   if i > rankHHM;
      s[i,i] = 0;
   else;
      s[i,i]    = 1/s[i,i];
      invHHMdet = invHHMdet*s[i,i];
   endif;
   i = i+1;
endo;

invHHM  = u*s*u';
sigmult = u*sqrt(s);

"Determinant of minus Hessian";
invHHMdet;
"sqrt(Diagonal of Inverse Hessian)";
sqrt(diag(invHHM));

"Post Mode Penalty";
penalt = (rankHHM/2)*ln(2*pi) + 0.5*ln(invHHMdet);
penalt;

/* Initialize Output files
*/
omult = postpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "mu";
create fhmult=^omult with MULT, npara, 8;
wr = writer(fhmult,sigmult);

closeall fhmult;
end;


/****************************************************/
/*                 PROCEDURES                       */
/****************************************************/

proc (1) = fcn(para);
local lnpY, lnprio, obsmean, obsvar, retcode, modelpara,shock;

modelpara = para;

@ modelpara = modelpara.*para_maskinv[.,mspec] + para_fix[.,mspec].*para_mask[.,mspec]; @

{lnpY,retcode,obsmean,obsvar,shock} = evaldsge(modelpara,mspec,YY); 

/* Evaluate the Prior distribution
*/
lnprio = priodens(modelpara, pmean, pstdd, pshape);

retp(real(lnpY+lnprio));  
endp;

#include c:\dge_sw\CEEmod.src;    
#include c:\dge_sw\CEElh.src;      
