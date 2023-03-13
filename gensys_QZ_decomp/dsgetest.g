/*************************************************************************
**  
**  This program solves the DSGE model via GENSYS
**  Author: Frank Schorfheide
**
**************************************************************************/

/* filename: lbdtest.g
** created:  5/1/01
** Modifications:
**
*/
library pgraph, user;
dlibrary -a PszTgSen,PsdTgSen;
#include c:\cee-dge\gensys_qz_decomp\dsgeproc.g

cls;

/* Define model parameters
*/
alp = 0.7267;
bet = 0.99;
gam = 0.0052;
del = 0.025;
nu  = 2;
rho = 0.926;
mu  = 0.33;
phi = 0.90;
chi = 2;

gst = exp(gam);

para = alp | bet | gam | del | nu |rho | mu | phi | chi ;
spec = 2;

{T1,TC,T0,RC} = dsgesolv(para,spec);

format /rdn 12,7;
"T1";
T1;

"T0";
T0;




