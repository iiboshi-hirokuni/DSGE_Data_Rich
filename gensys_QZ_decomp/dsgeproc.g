/**********************************************************
**
**  Procedure to compute the steady state
**  for the DSGE specification 
**
**********************************************************/

proc (2) = dsgesst(para, spec);
/* INPUT:
**   spec: 1 and 2
**     para is 8*1 vector of structural parameters (spec 1,2)
**            11*1 vector of structural parameters (spec 3)
**
**   spec is scalar that contains model specification
**     1: RBC
**     2: LBD
**     3: LBD and effort
**
** OUTPUT: sst is 5*1 vector containing steady state values
**         yst, cst, ist, kst, rst
**         valid = 1: steady state values are valid
**         valid = 0: steady state values are invalid
*/
local alp, bet, gam, del, nu, rho, mu , phi, chi, gst,
      valid, sst, rst, kst, yst, ist, cst,
      k_y, im_y, ih_y, c_y, shm;


/* Initialize Output Variables
*/
sst = zeros(5,1);
valid = 1;

/* assign names to parameters
*/
if (spec == 1) OR (spec == 2);
   alp = para[1];
   bet = para[2];
   gam = para[3];
   del = para[4];
   nu  = para[5];
   rho = para[6];
   mu  = para[7];
   phi = para[8];
elseif (spec == 3);
   alp = para[1];
   bet = para[2];
   gam = para[3];
   del = para[4];
   nu  = para[5];
   rho = para[6];
   mu  = para[7];
   phi = para[8];
   chi = para[9];
else;
   valid = 0;
   goto eosstdsge;
endif;

gst = exp(gam);

rst = gst/bet - 1 + del;
if rst < 1E-7;
   valid = 0;
   goto eosstdsge;
endif;

if spec == 1;
   
   kst = ( (1-alp)/Rst )^(1/alp) *gst; /* *lm */
   yst = kst*rst/( (1-alp)*gst );
   ist = kst*(1 - (1-del)/gst );
   cst = yst - ist;
   sst = yst | cst | ist | kst | rst;

elseif (spec == 2) OR (spec == 3);

   kst = ( (1-alp)/Rst )^(1/alp)*gst; /* *lm^( 1+mu/(1-phi) ) */
   yst = kst*rst/( (1-alp)*gst );
   ist = kst*(1 - (1-del)/gst );
   cst = yst - ist;
   sst = yst | cst | ist | kst | rst;

endif;
eosstdsge:

retp(sst,valid);
endp;


/**********************************************************
**
**  Procedure to solve the DSGE model with GENSYS
**
**********************************************************/


proc (4) = dsgesolv(para,spec);
local alp, bet, gam, del, nu, rho, mu, phi, chi, gst,
      sst, valid, yst, cst, ist, kst, rst,
      et, sc, sk, kap, xi, 
      GAM0, GAM1, C, PSI, PPI, 
      T1, TC, T0, TY, M, TZ, TETA, GEV, RC;

/* Define the model parameters
*/

if (spec == 1) OR (spec == 2);

   alp = para[1];
   bet = para[2];
   gam = para[3];
   del = para[4];
   nu  = para[5];
   rho = para[6];
   mu  = para[7];
   phi = para[8];

elseif (spec == 3);

   alp = para[1];
   bet = para[2];
   gam = para[3];
   del = para[4];
   nu  = para[5];
   rho = para[6];
   mu  = para[7];
   phi = para[8];
   chi = para[9];
   
endif;

gst  = exp(gam);

/* Determine the steady state of the stochastically
** detrended CIA model
*/
if (spec == 1) OR (spec == 2);
   {sst, valid} = dsgesst( (alp|bet|gam|del|nu|rho|mu|phi), spec );
elseif (spec ==3);
   {sst, valid} = dsgesst( (alp|bet|gam|del|nu|rho|mu|phi|chi), spec );
endif;
   
if valid == 0;
   RC = -5;
   T1 = 0;
   TC = 0; 
   T0 = 0;
   goto eofproc;
endif;

yst = sst[1];
cst = sst[2];
ist = sst[3];
kst = sst[4];
rst = sst[5];

/* Define some additional constants
*/

et  = rst/(rst+1-del);
sc  = cst/yst;
sk  = kst/yst;

if (spec == 2);
    kap = mu*bet/(1-phi*bet+mu*bet); /* LBD is considered in labor supply and inv decision */
  /*  kap = alp*mu*bet/(1-phi*bet+alp*mu*bet); */ /* LBD is only considered in labor supply decision */
    xi  = 0;
elseif (spec == 3);
    kap = mu*bet/(1-phi*bet+mu*bet);
    xi  = chi*(1-phi*bet)/( (1+chi)*(1-phi*bet+bet*mu) );

endif;

/* Define Matrices of canonical system  
*/

if spec == 1;
   /*  RBC Model
   **
   **  (1) Intratemporal constraint
   **  (2) Aggregate Resource Constraint
   **  (3) preference shock
   **  (4) Euler Equation
   **
   ** [dl(t), dk(t+1), dc(t), db(t)]
   */

      GAM0 = ( (alp-1) - 1/nu ~ 0   ~ -1 ~ -1 )|
             ( alp            ~ -sk ~ -sc ~ 0 )|
             ( 0              ~ 0     ~ 0 ~  1)|
             ( et*alp ~ 0  ~  -1 ~ 0 );

      GAM1 = ( 0                   ~ alp-1 ~ 0 ~ 0)|
             ( 0 ~ alp-1 - sk*(1-del)/gst  ~ 0 ~ 0)|
             ( 0                   ~ 0     ~ 0 ~ rho)|
             ( 0        ~ alp*et ~ -1 ~ 0);

      C    = zeros(rows(GAM0),1);

      PSI  = ( 1-alp ~ 0 )|
             ( 1-alp + sk*(1-del)/gst ~ 0 )|
             ( 0     ~ 1 )|
             ( 0     ~ 0 );

      PPI  = ( 0 |
               0 |
               0 |
               1 );

elseif (spec == 2);
   /*   LBD
   **
   **  (1) Intratemporal optimality
   **  (2) Aggregate resource constraint
   **  (3) learning by doing
   **  (4) preference shock
   **  (5) Euler equation 1
   **  (6) Euler Equation 2
   **
   ** [dl(t), dk(t+1),  dc(t), db(t)] [dx(t+1), dd(t)]
   */
 

      GAM0 = ( ( (1-kap)*(alp-1)+kap*(mu-1) - 1/nu ) ~  0   ~ kap-1 ~ -1 ~ 0 ~ kap)|
             ( alp ~ -sk  ~ -sc  ~ 0 ~ 0 ~ 0 )|
             ( -mu ~ 0    ~ 0    ~ 0 ~ 1 ~ 0 )|
             ( 0   ~ 0    ~ 0    ~ 1 ~ 0 ~ 0 )|
             (et*alp ~ 0          ~ -1  ~ 0 ~ 0 ~ 0 )|
             (-mu*phi*bet/(1-phi*bet) - alp  ~  0  ~ 1 ~ 0 ~ 0 ~ -phi*bet/(1-phi*bet));
     

      GAM1 = ( 0 ~ (kap-1)*(1-alp)  ~  0  ~ 0 ~ (kap-1)*alp - kap*phi ~ 0 )|
             ( 0 ~ alp-1 - sk*(1-del)/gst ~  0  ~ 0  ~ - alp ~ 0 )|
             ( 0 ~ 0                      ~  0  ~ 0  ~ phi   ~ 0 )|
             ( 0 ~ 0                      ~  0  ~ rho~ 0     ~ 0 )|
             ( 0       ~ alp*et  ~ -1  ~ 0 ~  -et*alp ~  0)|
             ( 0       ~ (1-alp)  ~  0  ~ 0 ~ phi*bet/(1-phi*bet)*(phi-1) + alp-1 ~ - 1/(1-phi*bet));

      C    = zeros(rows(GAM0),1);

      PSI  = ( (1-kap)*(1-alp) ~ 0 )|
             ( 1-alp + sk*(1-del)/gst ~ 0 )|
             ( 0     ~ 0 )|
             ( 0     ~ 1 )|
             ( 0     ~ 0 )|
             ( 0     ~ 0 );

      PPI  = zeros(4,2) | eye(2) ;
            
elseif (spec == 3);

   /*   LBD and Effort
   **
   **  (1) Intratemporal optimality
   **  (2) Aggregate resource constraint
   **  (3) learning by doing
   **  (4) preference shock
   **  (5) Effort
   **  (6) Euler equation 1
   **  (7) Euler Equation 2
   **
   ** [dl(t), dk(t+1),  dc(t), db(t)] [dx(t+1), dd(t)] [e(t)]
   */
 
      GAM0 = ( ( (1-kap)*(alp-1)+kap*(mu-1) - (1-xi)/nu ) ~  0   ~ kap-1 ~ -1 ~ 0 ~ kap 
                 ~ -(1-phi*bet)/(1-phi*bet+bet*mu) +(1-kap)*(alp-1))|
             ( alp ~ -sk  ~ -sc  ~ 0 ~ 0 ~ 0 ~ alp)|
             ( -mu ~ 0    ~ 0    ~ 0 ~ 1 ~ 0 ~ 0  )|
             ( 0   ~ 0    ~ 0    ~ 1 ~ 0 ~ 0 ~ 0  )|
             ( chi*(alp-1) ~ 0  ~ -chi   ~ -chi ~ 0 ~ 0 ~ chi*(alp-1)-1 )|
             (et*alp ~ 0          ~ -1  ~ 0 ~ 0 ~ 0 ~ et*alp)|
             (-mu*phi*bet/(1-phi*bet) - alp  ~  0  ~ 1 ~ 0 ~ 0 ~ -phi*bet/(1-phi*bet)~ -alp);
     

      GAM1 = ( 0 ~ (kap-1)*(1-alp)  ~  0  ~ 0 ~ (kap-1)*alp - kap*phi ~ 0 ~ 0)|
             ( 0 ~ alp-1 - sk*(1-del)/gst ~  0  ~ 0  ~ - alp ~ 0 ~ 0)|
             ( 0 ~ 0                      ~  0  ~ 0  ~ phi   ~ 0 ~ 0)|
             ( 0 ~ 0                      ~  0  ~ rho~ 0     ~ 0 ~ 0)|
             ( 0 ~ chi*(alp-1)            ~  0  ~ 0  ~ chi*(-alp) ~ 0 ~ 0)|
             ( 0       ~ alp*et  ~ -1  ~ 0 ~  -et*alp ~  0 ~ 0)|
             ( 0       ~ (1-alp)  ~  0  ~ 0 ~ phi*bet/(1-phi*bet)*(phi-1) + alp-1 ~ - 1/(1-phi*bet) ~ 0);

      C    = zeros(rows(GAM0),1);

      PSI  = ( (1-kap)*(1-alp) ~ 0 )|
             ( 1-alp + sk*(1-del)/gst ~ 0 )|
             ( 0     ~ 0 )|
             ( 0     ~ 1 )|
             ( chi*(1-alp) ~ 0 )|
             ( 0     ~ 0 )|
             ( 0     ~ 0 );

      PPI  = zeros(5,2) | eye(2);

endif;

{T1,TC,T0,TY,M,TZ,TETA,GEV,RC} = gensys(GAM0, GAM1, C, PSI, PPI, 1, 1);

eofproc:

retp(T1,TC,T0,RC);
endp;
