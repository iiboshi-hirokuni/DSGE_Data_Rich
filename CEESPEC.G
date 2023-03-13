/* Model Specification
*/

mspec  = 1;
lmodel = "m" $+ ftos(mspec,"%*.*lf",1,0) ;

mprior = 0;
lprior = ftos(mprior,"%*.*lf",1,0);


nvar   = 9;  nshock = 9;
npara  = 27;
ngrid  = 100;  /* number of grid points for density plots */


             
/* Selected data set

*/
datasel    = 2;    @3:sudo and ueda  2:Japan,  0:USA, 1:Simulation @
 dataselstr = ftos(datasel,"%*.*lf",1,0);


if datasel ne 1; 

   subT = 1;   
   
   subTstr = ftos(subT,"%*.*lf",1,0);
   
elseif datasel == 1; /* simulated data */

   Tsimtotal = 200;
   Tsimburn  = 200;
   Tsim = 200 | 500;
   subT = 1;   
   subTstr = ftos(subT,"%*.*lf",1,0);
     
endif;
