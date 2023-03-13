Solving Linear Rational Expectations Models 
            Using GAUSS - GENSYS
===========================================

A translation of Chris Sims' program from
MATLAB into GAUSS by Frank Schorfheide.
Email: schorf@ssc.upenn.edu
Web:   www.econ.upenn.edu/~schorf

The procedure uses the QZ decomposition
provided by Paul Soderlind
Email: Paul.Soderlind@hhs.se

The program has only been tested on the 
Learning-by-Doing Models, see
Chang, Gomes, Schorfheide "Learning-by-Doing
as a Propagation Mechanism".

1) Copy files in the following directories:

C:\gauss        : dforrt.dll
C:\gauss\dlib   : PsDgees.dll, PsDtgsen.dll, PsZgees.dll, PsZtgsen.dll
C:\gauss\usersrc: gensys.src

C:\projects\active\poleval\gauss:    dsgeproc.g, dsgetest.g

C:\gauss9.0\user\ dsgesel.src
C:\gauss9.0\user\matop.src

2) At the GAUSS-prompt type:
   lib user c:\gauss\usersrc\gensys.src

   lib user c:\gauss9.0\user\dsgesel.src

   lib user c:\gauss9.0\user\matop.src
   
   
3) test  dsgetest.g 


 
                