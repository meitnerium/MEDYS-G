

     ******************************************
     **    PROGRAM:              MCSCF       **
     **    PROGRAM VERSION:      5.5         **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************

 This program allows the csf mixing coefficient and orbital expansion coefficient
 optimization using the graphical unitary group approach and the exponential
 operator mcscf method.
 references:  r. shepard and j. simons, ' int. j. quantum chem. symp. 14, 211 (1980).
              r. shepard, i. shavitt, and j. simons, j. chem. phys. 76, 543 (1982).
              r. shepard in "ab initio methods in quantum chemistry ii" advances in chemical
                  physics 69, edited by k. p. lawley (wiley, new york, 1987) pp. 63-200.
 Original autor: Ron Shepard, ANL
 Later revisions: Michal Dallos, University Vienna

 This Version of Program MCSCF is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de



     ******************************************
     **    PROGRAM:              MCSCF       **
     **    PROGRAM VERSION:      5.4.0.2     **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************

 Workspace allocation information:
       131072000 of real*8 words ( 1000.00 MB) of work space has been allocated.

 user input information:

 ======== echo of the mcscf input ========
 ------------------------------------------------------------------------
  &input
   niter=100,
   nmiter=50,
   nciitr=300,
   tol(3)=1.e-4,
   tol(2)=1.e-4,
   tol(1)=1.e-8,
   NSTATE=0
   npath=1,3,9,10,13,17,19,21,-11,12, 2,
   ncoupl=5,
   tol(9)=1.e-3,
   FCIORB=  1,1,20,1,2,20
  &end
 ------------------------------------------------------------------------


 ***  Integral file informations  ***


 input integral file : /home/francois/workspace/medys-2016/BASIS/H2_R_1_40/1/WOR
 K/a

 Integral file header information:
                                                                                 
 aoints SIFS file created by argos.      franktoo          12:59:06.238 18-Jan-16

 Core type energy values:
 energy( 1)=  7.142857142857E-01, ietype=   -1,    core energy of type: Nuc.Rep.
 total ao core energy =    0.714285714


   ******  Basis set information:  ******

 Number of irreps:                  1
 Total number of basis functions:  10

 irrep no.              1
 irrep label             a
 no. of bas.fcions.    10
 inconsistent sifs parameter
 i,info(i),infoloc(i):     5  2700  2730
 bummer (warning):inconsistent sifs parameter 5


 ***  MCSCF optimization procedure parmeters:  ***


 maximum number of mcscf iterations:        niter=   100

 maximum number of psci micro-iterations:   nmiter=   50
 maximum r,s subspace dimension allowed:    nvrsmx=   30

 tol(1)=  1.0000E-08. . . . delta-emc convergence criterion.
 tol(2)=  1.0000E-04. . . . wnorm convergence criterion.
 tol(3)=  1.0000E-04. . . . knorm convergence criterion.
 tol(4)=  1.0000E-08. . . . apxde convergence criterion.
 tol(5)=  1.0000E-04. . . . small diagonal matrix element tolerance.
 tol(6)=  1.0000E-06. . . . minimum ci-psci residual norm.
 tol(7)=  1.0000E-05. . . . maximum ci-psci residual norm.
 tol(8)=  1.0000E+00. . . . maximum abs(k(xy)) allowed.
 tol(9)=  1.0000E-03. . . . wnorm coupling tolerance.
 tol(10)= 0.0000E+00. . . . maximum psci emergency shift parameter.
 tol(11)= 0.0000E+00. . . . minimum psci emergency shift parameter.
 tol(12)= 0.0000E+00. . . . increment of psci emergency shift parameter.


 *** State averaging informations: ***


 MCSCF calculation performed for  1 DRT.

 DRT  first state   no.of aver.states   weights
  1   ground state          1             1.000

 The number of hmc(*) eigenvalues and eigenvectors calculated each iteration per DRT:
 DRT.   no.of eigenv.(=ncol)
    1        2

 orbital coefficients are optimized for the ground state (nstate=0).

 Orbitals included in invariant subspaces:
   symmetry   orbital   mask
       1       1(  1)    20
       1       2(  2)    20

 npath(*) options:
  2:  orbital-state coupling terms will be included beginning on iteration ncoupl=  5
  3:  print intermediate timing information.
  9:  suppress the drt listing.
 10:  suppress the hmc(*) eigenvector listing.
 12:  diagonalize the hmc(*) matrix iteratively.
        nunitv= 1 nciitr=** mxvadd=20 nvcimx=20
       rtolci(*),wnorm=     1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 0.0000E+00
   noldv =   0
 13:  get initial orbitals from the formatted file, mocoef.
 17:  print the final natural orbitals and occupations.
 19:  transform the virtual orbitals to diagonalize qvv(*).
 21:  write out the one- and two- electron density for further use (files:mcd1fl, mcd2fl).


   ******  DRT info section  ******


 Informations for the DRT no.  1

 DRT file header:
  title                                                                          
 Molecular symmetry group:    a  
 Total number of electrons:    2
 Spin multiplicity:            1
 Number of active orbitals:    2
 Number of active electrons:   2
 Total number of CSFs:         3
 

 faar:   0 active-active rotations allowed out of:   1 possible.


 Number of active-double rotations:         0
 Number of active-active rotations:         0
 Number of double-virtual rotations:        0
 Number of active-virtual rotations:       16
 lenbfsdef=                131071  lenbfs=                    64
  number of integrals per class 1:11 (cf adda 
 class  1 (pq|rs):         #           6
 class  2 (pq|ri):         #           0
 class  3 (pq|ia):         #           0
 class  4 (pi|qa):         #           0
 class  5 (pq|ra):         #          48
 class  6 (pq|ij)/(pi|qj): #           0
 class  7 (pq|ab):         #         108
 class  8 (pa|qb):         #         192
 class  9 p(bp,ai)         #           0
 class 10p(ai,jp):        #           0
 class 11p(ai,bj):        #           0

 Size of orbital-Hessian matrix B:                      192
 Size of the orbital-state Hessian matrix C:             48
 Total size of the state Hessian matrix M:                0
 Size of HESSIAN-matrix for quadratic conv.:            240


 Source of the initial MO coeficients:

 Input MO coefficient file: /home/francois/workspace/medys-2016/BASIS/H2_R_1_40/1/WORK/m
 

               starting mcscf iteration...   1

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    10, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 131071506

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130905714
 address segment size,           sizesg = 130754228
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:         74 transformed 1/r12    array elements were written in       1 records.


 mosort: allocated sort2 space, avc2is=   130933359 available sort2 space, avcisx=   130933611

 trial vectors are generated internally.

 trial vector  1 is unit matrix column     1
 ciiter=   2 noldhv=  2 noldv=  2

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*       -1.1363850291       -1.8506707433        0.0000000000        0.0000010000
    2         0.0091880711       -0.7050976432        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.490420079330226E-003
 Total number of micro iterations:    4

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.74216745 pnorm= 0.0000E+00 rznorm= 1.9086E-15 rpnorm= 0.0000E+00 noldr=  4 nnewr=  4 nolds=  0 nnews=  0
 

 qvv(*) eigenvalues. symmetry block  1
     1.540928    2.616720    3.914762    3.914762    5.405447    5.856738    5.856738    9.056682

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    1 emc=     -1.1363850291 demc= 1.1364E+00 wnorm= 2.7923E-02 knorm= 6.7021E-01 apxde= 1.0406E-02    *not conv.*     

               starting mcscf iteration...   2

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    10, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 131071506

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130905714
 address segment size,           sizesg = 130754228
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:         80 transformed 1/r12    array elements were written in       1 records.


 mosort: allocated sort2 space, avc2is=   130933359 available sort2 space, avcisx=   130933611

   2 trial vectors read from nvfile (unit= 29).
 ciiter=   1 noldhv=  2 noldv=  2

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*       -1.1493815347       -1.8636672489        0.0000000000        0.0000010000
    2         0.5319949822       -0.1822907321        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.492705475920574E-003
 Total number of micro iterations:    4

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99955476 pnorm= 0.0000E+00 rznorm= 1.2025E-16 rpnorm= 0.0000E+00 noldr=  4 nnewr=  4 nolds=  0 nnews=  0
 

 qvv(*) eigenvalues. symmetry block  1
     1.545315    1.808486    3.922248    3.922248    5.413293    5.863606    5.863606    8.979150

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    2 emc=     -1.1493815347 demc= 1.2997E-02 wnorm= 2.7942E-02 knorm= 2.9838E-02 apxde= 1.4022E-04    *not conv.*     

               starting mcscf iteration...   3

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    10, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 131071506

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130905714
 address segment size,           sizesg = 130754228
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:         80 transformed 1/r12    array elements were written in       1 records.


 mosort: allocated sort2 space, avc2is=   130933359 available sort2 space, avcisx=   130933611

   2 trial vectors read from nvfile (unit= 29).
 ciiter=   1 noldhv=  2 noldv=  2

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*       -1.1495337540       -1.8638194683        0.0000000000        0.0000010000
    2         0.4812048319       -0.2330808823        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.019350532485659E-004
 Total number of micro iterations:    3

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99997695 pnorm= 0.0000E+00 rznorm= 8.9068E-07 rpnorm= 0.0000E+00 noldr=  3 nnewr=  3 nolds=  0 nnews=  0
 

 qvv(*) eigenvalues. symmetry block  1
     1.546809    1.836975    3.926486    3.926486    5.418039    5.867837    5.867837    9.022068

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    3 emc=     -1.1495337540 demc= 1.5222E-04 wnorm= 1.6155E-03 knorm= 6.7900E-03 apxde= 2.0936E-06    *not conv.*     

               starting mcscf iteration...   4

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    10, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 131071506

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130905714
 address segment size,           sizesg = 130754228
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:         80 transformed 1/r12    array elements were written in       1 records.


 mosort: allocated sort2 space, avc2is=   130933359 available sort2 space, avcisx=   130933611

   2 trial vectors read from nvfile (unit= 29).
 ciiter=   1 noldhv=  2 noldv=  2

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*       -1.1495361920       -1.8638219063        0.0000000000        0.0000010000
    2         0.4681948155       -0.2460908988        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.568580343136459E-005
 Total number of micro iterations:    3

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99999939 pnorm= 0.0000E+00 rznorm= 2.0783E-07 rpnorm= 0.0000E+00 noldr=  3 nnewr=  3 nolds=  0 nnews=  0
 

 qvv(*) eigenvalues. symmetry block  1
     1.546779    1.850581    3.926447    3.926447    5.418022    5.867818    5.867818    9.024121

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    4 emc=     -1.1495361920 demc= 2.4379E-06 wnorm= 2.8549E-04 knorm= 1.1035E-03 apxde= 5.4915E-08    *not conv.*     

               starting mcscf iteration...   5

 orbital-state coupling will be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    10, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 131071506

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130905714
 address segment size,           sizesg = 130754228
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:         80 transformed 1/r12    array elements were written in       1 records.


 mosort: allocated sort2 space, avc2is=   130933359 available sort2 space, avcisx=   130933611

   2 trial vectors read from nvfile (unit= 29).
 ciiter=   1 noldhv=  2 noldv=  2

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*       -1.1495362556       -1.8638219699        0.0000000000        0.0000010000
    2         0.4661133808       -0.2481723335        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  5.749434997093869E-006
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 Total number of micro iterations:    3

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99999998 pnorm= 1.5305E-05 rznorm= 3.8515E-08 rpnorm= 6.6144E-09 noldr=  3 nnewr=  3 nolds=  1 nnews=  1
 

 qvv(*) eigenvalues. symmetry block  1
     1.546775    1.852811    3.926443    3.926443    5.418021    5.867816    5.867816    9.024394

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    5 emc=     -1.1495362556 demc= 6.3681E-08 wnorm= 4.5995E-05 knorm= 2.0849E-04 apxde= 1.6536E-09    *not conv.*     

               starting mcscf iteration...   6

 orbital-state coupling will be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    10, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 131071506

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130905714
 address segment size,           sizesg = 130754228
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:         80 transformed 1/r12    array elements were written in       1 records.


 mosort: allocated sort2 space, avc2is=   130933359 available sort2 space, avcisx=   130933611

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   1 noldhv=  2 noldv=  2

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*       -1.1495362573       -1.8638219716        0.0000000072        0.0000010000
    2        -0.1976865417       -0.9119722559        0.0000000159        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  4.943990634410467E-009
 performing all-state projection
 Total number of micro iterations:    1

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 1.00000000 pnorm= 0.0000E+00 rznorm= 2.5414E-08 rpnorm= 2.6967E-09 noldr=  1 nnewr=  1 nolds=  0 nnews=  0
 

 qvv(*) eigenvalues. symmetry block  1
     1.546774    1.853233    3.926442    3.926442    5.418021    5.867816    5.867816    9.024446

 restrt: restart information saved on the restart file (unit= 13).

 all mcscf convergence criteria are satisfied.

 final mcscf convergence values:
 iter=    6 emc=     -1.1495362573 demc= 1.6534E-09 wnorm= 3.9552E-08 knorm= 5.8775E-09 apxde= 1.7998E-16    *converged*     




   ---------Individual total energies for all states:----------
   DRT #1 state # 1 wt 1.000 total energy=       -1.149536257, rel. (eV)=   0.000000
   ------------------------------------------------------------



          mcscf orbitals of the final iteration,    a block   1

               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
   1H__1s     0.32520646     0.73569641    -0.74410767     0.64946994    -0.00000000     0.00000000     0.38413453     0.00000000
   2H__1s     0.26987860     0.55241053     0.67935327    -2.14306234    -0.00000000    -0.00000000    -0.18757311    -0.00000000
   3H__2p    -0.01531090     0.04340487     0.03532403     0.20911429    -0.00000000    -0.00000000     0.64759192     0.00000000
   4H__2p    -0.00000000    -0.00000000    -0.00000000     0.00000000    -0.59854710    -0.12165123     0.00000000    -0.75260972
   5H__2p     0.00000000     0.00000000    -0.00000000    -0.00000000     0.12165123    -0.59854710    -0.00000000     0.43757356
   6H__1s     0.32520646    -0.73569641    -0.74410767    -0.64946994    -0.00000000     0.00000000     0.38413453     0.00000000
   7H__1s     0.26987860    -0.55241053     0.67935327     2.14306234     0.00000000    -0.00000000    -0.18757311    -0.00000000
   8H__2p     0.01531090     0.04340487    -0.03532403     0.20911429     0.00000000     0.00000000    -0.64759192    -0.00000000
   9H__2p     0.00000000     0.00000000    -0.00000000     0.00000000    -0.59854710    -0.12165123    -0.00000000     0.75260972
  10H__2p    -0.00000000     0.00000000    -0.00000000     0.00000000     0.12165123    -0.59854710     0.00000000    -0.43757356

               MO    9        MO   10
   1H__1s    -0.00000000    -1.38526395
   2H__1s    -0.00000000    -0.14358453
   3H__2p    -0.00000000     1.58716672
   4H__2p     0.43757356     0.00000000
   5H__2p     0.75260972    -0.00000000
   6H__1s    -0.00000000     1.38526395
   7H__1s     0.00000000     0.14358453
   8H__2p     0.00000000     1.58716672
   9H__2p    -0.43757356    -0.00000000
  10H__2p    -0.75260972     0.00000000

          natural orbitals of the final iteration,block  1    -    a

               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8

  occ(*)=     1.97710174     0.02289826     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
 
   1H__1s     0.32520646     0.73569641    -0.74410767     0.64943811    -0.00000000     0.00000000     0.38413451     0.00000000
   2H__1s     0.26987860     0.55241053     0.67935325    -2.14306565    -0.00000000     0.00000000    -0.18757309    -0.00000000
   3H__2p    -0.01531090     0.04340487     0.03532401     0.20915077     0.00000000    -0.00000000     0.64759191     0.00000000
   4H__2p    -0.00000000    -0.00000000    -0.00000000     0.00000000    -0.49889015    -0.35237235     0.00000000    -0.77171617
   5H__2p     0.00000000     0.00000000    -0.00000000     0.00000000     0.35237235    -0.49889015    -0.00000000     0.40292203
   6H__1s     0.32520646    -0.73569641    -0.74410769    -0.64943809     0.00000000    -0.00000000     0.38413450     0.00000000
   7H__1s     0.26987860    -0.55241053     0.67935329     2.14306563     0.00000000    -0.00000000    -0.18757309    -0.00000000
   8H__2p     0.01531090     0.04340487    -0.03532401     0.20915077    -0.00000000     0.00000000    -0.64759192    -0.00000000
   9H__2p     0.00000000     0.00000000     0.00000000    -0.00000000    -0.49889016    -0.35237235    -0.00000000     0.77171617
  10H__2p    -0.00000000     0.00000000    -0.00000000    -0.00000000     0.35237235    -0.49889016     0.00000000    -0.40292203

               MO    9        MO   10

  occ(*)=     0.00000000     0.00000000
 
   1H__1s    -0.00000000    -1.38527888
   2H__1s     0.00000000    -0.14353528
   3H__2p     0.00000000     1.58716192
   4H__2p     0.40292203    -0.00000000
   5H__2p     0.77171617     0.00000000
   6H__1s     0.00000000     1.38527888
   7H__1s    -0.00000000     0.14353528
   8H__2p    -0.00000000     1.58716192
   9H__2p    -0.40292203     0.00000000
  10H__2p    -0.77171617    -0.00000000
 d1(*), fmc(*), and qmc(*) written to the 1-particle density matrix file.
          5 d2(*) elements written to the 2-particle density matrix file: mcd2fl                                                      


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                          a partial gross atomic populations
   ao class       1  a       2  a       3  a       4  a       5  a       6  a
      _ s       0.982008   0.010994   0.000000   0.000000   0.000000   0.000000
      _ p       0.006543   0.000455   0.000000   0.000000   0.000000   0.000000
     1_ s       0.982008   0.010994   0.000000   0.000000   0.000000   0.000000
     1_ p       0.006543   0.000455   0.000000   0.000000   0.000000   0.000000
 
   ao class       7  a       8  a       9  a      10  a


                        gross atomic populations
     ao             _         1_
      s         0.993002   0.993002
      p         0.006998   0.006998
    total       1.000000   1.000000
 

 Total number of electrons:    2.00000000

 !timer: mcscf                           cpu_time=     0.020 walltime=     0.023
