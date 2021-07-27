                                program "scfpq"
                            columbus program system
             restricted hartree-fock scf and two-configuration mcscf

                   programmed (in part) by russell m. pitzer
                            version date: 30-sep-00

 This Version of Program SCFPQ is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de



     ******************************************
     **    PROGRAM:              SCFPQ       **
     **    PROGRAM VERSION:      5.5         **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************

 echo of the input file:
 ------------------------------------------------------------------------
 Default SCF Title
  -15   9   0  -1   0   8   0   0   0  40  20   1   1   0   0   0  -1   0  11   0
    1
    2
    0.90000   0.10000   0.05000   0.00000    0    0
  (4f18.12)
    1
 ------------------------------------------------------------------------

workspace parameters: lcore= 131072000
 mem1=          0 ifirst=          1

Default SCF Title                                                               

scf flags
 -15   9   0  -1   0   8   0   0   0  40  20   1   1   0   0   0  -1   0  11   0
   1
input orbital labels, i:bfnlab(i)=
   1:  1H__1s   2:  2H__1s   3:  3H__2p   4:  4H__2p   5:  5H__2p   6:  6H__1s
   7:  7H__1s   8:  8H__2p   9:  9H__2p  10: 10H__2p
bfn_to_center map(*), i:map(i)
   1:  1   2:  1   3:  1   4:  1   5:  1   6:  2   7:  2   8:  2   9:  2  10:  2
bfn_to_orbital_type map(*), i:map(i)
   1:  1   2:  1   3:  2   4:  2   5:  2   6:  1   7:  1   8:  2   9:  2  10:  2

                                                                                

core potentials used

normalization threshold = 10**(-20)
one- and two-electron energy convergence criterion = 10**(- 9)

nuclear repulsion energy =     0.7142857142857
 
 fock matrix built from AO INTEGRALS
                        ^^^^^^^^^^^^
  ao integrals in SIFS format 
 in-core processing switched off

DIIS SWITCHED ON (error vector is FDS-SDF)

total number of SCF basis functions:  10

                                                                                
Default SCF Title                                                               

                 55 s integrals
                 55 t integrals
                 55 v integrals
                 55 c integrals
 check on positive semi-definiteness of S passed .
                 55 h integrals

           starting vectors from diagonalization of one-electron terms

     occupied and virtual orbitals
                                                                                
Default SCF Title                                                               

                           a molecular orbitals
      sym. orb.       1  a        2  a        3  a        4  a        5  a
   1  a,  1H__1s    0.425069   -0.369435   -0.662168    0.889288    0.000000
   2  a,  2H__1s    0.167559   -1.293279    0.701193   -1.786631    0.000000
   3  a,  3H__2p   -0.037497   -0.048948    0.068984    0.193981    0.000000
   4  a,  4H__2p    0.000000    0.000000    0.000000    0.000000    0.610784
   5  a,  5H__2p    0.000000    0.000000    0.000000    0.000000    0.000000
   6  a,  6H__1s    0.425069    0.369435   -0.662168   -0.889288    0.000000
   7  a,  7H__1s    0.167560    1.293279    0.701193    1.786631    0.000000
   8  a,  8H__2p    0.037497   -0.048948   -0.068984    0.193981    0.000000
   9  a,  9H__2p    0.000000    0.000000    0.000000    0.000000    0.610784
  10  a, 10H__2p    0.000000    0.000000    0.000000    0.000000    0.000000

      orb. en.     -1.279908   -0.605169   -0.171773    0.223678    0.677101

      occ. no.      2.000000    0.000000    0.000000    0.000000    0.000000

      sym. orb.       6  a        7  a        8  a        9  a       10  a
   1  a,  1H__1s    0.000000    0.433433    0.000000    0.000000   -1.398106
   2  a,  2H__1s    0.000000   -0.223138    0.000000    0.000000   -0.232101
   3  a,  3H__2p    0.000000    0.643966    0.000000    0.000000    1.588926
   4  a,  4H__2p    0.000000    0.000000   -0.870570    0.000000    0.000000
   5  a,  5H__2p    0.610784    0.000000    0.000000   -0.870570    0.000000
   6  a,  6H__1s    0.000000    0.433433    0.000000    0.000000    1.398106
   7  a,  7H__1s    0.000000   -0.223138    0.000000    0.000000    0.232101
   8  a,  8H__2p    0.000000   -0.643966    0.000000    0.000000    1.588926
   9  a,  9H__2p    0.000000    0.000000    0.870570    0.000000    0.000000
  10  a, 10H__2p    0.610784    0.000000    0.000000    0.870570    0.000000

      orb. en.      0.677101    1.412124    1.662847    1.662847    3.195296

      occ. no.      0.000000    0.000000    0.000000    0.000000    0.000000


mo coefficients will be saved after each iteration

iteration       energy           one-electron energy       two-electron energy
         ttr          tmr           tnr

 total 2e-integrals processed:                   512
    1     -1.0691853341566         -2.5598157903865         0.77634474194424    
       0.9000       0.1000       0.5000E-01
diis info stored; extrapolation off
 total 2e-integrals processed:                   512
    2     -1.1153850002465         -2.5454007006879         0.71572998615571    
       0.9000       0.1000       0.5000E-01
diis info stored; extrapolation off
 total 2e-integrals processed:                   512
    3     -1.1269285268867         -2.5281708598803         0.68695661870782    
       0.9000       0.1000       0.5000E-01
 total 2e-integrals processed:                   512
    4     -1.1271190439708         -2.5276800112223         0.68627525296579    
       0.8000       0.1000       0.5000E-01
 total 2e-integrals processed:                   512
    5     -1.1290347036456         -2.5215418099965         0.67822139206523    
       0.7000       0.1000       0.5000E-01
 total 2e-integrals processed:                   512
    6     -1.1309876623146         -2.5094398223317         0.66416644573138    
       0.6000       0.1000       0.5000E-01
 total 2e-integrals processed:                   512
    7     -1.1312522295089         -2.5043846319540         0.65884668815942    
       0.5000       0.1000       0.5000E-01
 total 2e-integrals processed:                   512
    8     -1.1312817519834         -2.5025552134805         0.65698774721145    
       0.4000       0.1000       0.5000E-01
 total 2e-integrals processed:                   512
    9     -1.1312842119007         -2.5019899761919         0.65642005000547    
       0.3000       0.1000       0.5000E-01
 total 2e-integrals processed:                   512
   10     -1.1312843460273         -2.5018484959092         0.65627843559621    
       0.2000       0.1000       0.5000E-01
 total 2e-integrals processed:                   512
   11     -1.1312843498241         -2.5018227687189         0.65625270460905    
       0.1000       0.1000       0.5000E-01
 total 2e-integrals processed:                   512
   12     -1.1312843498555         -2.5018204295931         0.65625036545185    
       0.1000       0.1000       0.5000E-01
 total 2e-integrals processed:                   512
   13     -1.1312843498558         -2.5018202168876         0.65625015274605    
       0.1000       0.1000       0.5000E-01
 total 2e-integrals processed:                   512
   14     -1.1312843498558         -2.5018201975421         0.65625013340061    
       0.1000       0.1000       0.5000E-01

calculation has *converged*
scf gradient information written to file vectgrd

     total energy =         -1.1312843499
     kinetic energy =        1.1275379882
     potential energy =     -2.2588223380
     virial theorem =        1.9966883996
     wavefunction norm =     1.0000000000

     occupied and virtual orbitals
                                                                                
Default SCF Title                                                               

                           a molecular orbitals
      sym. orb.       1  a        2  a        3  a        4  a        5  a
   1  a,  1H__1s    0.318259    0.120418   -0.746082    0.882010    0.000000
   2  a,  2H__1s    0.275398    1.718477    0.677406   -1.392469    0.000000
   3  a,  3H__2p   -0.018262   -0.002371    0.035784    0.281534    0.000000
   4  a,  4H__2p    0.000000    0.000000    0.000000    0.000000    0.000000
   5  a,  5H__2p    0.000000    0.000000    0.000000    0.000000    0.610784
   6  a,  6H__1s    0.318259   -0.120418   -0.746082   -0.882010    0.000000
   7  a,  7H__1s    0.275398   -1.718477    0.677406    1.392469    0.000000
   8  a,  8H__2p    0.018262   -0.002371   -0.035784    0.281534    0.000000
   9  a,  9H__2p    0.000000    0.000000    0.000000    0.000000    0.000000
  10  a, 10H__2p    0.000000    0.000000    0.000000    0.000000    0.610784

      orb. en.     -0.594660    0.239318    0.771532    1.309808    1.959217

      occ. no.      2.000000    0.000000    0.000000    0.000000    0.000000

      sym. orb.       6  a        7  a        8  a        9  a       10  a
   1  a,  1H__1s    0.000000    0.386119    0.000000    0.000000   -1.445536
   2  a,  2H__1s    0.000000   -0.186592    0.000000    0.000000   -0.162355
   3  a,  3H__2p    0.000000    0.647490    0.000000    0.000000    1.576529
   4  a,  4H__2p    0.610784    0.000000   -0.870570    0.000000    0.000000
   5  a,  5H__2p    0.000000    0.000000    0.000000   -0.870570    0.000000
   6  a,  6H__1s    0.000000    0.386119    0.000000    0.000000    1.445536
   7  a,  7H__1s    0.000000   -0.186592    0.000000    0.000000    0.162355
   8  a,  8H__2p    0.000000   -0.647490    0.000000    0.000000    1.576529
   9  a,  9H__2p    0.610784    0.000000    0.870570    0.000000    0.000000
  10  a, 10H__2p    0.000000    0.000000    0.000000    0.870570    0.000000

      orb. en.      1.959217    2.704590    2.930152    2.930152    4.530253

      occ. no.      0.000000    0.000000    0.000000    0.000000    0.000000


     population analysis
                                                                                
Default SCF Title                                                               
  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                          a partial gross atomic populations
   ao class       1  a
      5 s       0.992043
      5 p       0.007957
     10 s       0.992043
     10 p       0.007957


                        gross atomic populations
     ao             5         10
      s         0.992043   0.992043
      p         0.007957   0.007957
    total       1.000000   1.000000

 Total number of electrons:    2.00000000

