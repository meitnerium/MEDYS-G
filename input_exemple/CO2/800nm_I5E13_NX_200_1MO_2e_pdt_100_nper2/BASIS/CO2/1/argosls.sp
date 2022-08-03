echo of the argos input file:
 ------------------------------------------------------------------------
 CO2_6-31G
   0  3  2 15  2 20  9  0  0  0  0  0  0  0  0
   1  1  a
   0
     1    1
     3    1    1    1
     1    1    1
     1
     3    3    2
     1    0    0
     0    1    0
     0    0    1
     6    1    1
    3047.5249000        0.0018347
     457.3695100        0.0140373
     103.9486900        0.0688426
      29.2101550        0.2321844
       9.2866630        0.4679413
       3.1639270        0.3623120
     3    1    1
       7.8682724       -0.1193324
       1.8812885       -0.1608542
       0.5442493        1.1434564
     1    1    1
       0.1687144        1.0000000
     3    2    1
       7.8682724        0.0689991
       1.8812885        0.3164240
       0.5442493        0.7443083
     1    2    1
       0.1687144        1.0000000
     6    1    1
    5484.6717000        0.0018311
     825.2349500        0.0139501
     188.0469600        0.0684451
      52.9645000        0.2327143
      16.8975700        0.4701930
       5.7996353        0.3585209
     3    1    1
      15.5396160       -0.1107775
       3.5999336       -0.1480263
       1.0137618        1.1307670
     1    1    1
       0.2700058        1.0000000
     3    2    1
      15.5396160        0.0708743
       3.5999336        0.3397528
       1.0137618        0.7271586
     1    2    1
       0.2700058        1.0000000
     6    1    1
    5484.6717000        0.0018311
     825.2349500        0.0139501
     188.0469600        0.0684451
      52.9645000        0.2327143
      16.8975700        0.4701930
       5.7996353        0.3585209
     3    1    1
      15.5396160       -0.1107775
       3.5999336       -0.1480263
       1.0137618        1.1307670
     1    1    1
       0.2700058        1.0000000
     3    2    1
      15.5396160        0.0708743
       3.5999336        0.3397528
       1.0137618        0.7271586
     1    2    1
       0.2700058        1.0000000
 C    5  1 6.
     0.00000000    0.00000000    0.00000000
   1  1
   2  1
   3  1
   4  2
   5  2
 O    5  1 8.
     0.00000000    0.00000000   -2.19208232
   6  1
   7  1
   8  1
   9  2
  10  2
 O    5  1 8.
     0.00000000    0.00000000    2.19208232
  11  1
  12  1
  13  1
  14  2
  15  2
 ------------------------------------------------------------------------
                              program "argos" 5.9
                            columbus program system

             this program computes integrals over symmetry orbitals
               of generally contracted gaussian atomic orbitals.
                        programmed by russell m. pitzer

                           version date: 20-aug-2001

references:

    symmetry analysis (equal contributions):
    r. m. pitzer, j. chem. phys. 58, 3111 (1973).

    ao integral evaluation (hondo):
    m. dupuis, j. rys, and h. f. king, j. chem. phys. 65, 111 (1976).

    general contraction of gaussian orbitals:
    r. c. raffenetti, j. chem. phys. 58, 4452 (1973).

    core potential ao integrals (meldps):
    l. e. mcmurchie and e. r. davidson, j. comput. phys. 44, 289 (1981)

    spin-orbit and core potential integrals:
    r. m. pitzer and n. w. winter, int. j. quantum chem. 40, 773 (1991)

 This Version of Program ARGOS is Maintained by:
    Hans Lischka
    Institute for Theoretical Chemistry
    University of Vienna
    Waeringerstr 17, A-1090 Wien, Austria
    Internet: hans.lischka@univie.ac.at

    Thomas Mueller
    Central Institute of Applied Mathematics
    Research Centre Juelich
    D-52425 Juelich, Germany
    Internet: th.mueller@fz-juelich.de


 workspace allocation parameters: lcore=   3000000 mem1=         0 ifirst=         1

 filenames and unit numbers:
 unit description                   filename
 ---- -----------                   ----------
   6  listing file:                 argosls                                                     
   4  1-e integral file:            aoints                                                      
   8  2-e integral file [fsplit=2]: aoints2                                                     
   5  input file:                   argosin                                                     

 argos input parameters and titles:
 ngen   =  0 ns     =  3 naords =  2 ncons  = 15 ngcs   =  2 itol   = 20
 icut   =  9 aoints =  4 only1e =  0 inrm   =  0 ncrs   =  0
 l1rec  =         0      l2rec  =         0      aoint2 =  8 fsplit =  2

CO2_6-31G                                                                       
aoints SIFS file created by argos.      colosse1          12:04:00.930 23-Jun-14


irrep            1
degeneracy       1
label             a


direct product table
   (  a) (  a) =   a


                     nuclear repulsion energy   58.39196769


primitive ao integrals neglected if exponential factor below 10**(-20)
contracted ao and so integrals neglected if value below 10**(- 9)
symmetry orbital integrals written on units  4  8


                                    C   atoms

                              nuclear charge   6.00

           center            x               y               z
             1             0.00000000      0.00000000      0.00000000

                    1s orbitals

 orbital exponents  contraction coefficients
    3047.525       1.8347002E-03
    457.3695       1.4037301E-02
    103.9487       6.8842606E-02
    29.21016       0.2321844    
    9.286663       0.4679413    
    3.163927       0.3623120    

                     symmetry orbital labels
                     1  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    1s orbitals

 orbital exponents  contraction coefficients
    7.868272      -0.1193324    
    1.881288      -0.1608542    
   0.5442493        1.143456    

                     symmetry orbital labels
                     2  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    1s orbitals

 orbital exponents  contraction coefficients
   0.1687144        1.000000    

                     symmetry orbital labels
                     3  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    2p orbitals

 orbital exponents  contraction coefficients
    7.868272       6.8999096E-02
    1.881288       0.3164240    
   0.5442493       0.7443083    

                     symmetry orbital labels
                     4  a1
                     5  a1
                     6  a1

           symmetry orbitals
 ctr, ao     a1     a1     a1
  1, 100  1.000  0.000  0.000
  1, 010  0.000  1.000  0.000
  1, 001  0.000  0.000  1.000

                    2p orbitals

 orbital exponents  contraction coefficients
   0.1687144        1.000000    

                     symmetry orbital labels
                     7  a1
                     8  a1
                     9  a1

           symmetry orbitals
 ctr, ao     a1     a1     a1
  1, 100  1.000  0.000  0.000
  1, 010  0.000  1.000  0.000
  1, 001  0.000  0.000  1.000


                                    O   atoms

                              nuclear charge   8.00

           center            x               y               z
             1             0.00000000      0.00000000     -2.19208232

                    1s orbitals

 orbital exponents  contraction coefficients
    5484.672       1.8310998E-03
    825.2350       1.3950099E-02
    188.0470       6.8445093E-02
    52.96450       0.2327143    
    16.89757       0.4701930    
    5.799635       0.3585209    

                     symmetry orbital labels
                    10  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    1s orbitals

 orbital exponents  contraction coefficients
    15.53962      -0.1107775    
    3.599934      -0.1480263    
    1.013762        1.130767    

                     symmetry orbital labels
                    11  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    1s orbitals

 orbital exponents  contraction coefficients
   0.2700058        1.000000    

                     symmetry orbital labels
                    12  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    2p orbitals

 orbital exponents  contraction coefficients
    15.53962       7.0874300E-02
    3.599934       0.3397528    
    1.013762       0.7271586    

                     symmetry orbital labels
                    13  a1
                    14  a1
                    15  a1

           symmetry orbitals
 ctr, ao     a1     a1     a1
  1, 100  1.000  0.000  0.000
  1, 010  0.000  1.000  0.000
  1, 001  0.000  0.000  1.000

                    2p orbitals

 orbital exponents  contraction coefficients
   0.2700058        1.000000    

                     symmetry orbital labels
                    16  a1
                    17  a1
                    18  a1

           symmetry orbitals
 ctr, ao     a1     a1     a1
  1, 100  1.000  0.000  0.000
  1, 010  0.000  1.000  0.000
  1, 001  0.000  0.000  1.000


                                    O   atoms

                              nuclear charge   8.00

           center            x               y               z
             1             0.00000000      0.00000000      2.19208232

                    1s orbitals

 orbital exponents  contraction coefficients
    5484.672       1.8310998E-03
    825.2350       1.3950099E-02
    188.0470       6.8445093E-02
    52.96450       0.2327143    
    16.89757       0.4701930    
    5.799635       0.3585209    

                     symmetry orbital labels
                    19  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    1s orbitals

 orbital exponents  contraction coefficients
    15.53962      -0.1107775    
    3.599934      -0.1480263    
    1.013762        1.130767    

                     symmetry orbital labels
                    20  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    1s orbitals

 orbital exponents  contraction coefficients
   0.2700058        1.000000    

                     symmetry orbital labels
                    21  a1

           symmetry orbitals
 ctr, ao     a1
  1, 000  1.000

                    2p orbitals

 orbital exponents  contraction coefficients
    15.53962       7.0874300E-02
    3.599934       0.3397528    
    1.013762       0.7271586    

                     symmetry orbital labels
                    22  a1
                    23  a1
                    24  a1

           symmetry orbitals
 ctr, ao     a1     a1     a1
  1, 100  1.000  0.000  0.000
  1, 010  0.000  1.000  0.000
  1, 001  0.000  0.000  1.000

                    2p orbitals

 orbital exponents  contraction coefficients
   0.2700058        1.000000    

                     symmetry orbital labels
                    25  a1
                    26  a1
                    27  a1

           symmetry orbitals
 ctr, ao     a1     a1     a1
  1, 100  1.000  0.000  0.000
  1, 010  0.000  1.000  0.000
  1, 001  0.000  0.000  1.000

lx:   a          ly:   a          lz:   a

output SIFS file header information:
CO2_6-31G                                                                       
aoints SIFS file created by argos.      colosse1          12:04:00.930 23-Jun-14

output energy(*) values:
 energy( 1)=  5.839196768851E+01, ietype=   -1,    core energy of type: Nuc.Rep.

total core energy =  5.839196768851E+01

nsym = 1 nbft=  27

symmetry  =    1
slabel(*) =    a
nbpsy(*)  =   27

info(*) =         2      4096      3272      4096      2700

output orbital labels, i:bfnlab(i)=
   1:  1C__1s   2:  2C__1s   3:  3C__1s   4:  4C__2p   5:  5C__2p   6:  6C__2p
   7:  7C__2p   8:  8C__2p   9:  9C__2p  10: 10O__1s  11: 11O__1s  12: 12O__1s
  13: 13O__2p  14: 14O__2p  15: 15O__2p  16: 16O__2p  17: 17O__2p  18: 18O__2p
  19: 19O__1s  20: 20O__1s  21: 21O__1s  22: 22O__2p  23: 23O__2p  24: 24O__2p
  25: 25O__2p  26: 26O__2p  27: 27O__2p

bfn_to_center map(*), i:map(i)
   1:  1   2:  1   3:  1   4:  1   5:  1   6:  1   7:  1   8:  1   9:  1  10:  2
  11:  2  12:  2  13:  2  14:  2  15:  2  16:  2  17:  2  18:  2  19:  3  20:  3
  21:  3  22:  3  23:  3  24:  3  25:  3  26:  3  27:  3

bfn_to_orbital_type map(*), i:map(i)
   1:  1   2:  1   3:  1   4:  2   5:  2   6:  2   7:  2   8:  2   9:  2  10:  1
  11:  1  12:  1  13:  2  14:  2  15:  2  16:  2  17:  2  18:  2  19:  1  20:  1
  21:  1  22:  2  23:  2  24:  2  25:  2  26:  2  27:  2


       27 symmetry orbitals,        a:  27
 !timer: syminp required                 cpu_time=     0.026 walltime=     0.362

 socfpd: mcxu=      174 mcxu2=      137 left=  2999826
 !timer: socfpd required                 cpu_time=     0.000 walltime=     0.012
 
oneint:   143 S1(*)    integrals were written in  1 records.
oneint:   143 T1(*)    integrals were written in  1 records.
oneint:   155 V1(*)    integrals were written in  1 records.
oneint:    78 X(*)     integrals were written in  1 records.
oneint:    78 Y(*)     integrals were written in  1 records.
oneint:   138 Z(*)     integrals were written in  1 records.
oneint:    78 Im(px)   integrals were written in  1 records.
oneint:    78 Im(py)   integrals were written in  1 records.
oneint:   116 Im(pz)   integrals were written in  1 records.
oneint:    72 Im(lx)   integrals were written in  1 records.
oneint:    72 Im(ly)   integrals were written in  1 records.
oneint:    36 Im(lz)   integrals were written in  1 records.
oneint:   143 XX(*)    integrals were written in  1 records.
oneint:    36 XY(*)    integrals were written in  1 records.
oneint:    82 XZ(*)    integrals were written in  1 records.
oneint:   143 YY(*)    integrals were written in  1 records.
oneint:    82 YZ(*)    integrals were written in  1 records.
oneint:   155 ZZ(*)    integrals were written in  1 records.
 
 !timer: oneint required                 cpu_time=     0.006 walltime=     0.087
 !timer: seg1mn required                 cpu_time=     0.032 walltime=     0.461

twoint:       20854 1/r12    integrals and     2959 pk flags
                                 were written in     8 records.

 twoint: maximum mblu needed =       682
 !timer: twoint required                 cpu_time=     0.105 walltime=     0.161
 
driver: 1-e  integral workspace high-water mark =     21169
driver: 2-e  integral workspace high-water mark =     14623
driver: overall argos workspace high-water mark =     21169
 !timer: argos required                  cpu_time=     0.138 walltime=     0.627
