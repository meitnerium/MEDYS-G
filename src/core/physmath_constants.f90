Module physmath_constants_module

use precision_module
Implicit none
Public

Save

Real(kind=real_8),parameter :: pi=3.141592653589793_real_8

Real(kind=real_8),parameter :: twopi=2.0_real_8*pi
Real(kind=real_8),parameter :: pio2=pi/2.0_real_8

Complex(kind=comp_16),parameter:: one_imag=(0._real_8,1._real_8)
End module
