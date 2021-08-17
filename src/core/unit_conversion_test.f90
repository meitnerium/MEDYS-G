Program unit_conversion_test

use basics
use unit_conversion
Implicit None
type(frequency)::freq1

logical:: test_conversion_factors=.false.!.False.
logical:: test_exponent_mantissa=.false.
logical:: test_dimensional_value=.True.

if(test_dimensional_value)then
write(*,*) "TEST DIMENSIONAL DERIVED TYPES"
freq1=frequency(1., "Hz")
write(*,*) freq1%value, " ", freq1%unit_
freq1=frequency_to_au(freq1)
write(*,*) freq1%value, " ", freq1%unit_
end if


end program
