module deriv
  use cdata
  implicit none
  contains
  real function first_deriv(vect,h) !real function as used for gradient of phase
    !TAKE THE 1ST DERIVATIVE OF THE VECT ARRAY
    implicit none
    real, intent(IN) :: vect(-3:3)
    real,intent(IN) :: h
    first_deriv=(-vect(-3)+9*vect(-2)-45*vect(-1)+45*vect(1)-9*vect(2)+vect(3))/(60*h)
  end function
  complex function second_deriv(vect,h) !complex function 
    !TAKE THE 2ND DERIVATIVE OF THE VECT ARRAY
    implicit none
    complex, intent(IN) :: vect(-3:3)
    real,intent(IN) :: h
    second_deriv=(2*vect(-3)-27*vect(-2)+270*vect(-1)-490*vect(0)+270*vect(1)- 27*vect(2)+2*vect(3))/(180*(h**2))
  end function
  real function real_second_deriv(vect,h) !real function
    !TAKE THE 2ND DERIVATIVE OF THE VECT ARRAY - ONLY 4th ORDER
    implicit none
    real, intent(IN) :: vect(-2:2)
    real,intent(IN) :: h
    real_second_deriv=(-vect(-2)+16*vect(-1)-30*vect(0)+16*vect(1)-vect(2))/(12*(h**2))
  end function
  real function low_first_deriv(vect,h) !real function
    !TAKE THE 1ST DERIVATIVE OF THE VECT ARRAY - ONLY 2ND ORDER
    implicit none
    real, intent(IN) :: vect(-1:1)
    real,intent(IN) :: h
    low_first_deriv=(-vect(-1)+vect(1))/(2*h)
  end function
end module
