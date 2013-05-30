module statistics
  use cdata
  contains
  !*********************************************************************
  !**************ROUTINES TO GENERATE RANDOM VARIABLES******************
  !*********************************************************************
  !used to draw random variables from a \f$N(0,1)\f$ distribution, 
  !uses box-muller algorithm
  real function rnorm(mu,sigma2)
    implicit none
    real, intent(IN) :: mu, sigma2
    real :: u1, u2
    call random_number(u1) ;  call random_number(u2)
    rnorm=mu+sqrt(sigma2)*sqrt(-2.*log(u1))*cos(2.*pi*u2)
  end function
  !*********************************************************************
  !used to draw random variables from laplace distribution
  !!http://en.wikipedia.org/wiki/Laplace_distributio
  real function rlaplace(mu,b) !NOT TESTED
    implicit none
    real, intent(IN) :: mu, b
    real :: u
    call random_number(u)
    u=u-0.5
    rlaplace=mu-b*sign(1.,u)*log(1.-2.*abs(u))
  end function
  !*********************************************************************
  !draw random variables from a \f$U(\alpha,\beta)\f$ distribution
  real function runif(alpha,beta)
    implicit none
    real, intent(IN) :: alpha, beta
    real :: u
    call random_number(u)
    runif=alpha+u*(beta-alpha)
  end function
end module