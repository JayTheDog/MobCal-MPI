module constants
   use xtb_mctc_accuracy, only : wp
   implicit none
   private

  ! from Handbook of Chemistry and Physics, 70th Edition
   real(wp),public,parameter :: pi = 3.1415926535897932384626433832795029_wp
  real(wp), public, parameter :: cang=180.d0/pi
  real(wp), public, parameter :: xe=1.60217733d-19
  real(wp), public, parameter :: xk=1.380658d-23
  real(wp), public, parameter :: xn=6.0221367d23
  real(wp), public, parameter :: xeo=8.854187817d-12

end module constants
