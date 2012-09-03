!--------------------------------------------------------------------
!	Joao Faria: 08/2012 
!--------------------------------------------------------------------
!	 Module that contains the constants that some subroutines
!	 need to share. 

module consts
	
	implicit none
	
	! real number precision options: single, double, quad
     integer, parameter :: sp = selected_real_kind(p=5)
     integer, parameter :: dp = selected_real_kind(p=15)
     integer, parameter :: qp = selected_real_kind(p=30)
     
     
	! math constants 
	real(dp), parameter :: pi = 3.1415926535897932384626433832795029D0
	real(dp), parameter :: pi4 = 4*pi
	real(dp), parameter :: eulercon = 0.577215664901532861d0
	real(dp), parameter :: ln10 = 2.3025850929940459d0
	real(dp), parameter :: a2rad = pi/180.0d0 ! angle to radians
	real(dp), parameter :: rad2a = 180.0d0/pi ! radians to angle
	real(dp), parameter :: one_third = 1d0/3d0
	real(dp), parameter :: two_thirds = 2d0/3d0
	real(dp), parameter :: ln4pi3 = 1.43241195830118d0 ! log(4*pi/3)
	real(dp), parameter :: pi43 = 4.d0*pi/3.d0 ! 4*pi/3
	
	! physical constants   
	real(dp) :: grav = 6.67428d-8 
		! gravitational constant (g^-1 cm^3 s^-2)

end module consts
