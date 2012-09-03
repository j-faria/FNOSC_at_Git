!--------------------------------------------------------------------
!	Joao Faria: 08/2012 
!--------------------------------------------------------------------
!	 Module that contains the variables that some subroutines
!	 need to share. 

module commonvar
	
	implicit none

	! real number precision options: single, double, quad
     integer, parameter :: sp = selected_real_kind(p=5)
     integer, parameter :: dp = selected_real_kind(p=15)
     integer, parameter :: qp = selected_real_kind(p=30)
     
     
     ! flags
     integer   :: iprnt
     logical   :: radial
     
     ! arrays
     real(dp), dimension(200)              :: r
     real(dp), dimension(200)              :: g
     real(dp), dimension(200)              :: rho
     real(dp), dimension(200)              :: x
     real(dp), dimension(2,200)            :: yeig


end module commonvar
