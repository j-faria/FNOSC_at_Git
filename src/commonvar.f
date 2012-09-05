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
     !! the following have dimensions depending on the
     !! number of points in the equilibrium model
     real(dp), dimension(:), allocatable    :: r
     real(dp), dimension(:), allocatable    :: g
     real(dp), dimension(:), allocatable    :: rho
     real(dp), dimension(1800)               :: x
!     real(dp), dimension(:), allocatable    :: x
!     real(dp), dimension(2,200)             :: yeig
     real(dp), dimension(:,:), allocatable  :: yeig
     
    
!     real(dp), dimension(300)    :: xi
!     real(dp), dimension(300)    :: y1
!     real(dp), dimension(300)    :: y2
!     real(dp), dimension(300)    :: y3
!     real(dp), dimension(300)    :: y4
!     real(dp), dimension(300)    :: y5 
     
     real(dp), dimension(:), allocatable    :: xi
     real(dp), dimension(:), allocatable    :: y1
     real(dp), dimension(:), allocatable    :: y2
     real(dp), dimension(:), allocatable    :: y3
     real(dp), dimension(:), allocatable    :: y4
     real(dp), dimension(:), allocatable    :: y5
     
     
!     real(dp), dimension(4,300)    :: c1
!     real(dp), dimension(4,300)    :: c2
!     real(dp), dimension(4,300)    :: c3
!     real(dp), dimension(4,300)    :: c4
!     real(dp), dimension(4,300)    :: c5
!     real(dp), dimension(4,300)    :: c6

     real(dp), dimension(:,:), allocatable    :: c1
     real(dp), dimension(:,:), allocatable    :: c2
     real(dp), dimension(:,:), allocatable    :: c3
     real(dp), dimension(:,:), allocatable    :: c4
     real(dp), dimension(:,:), allocatable    :: c5
     real(dp), dimension(:,:), allocatable    :: c6
     
!	real(dp), dimension(200)               :: r
!	real(dp), dimension(200)               :: g
!	real(dp), dimension(200)               :: rho


end module commonvar
