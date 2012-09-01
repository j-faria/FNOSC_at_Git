! *****************************************************
! *****************************************************
	subroutine read_sun (p)
!  this routine reads in oscillation frequencies from
!  observations of the sun (for l=0 modes) and stores
!  the periods in the array P

	  implicit none
	  
	  character(len=1)       :: iyorn
	  integer                :: i
	  integer, dimension(15) :: lobs
	  integer, dimension(15) :: nobs	
	  real, dimension(15)    :: nu  
	  real, dimension(15)    :: sigma
	  real, intent(inout)    :: p(*)

       character(len=20)      :: infile
!
	write (6,*) ' your frequency file is'
	read (5,3000) infile
 3000 format (a20)
	open (unit=12, file=infile, status='old', action='read')


	! read in frequencies
	do i=1,15
     	read (12,*) lobs(i),nobs(i),nu(i),sigma(i)
     	p(i) = 1.0d6 / nu(i)
	enddo

	return
	
	end subroutine read_sun
