! ****************************************************
! ****************************************************
	subroutine output(outfile, l, nodes1, modep, period)
!
! this routine outputs the mode frequencies to the
! file OUTFILE (unit 13)
!
	  use commonvar, only: radial
	  
	  implicit none
	  
	character(len=20), intent(in)  :: outfile	! output file
	integer, intent(in)            :: l 		! angular degree
	integer, intent(in)            :: nodes1	! radial degree (if l=0)
	integer, intent(in)            :: modep		! radial degree (if l>0)
	double precision, intent(in)   :: period	! mode period


	if( radial ) then
		write(13,'(i2,i4,f14.6)') l, nodes1+1, 1.d6/period
	else
		write(13,'(i2,i4,f14.6)') l, -modep, 1.d6/period
	end if
		 
			
	end subroutine output
