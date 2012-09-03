! ****************************************************
! ****************************************************
	subroutine output(outfile, l, nodes1, modep, period)
!
! this routine outputs the mode frequencies to the
! file OUTFILE (unit 13)
!
		  implicit double precision (a-h,o-z)
		  logical radial
		  common/rad/radial


		if( radial ) then
			write(13,'(i2,i4,f14.6)') l, nodes1+1, 1.d6/period
		else
			write(13,'(i2,i4,f14.6)') l, -modep, 1.d6/period
		end if
		 
			
	end subroutine output
