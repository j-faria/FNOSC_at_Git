! ****************************************************
! ****************************************************
	subroutine rkfrad(rt,y,yp)
!
! this routine supplies the derivatives for use in rkf
! for the radial case.
!

	  use consts
	  
	  implicit double precision (a-h,o-z)
	  common/misc/l,lhat,lindex,nsurf,period, &
          eps,verg,eig,eigt,nodes1,nodes2, &
          modep
	  double precision l,lhat,lindex
	  dimension y(2),yp(2)
	  common/splinq/vq(6)
!
	call splntl(rt,vq)

	yp(1)=-(3.0d0*y(1)+y(2)/vq(3))				! eq 8.11 ?
	yp(2)=vq(4)*(4.0d0*y(1)+eigt*y(1)/vq(1)+y(2))	! eq 8.12 ? 

	do i=1,2
		yp(i)=yp(i)*vq(5)
	enddo
	
	return
	
	end subroutine rkfrad
