! ***************************************************
! ***************************************************
	subroutine rkfcow(rt,y,yp)
!
! this routine supplies the derivatives for use in rkf
! for the nonradial case.
!
	  implicit double precision (a-h,o-z)
	  common/misc/l,lhat,lindex,nsurf,period,grav, &
          pi,pi4,p43,eps,verg,eig,eigt,nodes1,nodes2, &
          modep
	  double precision l,lhat,lindex
	  dimension y(2),yp(2)
	  common/splinq/vq(6)
	  common /surf/rsurf
	  call splntl(rt,vq)
!
	yp(1)=y(1)*(vq(2)-3.d0+lindex)+y(2)*(lhat*vq(1)/eigt-vq(2))
	yp(2)=y(1)*(eigt/vq(1)-vq(3))+y(2)*(1.d0-vq(4)+vq(3)+lindex)

	do i=1,2
		yp(i)=yp(i)*vq(5)
	enddo

	return

	end subroutine rkfcow
