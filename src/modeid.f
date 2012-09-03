! ***************************************************
! ***************************************************
	subroutine modeid
! this subroutine identifies the mode of the pulsation
! using a phase diagram method. 
! see unno et al. (1989).

	  use consts
	  
	  implicit double precision (a-h,o-z)
	  common g(200),rho(200),x(200),yeig(2,200)
	  double precision l,lhat,lindex
	  common/misc/l,lhat,lindex,nsurf,period, &
          eps,verg,eig,eigt,nodes1,nodes2, &
          modep
	  integer quad,quad2,quad1
!
	modep=0
	! first count the nodes in y1 and y2.
	nodes1=0
	nodes2=0
	do i=2,nsurf
		fit1=yeig(1,i)*yeig(1,i-1)
		fit2=yeig(2,i)*yeig(2,i-1)
		if (fit1.le.0) nodes1=nodes1+1
		if (fit2.le.0) nodes2=nodes2+1
	enddo

	! count crossings in the phase diagram.
	quad1=quad(yeig(1,1),yeig(2,1))
	do i=2,nsurf
		quad2=quad(yeig(1,i),yeig(2,i))
		idq=quad2-quad1
		if (idq.eq.0) cycle
		! see if quadrant chabge is a crossing in y1.
		if (iabs(idq).eq.3) go to 20
		if (quad1.eq.3.and.quad2.eq.2) go to 20
		if (quad1.eq.2.and.quad2.eq.3) go to 20
		! not a y1 crossing.
		quad1=quad2
		cycle
		! if crossing is clockwise, increment mode count.
		! if crossing is counterclockwise, decrement mode count.
   20	if (idq.eq.-3.or.idq.eq.1) kdmod=-1
		if (idq.eq.3.or.idq.eq.-1) kdmod=1
		modep=modep+kdmod
		quad1=quad2
		! get next point.
	enddo

	return

	end
