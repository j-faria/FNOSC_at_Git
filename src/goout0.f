! ***************************************************
! ***************************************************
	subroutine goout0(b2)
! this controls the integration for purely radial
! pulsations. It returns the value of the outer 
! boundary condition in B2

	  use consts
	  use commonvar
	  
	  implicit double precision (a-h,o-z)
	  common/misc/l,lhat,lindex,nsurf,period, &
          eps,verg,eig,eigt,nodes1,nodes2, &
          modep
	  double precision l,lhat,lindex
	  common/splinq/vq(6)
	  dimension y(2),work(50),iwork(6)
	  external rkfrad
!

	yeig(1,1) = 1.0d0
	rt=x(1)

	call splntl(rt,vq)
	
	
	! boundary condition at the center 
	yeig(2,1)=-3.0d0*vq(3)*yeig(1,1)	! eq 8.13
	y(1)=yeig(1,1)
	y(2)=yeig(2,1)
	iflag=1
	! number of equations
	neqn=2
	nsurf1=nsurf-1
	
	do j=1,nsurf1
		k=j+1
		xstart=x(j)
		xfin=x(k)
		relerr=1.d-8
		abserr=dabs(y(1))
		do i=2,neqn
			temp=dabs(y(i))
			abserr=dmin1(abserr,temp)
		enddo

		abserr=1.d-8*abserr
		abserr=dmax1(abserr,1.d-20)
		! integrate the eqs in RKFRAD from XSTART to XFIN
		! solution is returned in Y
		call rkf(rkfrad,neqn,y,xstart,xfin,relerr,abserr, &
				iflag,work,iwork)
		if (iflag.eq.3.or.iflag.eq.4.or.iflag.eq.5) &
                    write (11,1000) iflag,xfin
		iflag=1
		! update solution
		do i=1,2
			yeig(i,k)=y(i)
	   	enddo
	enddo

	rt=x(nsurf)
	call splntl(rt,vq)
	
	
	! the following is the outer boundary condition. 
	! we iterate until b2 is close enough to zero that
	! corrections to the eigenvalue are within a small
	! tolerance.
	b2=yeig(1,nsurf)*(4.0d0+eigt*r(nsurf)/g(nsurf))+yeig(2,nsurf)	! eq 8.14

 1000 format (17h watch out,iflag=,i4,6h xfin=,1pe10.2)

	return
	
	end subroutine goout0
