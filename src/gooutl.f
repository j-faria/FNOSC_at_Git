! ******************************************************
! ******************************************************
	subroutine gooutl(b2)
! this controls the integration for the cowling
! approximation.

	  use consts
	  
	  implicit double precision (a-h,o-z)
	  common g(200),rho(200),x(200),yeig(2,200)
	  common/misc/l,lhat,lindex,nsurf,period, &
          eps,verg,eig,eigt,nodes1,nodes2, &
          modep
	  common/rs/r(200)
	  double precision l,lhat,lindex
	  common/splinq/vq(6)
	  dimension y(2),work(50),iwork(6)
	  external rkfcow
!
	yeig(1,1)=1.0d0
	yeig(2,1)=eigt/(l*grav*pi43*rho(1))
	y(1)=yeig(1,1)
	y(2)=yeig(2,1)
	iflag=1
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
		call rkf(rkfcow,neqn,y,xstart,xfin,relerr,abserr, &
              iflag,work,iwork)
		if (iflag.eq.3.or.iflag.eq.4.or.iflag.eq.5) &
                    write (11,1000) iflag,xfin
		iflag=1
		do i=1,2
			yeig(i,k)=y(i)
		enddo
	enddo

	rt=x(nsurf)
	call splntl(rt,vq)
	temp=1.d0/vq(5)-1.d0
	
	! the following is the outer boundary condition. we
	! iterate until b2 is close enough to zero that
	! corrections to the eigenvalue are within a small
	! tolerance.
	b2=yeig(1,nsurf)*((4.d0+r(nsurf)*eigt/g(nsurf)) &
          /temp-1.d0)+(1.d0-lhat*g(nsurf)/r(nsurf) &
          /eigt/temp)*yeig(2,nsurf)
          
 1000 format (17h watch out,iflag=,i4,6h xfin=,1pe10.2)

	return

	end subroutine gooutl
