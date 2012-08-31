! ***************************************************
! ***************************************************
      subroutine consl(x,y,m,c)
! this subroutine sets up coefficient arrays for spline
! interpolation.
!
      implicit double precision (a-h,o-z)
      dimension c(4,m)
      dimension x(m),y(m)
      dimension a(200,3),d(200),b(200),z(200),p(200)
!
	mm=m-1
	do k=1,mm
		d(k)=x(k+1)-x(k)
      	p(k)=d(k)/6.d0
      	z(k)=(y(k+1)-y(k))/d(k)
	enddo
	
	do k=2,mm
      	b(k)=z(k)-z(k-1)
	enddo

	a(1,2)=-1.d0-d(1)/d(2)
	a(1,3)=d(1)/d(2)
	a(2,3)=p(2)-p(1)*a(1,3)
	a(2,2)=2.d0*(p(1)+p(2))-p(1)*a(1,2)
	a(2,3)=a(2,3)/a(2,2)
	b(2)=b(2)/a(2,2)

	do k=3,mm
		a(k,2)=2.d0*(p(k-1)+p(k))-p(k-1)*a(k-1,3)
		b(k)=b(k)-p(k-1)*b(k-1)
		a(k,3)=p(k)/a(k,2)
		b(k)=b(k)/a(k,2)
	enddo

	q=d(m-2)/d(m-1)
	a(m,1)=1.d0+q+a(m-2,3)
	a(m,2)=-q-a(m,1)*a(m-1,3)
	b(m)=b(m-2)-a(m,1)*b(m-1)
	z(m)=b(m)/a(m,2)
	mn=m-2

	do i=1,mn
		k=m-i
		z(k)=b(k)-a(k,3)*z(k+1)
	enddo
	
	z(1)=-a(1,2)*z(2)-a(1,3)*z(3)

	do k=1,mm
		q=1.d0/(6.d0*d(k))
		c(1,k)=z(k)*q
		c(2,k)=z(k+1)*q
		c(3,k)=y(k)/d(k)-z(k)*p(k)
		c(4,k)=y(k+1)/d(k)-z(k+1)*p(k)
	enddo

	return
	
	end subroutine consl
