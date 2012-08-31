! ****************************************************
! ****************************************************
	subroutine splntl(xint,yout)
! this subroutine interpolates model quantities for
! intermediate steps in the integration using cubic
! splines.
	  implicit double precision (a-h,o-z)
	  common/savngl/k,kq
	  common /savel/m,x(200),c1(4,200),y1(200),c2(4,200), &
         y2(200),c3(4,200),y3(200),c4(4,200),y4(200), &
         c5(4,200),y5(200),c6(4,200)
	  common/rs/y6(200)
	  dimension yout(6)
!
	mm=m-1
	if (xint.ge.x(1).and.xint.lt.x(m)) go to 20
	if (xint .lt. (1.+1.d-8)*x(1)) go to 60
	if (xint .ge. x(m)) go to 10
	k=1
	go to 50
	
   10 if (xint .gt. (1.+1.d-6)*x(m)) go to 60
      k=mm
      go to 50
      
   20 il=1
      ir=m
      
   30 k=il+((ir-il)/2)
      if (xint.ge.x(k)) go to 40
      ir=k
      go to 30
      
   40 if (xint.lt.x(k+1)) go to 50
      il=k
      go to 30
      
   50 continue
   
	x1=x(k+1)-xint
	xx=xint-x(k)
	x12=x1*x1
	xx2=xx*xx
	yout(1)=x1*(c1(1,k)*x12 +c1(3,k))+xx*(c1(2,k)*xx2+c1(4,k))
	yout(2)=x1*(c2(1,k)*x12 +c2(3,k))+xx*(c2(2,k)*xx2+c2(4,k))
	yout(3)=x1*(c3(1,k)*x12 +c3(3,k))+xx*(c3(2,k)*xx2+c3(4,k))
	yout(4)=x1*(c4(1,k)*x12 +c4(3,k))+xx*(c4(2,k)*xx2+c4(4,k))
	yout(5)=x1*(c5(1,k)*x12 +c5(3,k))+xx*(c5(2,k)*xx2+c5(4,k))
	yout(6)=x1*(c6(1,k)*x12 +c6(3,k))+xx*(c6(2,k)*xx2+c6(4,k))
	
	return
	
   60 continue
	 write (16,1000) xint
 1000 format (1h0,25hrange error in spline, x=,e15.6)
 
	stop
	
	end subroutine splntl
