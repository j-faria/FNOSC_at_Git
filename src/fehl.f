!-----------------------------------------------------
!
	subroutine fehl(f,neqn,y,t,h,yp,f1,f2,f3,f4,f5,s)
	  implicit double precision (a-h,o-z)
	  dimension y(neqn),yp(neqn),f1(neqn),f2(neqn), &
         f3(neqn),f4(neqn),f5(neqn),s(neqn)
	
	ch=0.25d0*h
      
	do k=1,neqn
		f5(k)=y(k)+ch*yp(k)
	enddo

	call f(t+0.25d0*h,f5,f1)
	ch=0.09375d0*h
	
	do k=1,neqn
		f5(k)=y(k)+ch*(yp(k)+3.d0*f1(k))
	enddo

	call f(t+0.375d0*h,f5,f2)
	ch=h/2197.d0

	do k=1,neqn
		f5(k)=y(k)+ch*(1932.d0*yp(k)+(7296.d0*f2(k) -7200.d0*f1(k)))
	enddo

	call f(t+12.d0/13.d0*h,f5,f3)
	ch=h/4104.d0

	do k=1,neqn
		f5(k)=y(k)+ch*((8341.d0*yp(k)-845.d0*f3(k)) + (29440.d0*f2(k)-32832.d0*f1(k)))
	enddo

	call f(t+h,f5,f4)
	ch=h/20520.d0

	do k=1,neqn
	    f1(k)=y(k)+ch*((-6080.d0*yp(k)+(9295.d0*f3(k) &
			-5643.d0*f4(k)))+(41040.d0*f1(k)-28352.d0 * f2(k)))
	enddo

	call f(t+0.5d0*h,f1,f5)
	ch=h/7618050.d0

	do k=1,neqn
		s(k)=y(k)+ch*((902880.d0*yp(k)+(3855735.d0 &
			*f3(k)-1371249.d0*f4(k)))+(3953664.d0*f2(k) &
			+277020.d0*f5(k)))
	enddo

	return

	end subroutine fehl     
