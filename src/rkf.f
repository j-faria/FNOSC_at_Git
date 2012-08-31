! ****************************************************
! ********************************************************
	subroutine rkf(f,neqn,y,t,tout,relerr,abserr, &
                    iflag,work,iwork)

	  implicit double precision (a-h,o-z)
	  dimension y(neqn),work(1),iwork(5)
	  external f

	k1m=neqn+1
	k1=k1m+1
	k2=k1+neqn
	k3=k2+neqn
	k4=k3+neqn
	k5=k4+neqn
	k6=k5+neqn
	call rkfs(f,neqn,y,t,tout,relerr,abserr,iflag, &
               work(1),work(k1m),work(k1),work(k2),work(k3), &
               work(k4),work(k5),work(k6),work(k6+1), &
               iwork(1),iwork(2),iwork(3),iwork(4),iwork(5))

	return

	end subroutine rkf
