!------------------------------------------------------
!
      subroutine rkfs(f,neqn,y,t,tout,relerr,abserr, &
          iflag,yp,h,f1,f2,f3,f4,f5,savre,savae,nfe, &
          kop,init,jflag,kflag)
      implicit double precision (a-h,o-z)
      logical hfaild,output
      dimension y(neqn),yp(neqn),f1(neqn),f2(neqn), &
          f3(neqn),f4(neqn),f5(neqn)
      external f
      data u26/2.d-13/ , remin/2.d-13/
      data maxnfe/3000/
      if (neqn .lt. 1) go to 10
      if ((relerr .lt. 0.d0)  .or.  (abserr .lt. 0.d0)) &
         go to 10
      mflag=iabs(iflag)
      if ((mflag .ge. 1) .and. (mflag .le. 7)) go to 20
   10 iflag=7
      return
   20 if (mflag .eq. 1) go to 50
      if (t .eq. tout) go to 10
      if(mflag .ne. 2) go to 25
      if (init .eq. 0) go to 45
      if (kflag .eq. 3) go to 40
      if ((kflag .eq. 4) .and.  (abserr .eq. 0.d0)) &
          go to 30
      if ((kflag .eq. 5)  .and. (relerr .le. savre) &
          .and. (abserr .le. savae)) go to 30
      go to 50
   25 if (iflag .eq. 3) go to 40
      if ((iflag .eq. 4) .and. (abserr .gt. 0.d0)) &
          go to 45
   30 write (6,1000) iflag,t
 1000 format (16h rkf says iflag=,i5,3h x=,1pd12.4)
      stop
   40 nfe=0
      if (mflag .eq. 2) go to 50
   45 iflag=jflag
   50 jflag=iflag
      kflag=0
      savre=relerr
      savae=abserr
      rer=dmax1(relerr,remin)
      dt=tout-t
      if (mflag .eq. 1) go to 60
      if (init .eq. 0) go to 65
      go to 80
   60 init=0
      kop=0
      a=t
      call f(a,y,yp)
      nfe=1
      if (t .ne. tout) go to 65
      iflag=2
      return
   65 init=1
      ymax=0.d0
      ypn=0.d0
      do 70 k=1,neqn
        ypn=dmax1(dabs(yp(k)),ypn)
70      ymax=dmax1(dabs(y(k)),ymax)
      etn=rer*ymax+abserr
      h=dabs(dt)
      if(etn.ge.ypn*h**5) go to 80
      h=dmax1((etn/ypn)**0.2d0,u26*dmax1(dabs(t),h))
80    h=dsign(h,dt)
      if (dabs(h).ge.dabs(dt)) kop=kop+1
      if (kop.ne.100) go to 85
      iflag=6
      return
85    if (dabs(dt).gt.u26*dabs(t))go to 95
      do 90 k=1,neqn
90      y(k)=y(k)+dt*yp(k)
      a=tout
      call f(a,y,yp)
      nfe=nfe+1
      go to 300
95    output=.false.
      scale=2.d0/rer
      ae=scale*abserr
100   hfaild=.false.
      hmin=u26*dabs(t)
      dt=tout-t
      if (dabs(dt).ge.2.d0*dabs(h)) go to 200
      if (dabs(dt).gt.dabs(h)/0.9d0) go to 150
      output=.true.
      h=dt
      go to 200
150   h=0.5d0*dt
200   if (nfe.le.maxnfe) go to 220
      iflag=3
      kflag=3
      return
220   call fehl(f,neqn,y,t,h,yp,f1,f2,f3,f4,f5,f1)
      nfe=nfe+5
      eeoet=0.d0
      do 250 k=1,neqn
        et=dabs(y(k))+dabs(f1(k))+ae
        if (et.gt.0.d0) go to 240
        iflag=4
        kflag=4
        return
240     ee=dabs((-2090.d0*yp(k)+(21970.d0*f3(k)-15048.d0 &
           *f4(k)))+(22528.d0*f2(k)-27360.d0*f5(k)))
250     eeoet=dmax1(eeoet,ee/et)
      esttol=dabs(h)*eeoet*scale/752400.d0
      if (esttol.le.1.d0) go to 260
      hfaild=.true.
      output=.false.
      s=0.1d0
      if (esttol.lt.59049.d0)s=0.9d0/esttol**0.2d0
      h=s*h
      if (dabs(h).gt.hmin) go to 200
      iflag=5
      kflag=5
      return
260   t=t+h
      do 270 k=1,neqn
270     y(k)=f1(k)
      a=t
      call f(a,y,yp)
      nfe=nfe+1
      if (hfaild) go to 290
      s=5.d0
      if (esttol.gt.1.889568d-04) s=0.9d0/esttol**0.2d0
      h=dsign(dmax1(s*dabs(h),hmin),h)
290   if (output) go to 300
      if (iflag.gt.0) go to 100
      iflag=-2
      return
300   t=tout
      iflag=2
      return
      end
