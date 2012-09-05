! *****************************************************
	program puls
!******************************************************
! this program solves both the nonradial pulsation problem
! in the adiabatic cowling approximation and the
! purely radial adiabatic problem. we thank
! w. d. pesnell and p. jones for help in writing earlier
! versions of this code which solved the full fourth
! order system for adiabatic nonradial oscillations
! (plus quite a few other things).
! the input used is generated by the zams model builder.
! the method is to solve the two cowling approximation
! equations in their dziembowski form or the two
! purely radial equations by shooting from
! the center to the surface where, eventually and by
! iteration, the outer boundary conditions on the
! eigenfunctions are satisfied. the only guess you have
! to provide is that for the period (in seconds). you
! must also specify l which is angular momentum
! quantum number of the spherical harmonic function.
! enter a zero if you want the radial case.
! ******************************************************
	  use lib_array
	  use consts
	  use commonvar
	  
	  implicit double precision (a-h,o-z)
	  double precision l,lhat,lindex
	  common/misc/l,lhat,lindex,nsurf,period, &
          eps,verg,eig,eigt,nodes1,nodes2, &
          modep
	  real, dimension(16)     :: p   ! periods from obs file
	  integer                 :: iter = 1
	  character(len=20)       :: outputfile
	  character               :: ioeigfun
	
	
	
	!
	! read in model quantities.
	call read
	
	
	nfun = 2
	allocate( yeig(nfun,nsurf) )
	
	! do you want the eigenfunctions saved to file?
	write (6,2000)
 2000 format (' print out eigenfunctions (y/n)?')
	read (5,1100) ioeigfun
 1100 format (a1)
      if (ioeigfun.eq.'y') iprnt=1
      

	! read solar frequencies to set search
	p = 0.0	! initialize 
	call read_sun( p )
!	do i=1,size(p)
!		write(6,*) p(i)
!	enddo


   10 continue 
!	 write (6,2000)
! 2000 format (' enter guessed period (in seconds)')
!	 write (6,2001)
! 2001 format (' (enter a period of 0 to stop)')
!      read (5,*) period

	 period = p(iter)
	 iter = iter+1
	 if( period.eq.0.0d0 ) then
      	write (6,1000)
 1000	format (' calculation completed')
 		stop
	 else
	 	continue
	 end if

   20 write (11,1010) period
      write (6, 1010) period
 1010 format (' guessed period=',1pe14.6,4h sec)

	! eig is the first guess at the eigenvalue = sigma**2.
	eig = (2.0d0*pi/period)**2
	nconv = 0
	! maximum iterations allowed
	nserch=20



	! try to converge to a solution using newton's method.
	!
	do ntry=1,nserch
		eigt=eig
		if (radial) then
			call goout0(disc)
		else
			call gooutl(disc)
		endif
		!
		! perturb eigenvalue and find change in boundary condition.
		! disc is the boundary condition.
		!
		eigt=(1.d0+eps)*eig
		if (radial) then
			call goout0(dum1)
		else
			call gooutl(dum1)
		endif
		deigd=dum1-disc
		deig=-disc*eps*eig/deigd
		write (6,1060) ntry,eig,deig
		write (11,1060) ntry,eig,deig
		! deig is the predicted change in eig.
 1060 format(' iteration=',i3,' eigval=',1pe15.7,' diffeigval=', 1pe13.5)
		adeig=dabs(deig/eig)
		! check for convergence.
		if (adeig.lt.verg) nconv=1
		! compute new guess for eig.
		epseig=deig/eig
		if (adeig.gt.0.1d0) &
		     epseig=0.1d0*deig/dabs(deig)
		eig=eig*(1.d0+epseig)
		if (nconv.eq.1) go to 80
	enddo

	write (11,1080)
 1080 format (14h not converged)
	stop
	
	!
	! if converged.
	!
   80    eigt=eig
	if (radial) then
		call goout0(disc)
	else
		call gooutl(disc)
	endif
	
	period=(2.d0*pi)/dsqrt(eig)
	write ( 6,1090) period,eig
	write (11,1090) period,eig
 1090 format (' final values: period= ',1pe14.6, ' eig=',1pe17.9)
 
	! count nodes in y1.
	call modeid
	! in radial modes we introduce an additional (+1) 
	! cross in the phase diagram 
	if (radial) then
		write (6,1050) nodes1+1
		write (11,1050) nodes1+1
 1050 format (i3,' radial nodes in y1')
	else
		write (6,1110) nodes1,nodes2,modep
		write (11,1110) nodes1,nodes2,modep
 1110 format(i3,' nodes in y1   ',i3,' nodes in y2 '/'     phase diagram mode ',i3)
	endif

	!
	! output mode frequencies to file
	outputfile = 'freqs'
	open (unit=13, file=outputfile, status='unknown', action='write')
	call output(outputfile, int(l), nodes1, modep, period)
	
	
	if (iprnt.ne.1) go to 140
	
	! if user asked to print eigenfunctions
	ysurf=yeig(1,nsurf)

	do i=1,nsurf
		yeig(1,i)=yeig(1,i)/ysurf
		yeig(2,i)=yeig(2,i)/ysurf
	enddo

	write (11,1130)
	write (11,1140) (i+1,r(i)/r(nsurf),yeig(1,i), &
             yeig(2,i),i=1,nsurf)
 1130 format('                   2(n, r/r* ,  y1  ,  y2)')
 1140 format(2(i4,0pf10.6,2(1pe12.4)))
  
  140 continue
      
      go to 10
      
!      deallocate(yeig)
      deallocate(r, g, rho)
     
      end
