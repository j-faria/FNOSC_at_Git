! *****************************************************
! *****************************************************
	subroutine read
!  this routine reads in model quantities.

	  use consts
	  use commonvar
	  
	  implicit double precision (a-h,o-z)
       common/misc/l,lhat,lindex,nsurf,period, &
          eps,verg,eig,eigt,nodes1,nodes2, &
          modep
	  common /savel/m
	  common /surf/rsurf
	  
	  dimension xc(1800), y(1800), c(4,1800)

	  
      
	double precision l,lhat,lindex
	character(len=20)    :: outfile, infile
!	

	
!
      radial=.false.
      !
      ! ask for input and output files
      write (6,*) ' your input filename is'
      read (5,3000) infile
 3000 format (a20)
      open (10,file=infile,status='unknown')
      write (6,*)  ' your output filename is '
      read (5,3000) outfile
      open (11,file=outfile,status='unknown')
      !
      ! ask for l
      write (6,1080)
 1080 format (' enter l')
      read (5,*) l
      
      lhat=l*(l+1.d0)
      lindex=2.0d0-l
      
	! setting up convergence criteria and eps for
	! newton's method.
	verg=1.00d-05
	eps=1.00d-03

	! read globals from model file
	read (10,*) amass	! mass of model
	read (10,*) xx,yy	! abundances
	write (11,1000) amass
!      write (6,1000) amass
 1000 format (' mass (in msun)=',0pf7.3)
	write (11,1001) xx,yy
!      write (6,1001) xx,yy
 1001 format (' x=',0pf6.3,' y=', f6.3)
      write (11,1090) l
 1090 format (3h l=,0pf7.0)
 
	! check to see if a purely radial calculation
	! is to be done.
	if (l .lt. 0.9d0) radial=.true.

	! read in model quantities from zams.
	read (10,*) nsurf
!	write(6,*) nsurf
	
	! allocate arrays
	allocate(r(nsurf), g(nsurf), rho(nsurf) ) 
	allocate( xi(nsurf), y1(nsurf), y2(nsurf), y3(nsurf), y4(nsurf), y5(nsurf) )

	
	read (10,*) (x(i),r(i),g(i),rho(i),i=1,nsurf)
	rsurf=r(nsurf)
	read (10,*) m
	do i=1,m
		read (10,*) xi(i),y1(i),y2(i),y3(i)
		read (10,*) y4(i),y5(i)
	enddo  
	
	

	if (radial) then
		! for the radial case y4 contains v (of the u-v plane)
		! and y3 contains gamma1.
		do i=1,m
			y4(i)=1.0d0/y5(i)-1.0d0
			y3(i)=y4(i)/y2(i)
		enddo
	endif
	
! set up spline coefficients for use with integrator.
	nc = 4
	allocate( c1(nc,nsurf), c2(nc,nsurf), c3(nc,nsurf), c4(nc,nsurf), c5(nc,nsurf), c6(nc,nsurf) )

	
	xc  = xi; y = y1; c = c1
	call consl(xc,y,m,c)
	c1 = c
	xc = xi; y = y2; c = c2
	call consl(xc,y,m,c)
	c2 = c
	xc = xi; y = y3; c = c3
	call consl(xc,y,m,c)
	c3 = c	
	xc = xi; y = y4; c = c4
	call consl(xc,y,m,c)
	c4 = c
	xc = xi; y = y5; c = c5
	call consl(xc,y,m,c)
	c5 = c
	xc = xi; y = r; c = c6
	call consl(xc,y,m,c)
	c6 = c
	

      return
      end


