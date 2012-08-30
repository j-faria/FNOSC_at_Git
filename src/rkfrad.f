! ****************************************************
! ****************************************************
      subroutine rkfrad(rt,y,yp)
!
! this routine supplies the derivatives for use in rkf
! for the radial case.
!
      implicit double precision (a-h,o-z)
      common/misc/l,lhat,lindex,nsurf,period,grav, &
          pi,pi4,p43,eps,verg,eig,eigt,nodes1,nodes2, &
          modep
      double precision l,lhat,lindex
      dimension y(2),yp(2)
      common/splinq/vq(6)
!
      call splntl(rt,vq)
      yp(1)=-(3.0d0*y(1)+y(2)/vq(3))
      yp(2)=vq(4)*(4.0d0*y(1)+eigt*y(1)/vq(1)+y(2))
      do 10 i=1,2
   10      yp(i)=yp(i)*vq(5)
      return
      end
