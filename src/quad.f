! ****************************************************
! ****************************************************
	function quad(y1,y2)
! identifies the quadrant in the phase diagram
! used in modeid.
	  implicit double precision (a-h,o-z)
	  integer quad
!
	if (y1.gt.0.and.y2.ge.0) quad=1
	if (y1.le.0.and.y2.gt.0) quad=2
	if (y1.lt.0.and.y2.le.0) quad=3
	if (y1.ge.0.and.y2.lt.0) quad=4

	return
	
	end function quad
