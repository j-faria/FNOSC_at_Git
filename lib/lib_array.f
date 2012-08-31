! MD5 of template: 874b367ed9a2602d10f066584f750d00
! Array related routines (Integration, Interpolation, etc.)
! Thomas Robitaille (c) 2009

module lib_array

  implicit none
  save


  integer,parameter :: idp = selected_int_kind(13)
  integer,parameter :: sp = selected_real_kind(p=6,r=37)
  integer,parameter :: dp = selected_real_kind(p=15,r=307)


contains



  subroutine linspace_dp(xmin,xmax,x)
    implicit none
    real(dp),intent(in) :: xmin,xmax
    real(dp),intent(out) :: x(:)
    integer :: i,n
    n = size(x)
    if (n == 1) then
       if(xmin /= xmax) then
          write(0,'("ERROR: Cannot call linspace with n=1 and xmin /= xmax")')
          stop
       else
          x = xmin
       end if
    else
       do i=1,n
          x(i) = (xmax-xmin) * real(i-1,dp) / real(n-1,dp) + xmin
       end do
    end if
  end subroutine linspace_dp


  subroutine logspace_dp(xmin,xmax,x)
    implicit none
    real(dp),intent(in) :: xmin,xmax
    real(dp),intent(out) :: x(:)
    if (size(x) == 1 .and. xmin /= xmax) then
       write(0,'("ERROR: Cannot call logspace with n=1 and xmin /= xmax")')
       stop
    end if
    call linspace_dp(log10(xmin),log10(xmax),x)
    x = 10._dp**x
  end subroutine logspace_dp


  recursive subroutine quicksort_dp(array, left, right)
    implicit none
    real(dp),intent(inout) :: array(:)
    integer,intent(in) :: left, right
    integer :: pivot_index, new_pivot_index
    if(right > left) then
       pivot_index = left+(right-left)/2
       call partition_dp(array, left, right, pivot_index, new_pivot_index)
       call quicksort_dp(array, left, new_pivot_index - 1)
       call quicksort_dp(array, new_pivot_index + 1, right)
    end if
  end subroutine quicksort_dp

  recursive subroutine quicksort_all_dp(array)
    implicit none
    real(dp),intent(inout) :: array(:)
    call quicksort_dp(array, 1, size(array))
  end subroutine quicksort_all_dp

  subroutine swap_dp(array, i, j)
    implicit none
    real(dp),intent(inout) :: array(:)
    integer,intent(in) :: i, j
    real(dp) :: temp
    temp = array(j)
    array(j) = array(i)
    array(i) = temp
  end subroutine swap_dp

  subroutine partition_dp(array, left, right, pivot_index, store_index)
    implicit none
    real(dp),intent(inout) :: array(:)
    integer,intent(in) :: left, right, pivot_index
    integer,intent(out) :: store_index
    real(dp) :: pivot_value
    integer :: i
    pivot_value = array(pivot_index)
    call swap_dp(array, pivot_index, right)
    store_index = left
    do i=left, right-1
       if(array(i) <= pivot_value) then
          call swap_dp(array, i, store_index)
          store_index = store_index + 1
       end if
    end do
    call swap_dp(array, store_index, right)
  end subroutine partition_dp


end module lib_array
