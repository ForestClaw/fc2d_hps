!! # check to see if value exceeds threshold

integer function user_exceeds_th(blockno,& 
                                  qval,qmin,qmax,quad, & 
                                  dx,dy,xc,yc,threshold, &
                                  init_flag, is_ghost)
    implicit none
    
    double precision :: qval,qmin,qmax,threshold
    double precision :: quad(-1:1,-1:1)
    double precision :: dx,dy, xc, yc
    integer :: blockno, init_flag
    logical(kind=4) :: is_ghost
    integer :: refine

    refine = 0
    if (xc .gt. 0 .and. yc .gt. 0) then
        refine = 1
    endif

    user_exceeds_th = refine

end function user_exceeds_th
