subroutine simple_fort_rhs(blockno, mbc,mx,my,mfields, & 
                           xlower,ylower,dx,dy,rhs)
    IMPLICIT NONE

    INTEGER :: mbc,mx,my, mfields, blockno
    DOUBLE PRECISION :: xlower,ylower,dx,dy
    DOUBLE PRECISION :: rhs(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)

    INTEGER :: i,j, m
    DOUBLE PRECISION :: xc,yc, simple_qexact_rhs
    double precision :: qlap

    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            qlap = simple_qexact_rhs(xc,yc)
            do m = 1,mfields
                rhs(i,j,m) =  qlap
            end do
        end do
    end do

end subroutine simple_fort_rhs


