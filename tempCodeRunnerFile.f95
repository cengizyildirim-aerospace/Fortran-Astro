program d 
    implicit none 
    real(kind(1.d0)) :: rp , ra , ecc,m,e,x

    ra = 69.82E9
    rp = 46E9

    call BiSectionMethod(m,e,x)
    m = 77.44
    e = 0.076
    
    print *, x 

    




    contains

    subroutine BiSectionMethod(m, e , x1)
        implicit none 
        real(kind(1.d0)) , intent(in) ::  m , e
        real(kind(1.d0)) , intent(out)  :: x1 
        real(kind(1.d0)) :: x, x2 
        real(kind(1.d0)) :: error = 1d-10 , equation1 , equation2,equation
        integer :: i , n = 47
        equation1 = x1 - ( e*sin(x1) ) - m 
        equation2 = x2 - ( e*sin(x2) ) - m
        equation = x - ( e*sin(x) ) - m 
        

        x = 0.0
        x1 = 60
        x2 = 90


        do i = 1,n
            equation1 = x1 - ( e*sin(x1) ) - m 
            equation2 = x2 - ( e*sin(x2) ) - m
            equation = x - ( e*sin(x) ) - m 

            if (equation1*equation2  < 0) then
                x = x1
                x1 = (x1+x2)/real(2)
            else 
                x2 = x1 
                x1 = (x1+x)/real(2)
            end if 
        end do 
        print *, x1

        return
    
    end subroutine BiSectionMethod



end program d 