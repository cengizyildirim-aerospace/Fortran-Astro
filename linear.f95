program linear 
    
    ! EXAMPLE USAGE OF LINEAR OPERATIONS IN ASTRO MODULE 

    use astro

    implicit none 
    real(kind(1.d0)),dimension(3) :: vector1 , vector2 , resultant_vector , ecc_v
    real(kind(1.d0)) :: dot_result, angle_dot,angle_dot_degrees , mag_vector1 , mag_vector2 , k_
    integer :: i 

    write(*,"(a)",advance="no") "Enter your vector1 as x y z : "
    read(*,*) vector1

    write(*,"(a)",advance="no") "Enter your vector2 as x y z : "
    read(*,*) vector2


    call cross_product3D(vector1,vector2,resultant_vector)

    call dot_product3D(vector1,vector2,dot_result)

    
    
    print "(a,f30.7,/)", "The xyz components of the cross product:" 
    
    do i = 1,3
        if (i /= 3) then
            write (*,"(f30.7)",advance="no")  resultant_vector(i)
        else 
            write (*,"(f30.7,/)")  resultant_vector(i)
        end if 

    end do 


    mag_vector1 = sqrt( vector1(1)**2 + vector1(2)**2 + vector1(3)**2 )
    mag_vector2 = sqrt( vector2(1)**2 + vector2(2)**2 + vector2(3)**2 )

    call angle_between_vectors(vector1,vector2,angle_dot)
    angle_dot_degrees = 180.0*angle_dot/3.14159265358979323

    print "(a,f30.7,/)",  "This is the result of the dot product of v1 and v2 : " , dot_result

    write (*,"(a,f30.7,a/)") "This is the angle between the vector1 and vector2 in  3D:" , angle_dot_degrees , " degrees"

    call init_planets()

    k_ = Earth%k

    call EccentricityVector(vector1,vector2,k_,ecc_v)


    print *, ecc_v

end program linear 

