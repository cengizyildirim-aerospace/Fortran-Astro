program precesion 
    use astro 

    implicit none 
    
    character(1) :: gm 
    real(kind(1.q0)), dimension(3) :: r_vector_xyz , v_vector_xyz ,fr_vector,fv_vector, h_vector ,e_vector, node_vector
    real(kind(1.q0)), dimension(3) :: e_vector_new
    real(kind(1.q0)) :: gp , step,i,r_total, pop ,e , h , p ,a , diff,w , new_w,r_earth
    integer :: fake_step,wr,m,finalt
    real(kind(1.q0)) :: raan , ta ,h_total  , e_ang
    



    write(*,*) " You need to enter your r_vector and v_vector in xyz :  " 
    read (*,*) r_vector_xyz , v_vector_xyz


    write(*,*) " You need to enter your step size and final time at which you want calculation to stop: " 
    read (*,*) step , finalt 

    write(*,"(a/)") " You need to select your gravitational parameter here:"
    write(*,"(a/)") "DO NOT ENTER A NUMBER! You shall select your centeral body (e - the EARTH , s - the SUN, m - The MOON) "
    
    read (*,*) gm

    if (gm == "e") then 
        call init_planets()
        gp = Earth%k 
        r_earth = 6378.1366d3
    else if (gm =="s" ) then
        call init_planets()
        gp = Sun%k
    else if (gm =="m" ) then 
        call init_planets()
        gp = Moon%k 
    else 
        write(*,*) "You need to enter a character and it can only be s for the Sun - e for the Earth and m for the Moon"
    end if 

    !Precession calculation

    
    call EccentricityVector(r_vector_xyz,v_vector_xyz,gp,e_vector)
    call SpecificAngularMomentumVector(r_vector_xyz,v_vector_xyz,h_vector)
    call NodeVector(h_vector,node_vector)
    call inclination(h_vector,i)
    call ArgumentOfPeriapsis(node_vector,e_vector,w)
    call dot_product3D(e_vector,e_vector,e)
    call dot_product3D(h_vector,h_vector,h_total)


    e = sqrt(e)
    h_total = sqrt(h_total)
 
    call SemiLatusRectum_h_k(h_total,gp,p)
    call SemiMajorAxis_p_e (p,e,a)
    
    
    call ArgumentOfPeriapsis(node_vector,e_vector,w)

    call PrecessionOfThePeriapsis(e,a,i,gp,r_earth,pop)

    

    

    write(*,"(a,f30.13)") " The predicted precession of the periapsis : " , pop



    call LeapFrogPerturbed3D(r_vector_xyz,v_vector_xyz,step,finalt,gp,fr_vector,fv_vector)

    r_vector_xyz(1) = fr_vector(1)
    r_vector_xyz(2) = fr_vector(2)
    r_vector_xyz(3) = fr_vector(3)


    v_vector_xyz(1) = fv_vector(1)      
    v_vector_xyz(2) = fv_vector(2)
    v_vector_xyz(3) = fv_vector(3)



    fake_step= finalt/step 

    do m  = 1,fake_step
        call RLeapFrogPerturbed3D(r_vector_xyz,v_vector_xyz,step,finalt,gp,fr_vector,fv_vector)  ! After the first step size, each element will have the same step size. So we use Rleapfrog.

        r_vector_xyz(1) = fr_vector(1)
        r_vector_xyz(2) = fr_vector(2)
        r_vector_xyz(3) = fr_vector(3)


        v_vector_xyz(1) = fv_vector(1)      
        v_vector_xyz(2) = fv_vector(2)
        v_vector_xyz(3) = fv_vector(3) 

    end do 


    call EccentricityVector(r_vector_xyz,v_vector_xyz,gp,e_vector_new)
    call SpecificAngularMomentumVector(r_vector_xyz,v_vector_xyz,h_vector)
    call NodeVector(h_vector,node_vector)

    call inclination(h_vector,i)


    call dot_product3D(e_vector_new,e_vector_new,e)
    call dot_product3D(h_vector,h_vector,h_total)


    e = sqrt(e)
    h_total = sqrt(h_total)
 
    call SemiLatusRectum_h_k(h_total,gp,p)
    call SemiMajorAxis_p_e (p,e,a)
    
    

    call ArgumentOfPeriapsis(node_vector,e_vector_new,new_w)

        

    diff =(new_w - w) 

    write(*,"(/a,f30.16)") " The actual precession of the periapsis : " , diff

    call angle_between_vectors(e_vector,e_vector_new,e_ang)

    e_ang = e_ang*180.0/3.14159265358979323
    print *, e_ang
    

end program precesion