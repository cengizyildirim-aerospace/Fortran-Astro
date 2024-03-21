program keptoxyz 
    use astro
    implicit none 
    character(1) :: gm 
    real(kind(1.q0)) :: e , a , i , w , ra , m , gp ,ta 
    real(kind(1.q0)), dimension(3) :: r_vector , v_vector_xyz

    


    write(*,"(/a)",advance="no") "What is the eccentricity of the object? "
    read(*,*) e

    write(*,"(/a)",advance="no") "What is the semimajor axis of the object in meters? "  
    read(*,*) a

    write(*,"(/a)",advance="no") "What is the inclination of the object in degrees? "
    read(*,*) i

    write(*,"(/a)",advance="no") "What is the argument of periapsis of the object in degrees? "
    read(*,*) w
    
    write(*,"(/a)",advance="no") "What is the right ascension of the ascending node (RAAN) of the object in degrees? "
    read(*,*) ra

    write(*,"(/a)",advance="no") "What is the true anomaly of the object in degrees? "
    read(*,*) ta
    

    m = MeanAnomaly(ta,e)
    print*, m 

    write(*,"(a/)") " You need to select your gravitational parameter here:"
    write(*,"(a/)") "DO NOT ENTER A NUMBER! You shall select your centeral body (e - the EARTH , s - the SUN, m - The MOON) "
    
    read (*,*) gm

    if (gm == "e") then 
        call init_planets()
        gp = Earth%k 
    else if (gm =="s" ) then
        call init_planets()
        gp = Sun%k
    else if (gm =="m" ) then 
        call init_planets()
        gp = Moon%k 
    else 
        write(*,*) "You need to enter a character and it can only be s for the Sun - e for the Earth and m for the Moon"
    end if 


    call KeplarianElementsToXYZ(a,e,m,i,w,ra,ta,r_vector)

    write(*,*) "x: " ,r_vector(1)
    write(*,*) "y: " ,r_vector(2) 
    write(*,*) "z: " ,r_vector(3)

    call KeplarianElementsToVxyz(a,e,m,i,w,ra,ta,gp,r_vector,v_vector_xyz)

    write(*,*) "vx: " ,v_vector_xyz(1)
    write(*,*) "vy: " ,v_vector_xyz(2) 
    write(*,*) "vz: " ,v_vector_xyz(3)

    write(*,"(a,f30.7)") "The distance from the main body to the orbiting object: " , &
    sqrt(r_vector(1)**2 + r_vector(2)**2 + r_vector(3)**2)

    


end program 

