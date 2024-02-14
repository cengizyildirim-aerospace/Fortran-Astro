    program orbital
        use astro 
        implicit none
        character(1) :: gm 
        real(kind(1.d0)) :: raan , ta , i , w , a ,e  , k, gp ,h_total , p  
        real(kind(1.d0)) , dimension(3) ::  r_vector_xyz , v_vector_xyz , h_vector, node_vector , e_vector 

        write(*,"(a)",advance="no") " Enter your r_vector as xyz and your v_vector as xyz : "
        read(*,*) r_vector_xyz , v_vector_xyz

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



        ! First calculate the three essential vectors!

        call EccentricityVector(r_vector_xyz,v_vector_xyz,gp,e_vector)
        call SpecificAngularMomentumVector(r_vector_xyz,v_vector_xyz,h_vector)
        call NodeVector(h_vector,node_vector)


        call inclination(h_vector,i)

        call LongitudeOfTheAscendingNode(node_vector,raan)

        call ArgumentOfPeriapsis(node_vector,e_vector,w)

        call TrueAnomaly(e_vector,r_vector_xyz,v_vector_xyz,ta)

        call dot_product3D(e_vector,e_vector,e)
        call dot_product3D(h_vector,h_vector,h_total)


        e = sqrt(e)
        h_total = sqrt(h_total)
 
        call SemiLatusRectum_h_k(h_total,gp,p)
        call SemiMajorAxis_p_e (p,e,a)

         


        write(*,*) " All of the Keplarian elements are as follows : " 

        write(*,"(a,f23.7,/)") "Inclination : " , i 
        write(*,"(a,f23.7,/)") "Longititude of the ascending node : " , raan 
        write(*,"(a,f23.7,/)") "Argument of the periapsis : " , w 
        write(*,"(a,f23.7,/)") "True Anomaly : " , ta 
        write(*,"(a,f23.7,/)") "SemiMajorAxis : " , a 
        write(*,"(a,f23.7,/)") "Eccentricity : " , e 

        write(*,*) " These elements could be used in NASA GMAT! "


    end program orbital 

        