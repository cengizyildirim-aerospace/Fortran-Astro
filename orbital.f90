module orbital_elements

    type planet

        character(60) :: name
        real(kind(1.q0)) :: k , radius , mass 

    end type planet

    type (planet) :: Earth , Moon, Sun ,Mercury, Mars, Jupiter, Saturn,Neptune,Venus,Uranus
    real(kind(1.q0)),parameter ::  pi = 3.14159265358979323846264338327950288q0

    

        
    contains

    subroutine init_planets()
        ! Mean radius is used 
        ! k stands for gravitational parameter

        Earth%k = 3.98600435507q14
        Earth%radius = 6371.0084q3
        Earth%mass = 5.97217q24

        
        Moon%k = 4.902800118q12
        Moon%mass = 7.34578924831069q22
        Moon%radius = 1737.4q3

        
        Sun%k = 1.32712440041279419q20
        Sun%mass = 1.98840987131q30
        Sun%radius = 695700q3

        Mercury%k = 2.2031868551q13
        Mercury%mass = 0.330103q24
        Mercury%radius = 2439.4q3
        
        Venus%k = 3.24858592000q14
        Venus%mass = 4.86731q24
        Venus%radius =  6051.8q3

        Mars%k = 4.2828375816q13
        Mars%mass = 0.64169q24
        Mars%radius = 3389.50q3

        Jupiter%k = 1.26712764100000q17
        Jupiter%mass = 1898.125q24
        Jupiter%radius = 69911q3

        Saturn%k = 3.7940584841800q16
        Saturn%mass= 568.317q24
        Saturn%radius = 58232q3

        Uranus%k = 5.794556400000q15
        Uranus%mass = 86.8099q24
        Uranus%radius = 25362q3

        Neptune%k = 6.836527100580q15
        Neptune%mass = 102.4092q24
        Neptune%radius = 24622q3

        

    end subroutine init_planets



    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ---------------------------------- From Keplarian to x y z and vx vy vz ----------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------



    subroutine KeplarianElementsToXYZ(a,e,m,i,w,ra,ta,r_vector_xyz)   ! https://ntrs.nasa.gov/api/citations/19650015945/downloads/19650015945.pdf
        implicit none
        real(kind(1.q0)), intent(inout) :: a,e,m,i,w,ra ,ta
        real(kind(1.q0)),dimension(3), intent(out) :: r_vector_xyz 
        real(kind(1.q0)) :: ecc_anomaly 
        integer :: b 
        
        call init_planets()

        i = i*pi/180.0
        w = w*pi/180.0
        ra = ra*pi/180.0

        ecc_anomaly = EccentricAnomaly(ta,e)
        

        write(*,"(a , f30.7)") "Eccentric anomoly is : " , ecc_anomaly

        r_vector_xyz(1) = a * ((cos(ecc_anomaly)-e)*(cos(w)*cos(ra)-sin(w)*sin(ra)*cos(i))+ &
        sqrt(1-e**2) * sin(ecc_anomaly)* (-sin(w)*cos(ra)-cos(w)*sin(ra)*cos(i)))

        r_vector_xyz(2) = a * ((cos(ecc_anomaly)-e)*(cos(w)*sin(ra)+sin(w)*cos(ra)*cos(i))+  &
        sqrt(1-e**2) * sin(ecc_anomaly)* (-sin(w)*sin(ra)+cos(w)*cos(ra)*cos(i)))

        r_vector_xyz(3) = a * ((cos(ecc_anomaly)-e)*(sin(w)*sin(i))+  sqrt(1-e**2) * sin(ecc_anomaly)* (cos(w)*sin(i)))

    
        return

    end subroutine KeplarianElementsToXYZ




    subroutine KeplarianElementsToVxyz(a,e,m,i,w,ra,ta,k,r_vector_xyz,v_vector_xyz)

        implicit none 
        real(kind(1.q0)), intent(in) :: a,e,m,i,w,ra,ta, k
        real(kind(1.q0)) , dimension(3) , intent(in) :: r_vector_xyz
        real(kind(1.q0)) , dimension(3) , intent(out) :: v_vector_xyz
        real(kind(1.q0)) :: x_versor , y_versor, z_versor , v_total , gama , u, p , r_total

        

        call SemiLatusRectum_a_e(a,e,p)
        
        u = w + ta 
        r_total = sqrt(r_vector_xyz(1)**2 + r_vector_xyz(2)**2 + r_vector_xyz(3)**2)
        v_total = sqrt(k*(2/r_total -1/a))
        print *, v_total
        print *, r_total

        gama = asin(sqrt(k/p)*e*sin(ta)/v_total)
        print*, gama 

        x_versor = cos(ra) * cos(u) - sin(ra)*cos(i)*sin(u)
        y_versor = sin(ra)*cos(u) +cos(ra)*cos(i)*sin(u)
        z_versor = sin(i)*sin(u)

        v_vector_xyz(1) = v_total * (    x_versor*sin(gama) - cos(gama) * ( cos(ra)*sin(u)+sin(ra)*cos(i)*cos(u)  )   )
        v_vector_xyz(2) = v_total * ( y_versor*sin(gama) - cos(gama) * ( sin(ra)*sin(u)- cos(ra)*cos(i)*cos(u)  ) )
        v_vector_xyz(3) = v_total * ( z_versor*sin(gama)  + cos(gama)*cos(u)*sin(i)        )

        return

    

    end subroutine KeplarianElementsToVxyz



   

    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ---------------------------------- From x y z and vx vy vz to Keplarian ----------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------






    subroutine EccentricityVector(r_vector_xyz,v_vector_xyz,k,ecc_vector)
        real(kind(1.q0)) , dimension(3) ,intent(inout) :: r_vector_xyz , v_vector_xyz
        real(kind(1.q0)) , dimension(3) ,intent(out) :: ecc_vector
        real(kind(1.q0)) , dimension(3) :: h , j  ! specific angular momentum and  j= vxh  BATE, page 25 
        real(kind(1.q0)) , intent(in) :: k  ! gravitaonal parameter (GM)
        real(kind(1.q0)) :: v_total , r_total , v, r, s    !s = r * v

        call dot_product3D(v_vector_xyz,v_vector_xyz,v)

        call dot_product3D(r_vector_xyz,r_vector_xyz,r)

        call dot_product3D(r_vector_xyz,v_vector_xyz,s)

        call cross_product3D(r_vector_xyz,v_vector_xyz,h)

        call cross_product3D(v_vector_xyz,h,j)

        v_total = sqrt(v)
        r_total = sqrt(r)

        ecc_vector =  (j/k) - (r_vector_xyz/r_total)

        return 
        
    end subroutine EccentricityVector



    subroutine NodeVector(h,node_vector)  ! h is anglular momentum vector.
        implicit none
        real(kind(1.q0)) , dimension(3) ,intent(in) :: h
        real(kind(1.q0)) , dimension(3) ,intent(out) :: node_vector
        real(kind(1.q0)) , dimension(3)  :: k 

        k(1) = 0 
        k(2) = 0
        k(3) = 1

        call cross_product3D(k,h,node_vector)

        return
    end subroutine NodeVector

    subroutine SpecificAngularMomentumVector(r_vector_xyz,v_vector_xyz,h)
        implicit none
        real(kind(1.q0)) , dimension(3) ,intent(in) :: r_vector_xyz , v_vector_xyz
        real(kind(1.q0)) , dimension(3) ,intent(out) :: h 
    
        call cross_product3D(r_vector_xyz,v_vector_xyz,h)

        return
    
    end subroutine SpecificAngularMomentumVector

    ! -------------------------------------------------------------------------------------------------
    ! -----------------------------ORBITAL ELEMENTS (KEPLARIAN) ---------------------------------------
    ! -------------------------------------------------------------------------------------------------


    subroutine inclination(h,i)
        implicit none 
        real(kind(1.q0)) , dimension(3), intent(in) :: h 
        real(kind(1.q0)) , intent(out) :: i 
        real(kind(1.q0)) :: h_total 

        call dot_product3D(h,h,h_total)
        h_total = sqrt(h_total)
        i = acos(h(3)/ h_total)
        i = i*180.0/pi
        
        return


    end subroutine inclination


    subroutine LongitudeOfTheAscendingNode(node_vector,ascending_node)
        implicit none 
        real(kind(1.q0)) , dimension(3), intent(in) :: node_vector
        real(kind(1.q0)) , intent(out) :: ascending_node 
        real(kind(1.q0)) :: n_total 

        call dot_product3D(node_vector,node_vector,n_total)
        n_total = sqrt(n_total)
        ascending_node = acos(node_vector(1)/n_total)
        ascending_node = ascending_node*180/pi

        if (node_vector(2)<0) ascending_node = 360.0 - ascending_node

        
        return
    end subroutine LongitudeOfTheAscendingNode


    subroutine ArgumentOfPeriapsis(node_vector,EccentricityVector,w)
        real(kind(1.q0)) , dimension(3), intent(in) :: node_vector,EccentricityVector 
        real(kind(1.q0)) , intent(out) :: w
        real(kind(1.q0)) :: n_total , e_total , dot_ne 

        call dot_product3D(node_vector,node_vector,n_total)
        call dot_product3D(EccentricityVector,EccentricityVector,e_total)
        call dot_product3D(node_vector,EccentricityVector,dot_ne)


        n_total = sqrt(n_total)
        e_total = sqrt(e_total)

        w = acos(dot_ne/(n_total*e_total))
        w = w*180.0/pi

        if (EccentricityVector(3)<0) w = 360.0 - w 

        return
    
    end subroutine ArgumentOfPeriapsis


    subroutine TrueAnomaly(EccentricityVector,r_vector_xyz,v_vector_xyz,ta)
        implicit none 
        real(kind(1.q0)) , dimension(3), intent(in) :: r_vector_xyz,EccentricityVector ,v_vector_xyz
        real(kind(1.q0)) , intent(out) :: ta
        real(kind(1.q0)) :: r_total , e_total , dot_re, checker

        call dot_product3D(r_vector_xyz, r_vector_xyz ,r_total)
        call dot_product3D(EccentricityVector,EccentricityVector,e_total)

        call dot_product3D(EccentricityVector,r_vector_xyz,dot_re)
        call dot_product3D(r_vector_xyz,v_vector_xyz,checker)

        e_total =sqrt(e_total)
        r_total =sqrt(r_total)

        ta = acos(dot_re/(e_total*r_total))

        ta = ta*180.0/pi

        if (checker<0) ta =360.0 - ta 


        return
    
    end subroutine TrueAnomaly


 

    subroutine SemiMajorAxis_p_e(p,e,a)
        implicit none 
        real(kind(1.q0)), intent(in) :: p , e 
        real(kind(1.q0)), intent(out) :: a

        a = p / (1-e**2)

        return 

    end subroutine SemiMajorAxis_p_e


    subroutine PrecessionOfThePeriapsis(e,a,i,k,r,pop)   ! AVARAGE PRECESSION
        implicit none 
        real(kind(1.q0)), intent(in) ::  e , a ,k, r 
        real(kind(1.q0)) ,intent(out)::  pop
        real(kind(1.q0)), intent(inout) :: i
        real(kind(1.q0)) :: j2 

        j2= 1.08264q-3  

        i = i*pi/180.0

        pop = 5.0 * ((r/a)**3.5) * (5*((cos(i))**2) -1)/(1-e**2)**2

        

        return 

    end subroutine PrecessionOfThePeriapsis


    subroutine PrecessionOfTheNodeLine(e,a,i,k,r,pon) ! AVARAGE PRECESSION 
        implicit none 
        real(kind(1.q0)), intent(in) ::  e , a ,k, r 
        real(kind(1.q0)) ,intent(out)::  pon
        real(kind(1.q0)), intent(inout) :: i
        real(kind(1.q0)) :: j2

        call init_planets()

        j2= 1.08264q-3 

        i = i*pi/180.0

        pon = - (   (3.0/2.0)*sqrt(k)*j2*r**2 / (  ((1-e**2)**2) * a**3.5      )  * (cos(i))      )

        pon = pon*180.0/pi

        return 

    end subroutine PrecessionOfTheNodeLine


    
    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ---------------------------------- Other Orbital elements and more ---------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------


    ! -------------------------------------------------------------------------------------------------
    ! ----------------------------- Semi Latus Rectum -------------------------------------------------
    ! -------------------------------------------------------------------------------------------------


    subroutine SemiLatusRectum_a_e(a,e,p)
        implicit none 
        real(kind(1.q0)) , intent(in) :: a ,e 
        real(kind(1.q0)) , intent(out) :: p 

        p = a * (1-e**2)

        return 

    end subroutine SemiLatusRectum_a_e

    
    real function SemiLatusRectum(ra,rp)

        real(kind(1.q0)),intent(in) :: ra,rp 
        real(kind(1.q0)) :: sm , ecc

        sm = SemiMajorAxis(ra,rp)
        ecc = Eccentricity(ra,rp)
        SemiLatusRectum = sm * (1-ecc**2)
        return
    end function SemiLatusRectum
    

    subroutine SemiLatusRectum_h_k(h,k,p)
        implicit none 
        real(kind(1.q0)), intent(in) :: h , k 
        real(kind(1.q0)), intent(out) :: p

        p = h**2/k 
        
        return 
    end subroutine SemiLatusRectum_h_k





    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ---------------------------------------- LINEAR OPERATIONS  ----------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------





    subroutine cross_product3D(a,b,c)
        
        implicit none

        real(kind(1.q0)), dimension(3), intent(in) :: a,b  ! the vectors of a x b .. they should be placed in this order!
        real(kind(1.q0)), dimension(3), intent(out) :: c   ! the resulting vector 
        

        ! a(1) is x   , a(2) is y and a(3) is z  in xyz cartesian coordinate system.

        c(1) = a(2)*b(3) - b(2)*a(3)
        c(2) = -(a(1)*b(3) - b(1)*a(3))
        c(3) = a(1)*b(2) - b(1)*a(2)

        return 

    end subroutine cross_product3D


    subroutine dot_product3D(a,b,c)

        implicit none 

        real(kind(1.q0)), dimension(3), intent(in) :: a,b  ! the vectors of a * b .. they should be placed in this order!
        real(kind(1.q0)), intent(out) :: c   !the result 

        c = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

        return 
        
    end subroutine dot_product3D



    subroutine angle_between_vectors(v1,v2,angle)  ! It will print it in rad 
        implicit none
        real(kind(1.q0)) , dimension(3) , intent(in) :: v1 ,v2 
        real(kind(1.q0)) , intent(out) :: angle
        real(kind(1.q0)) :: dt_product , v1_m , v2_m

        call dot_product3D(v1,v2,dt_product)
        call dot_product3D(v1,v1,v1_m)
        call dot_product3D(v2,v2,v2_m)

        angle = acos(dt_product / (  sqrt(v1_m)*sqrt(v2_m)  )) 

        return
    end subroutine angle_between_vectors



end module orbital_elements
