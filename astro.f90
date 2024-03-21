module astro

    

    type planet

        character(60) :: name
        real(kind(1.q0)) :: k , radius , mass 

    end type planet

    type (planet) :: Earth , Moon, Sun ,Mercury, Mars, Jupiter, Saturn,Neptune,Venus,Uranus,Pluto, &
    Io , Europa , Ganymede, Callisto, Mimas , Tethys, Dione , Rhea , Titan , Iapetus, Triton
    real(kind(1.q0)),parameter ::  pi = 3.14159265358979323846264338327950288q0

    

        
    contains

    subroutine init_planets()
        ! Mean radius is used 
        ! k stands for gravitational parameter

        Earth%k = 3.98600435436q14
        Earth%radius = 6371.0084q3
        Earth%mass = 5.97217q24

        
        Moon%k = 4.902800066q12
        Moon%mass = 7.34578924831069q22
        Moon%radius = 1737.4q3

        
        Sun%k = 1.3271244004193938q20
        Sun%mass = 1.98840987131q30
        Sun%radius = 695700q3

        Mercury%k = 2.2031868551q13
        Mercury%mass = 0.330103q24
        Mercury%radius = 2439.4q3
        
        Venus%k = 3.24858592000q14
        Venus%mass = 4.86731q24
        Venus%radius =  6051.8q3

        Mars%k = 4.2828375214q13
        Mars%mass = 0.64169q24
        Mars%radius = 3389.50q3

        Jupiter%k = 1.26686531900q17
        Jupiter%mass = 1898.125q24
        Jupiter%radius = 69911q3

        Saturn%k = 3.7931206234q16
        Saturn%mass= 568.317q24
        Saturn%radius = 58232q3

        Uranus%k = 5.793951256q15
        Uranus%mass = 86.8099q24
        Uranus%radius = 25362q3

        Neptune%k = 6.83509997q15
        Neptune%mass = 102.4092q24
        Neptune%radius = 24622q3

        Pluto%k = 8.6961q11 
        Pluto%mass = 1.307q22
        Pluto%radius = 1188.3q3

        Io%k = 5.9599155q12
        Io%mass = 8.9296488q22
        Io%radius = 1188.3q3

        Europa%k = 3.2027121q12
        Europa%mass = 8.9296488q22
        Europa%radius = 1560.8q3

        Ganymede%k = 9.8878328q12
        Ganymede%mass = 8.9296488q22
        Ganymede%radius = 2631.2q3

        Callisto%k = 7.1792834q12
        Callisto%mass = 8.9296488q22
        Callisto%radius = 2410.3q3

        Mimas%k = 2.503489q9
        Mimas%mass = 8.9296488q22
        Mimas%radius = 198.8q3 

        Tethys%k = 41.21q9
        Tethys%mass = 8.9296488q22
        Tethys%radius = 536.3q3 
  
        Dione%k = 73.116q9
        Dione%mass = 8.9296488q22
        Dione%radius = 562.5q3 

        Rhea%k = 153.94q9
        Rhea%mass = 8.9296488q22
        Rhea%radius = 764.5q3  

        Titan%k = 8978.14q9
        Titan%mass = 8.9296488q22
        Titan%radius = 2575.5q3   

        Iapetus%k = 120.52q9
        Iapetus%mass = 8.9296488q22
        Iapetus%radius = 734.5q3   

        Triton%k = 1428.495q9
        Triton%mass = 8.9296488q22
        Triton%radius = 1352.6q3  


        
    end subroutine init_planets





    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! --------------------------------------  LEAP FROG IN 2D (No Pertubations)  -------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------

   

    subroutine LeapFrog2D(x,vx0,y,vy0,step,finalt,k,fx,fy,fvx,fvy)
        implicit none 
        integer , intent(in) :: finalt
        real(kind(1.q0)) , intent(inout) :: vx0,vy0 ,step , k 
        real(kind(1.q0)) , intent(inout) :: x ,y 
        real(kind(1.q0)) , intent(out) :: fx , fy ,fvx,fvy
        real(kind(1.q0)) :: ax , ay  , r  



        

        r = sqrt(x**2+y**2)

        ax = -k * x/(r**3)
        ay = -k *y/(r**3)

        fvx = vx0 + ax*step/real(2)
        fvy = vy0 + ay*step/real(2)

        fx = x + fvx *step 
        fy = y +fvy*step     

        return
    end subroutine LeapFrog2D 


    ! ----------------------------------------------------------------------------------------------


    subroutine RLeapFrog2D(x,vx0,y,vy0,step,finalt,k,fx,fy,fvx,fvy)
        implicit none
        integer , intent(in) :: finalt 
        real(kind(1.q0)) , intent(inout) :: vx0,vy0 ,step , k 
        real(kind(1.q0)) , intent(inout) :: x ,y 
        real(kind(1.q0)) , intent(out) :: fx , fy ,fvx,fvy
        real(kind(1.q0)) :: ax , ay  , r  



        

        r = sqrt(x**2+y**2)

        ax = -k * x/(r**3)
        ay = -k *y/(r**3)

        fvx = vx0 + ax*step
        fvy = vy0 + ay*step

        fx = x + fvx *step 
        fy = y +fvy*step     

        return
    end subroutine RLeapFrog2D  


    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! --------------------------------------  LEAP FROG IN 3D (No Pertubations)  -------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------

    subroutine LeapFrog3D(r_vector_xyz,v_vector_xyz,step,finalt,k,fr_vector,fv_vector)
        implicit none 
        integer , intent(in) :: finalt
        real(kind(1.q0)), intent(in) :: step , k
        real(kind(1.q0)), dimension(3), intent(inout) :: v_vector_xyz  , r_vector_xyz
        real(kind(1.q0)), dimension(3), intent(out) :: fr_vector ,fv_vector
        real(kind(1.q0)), dimension(3) :: a  
        real(kind(1.q0)) :: r 


        call dot_product3D(r_vector_xyz,r_vector_xyz,r)

        r = sqrt(r)

        a(1) = -k * r_vector_xyz(1)/(r**3)
        a(2) = -k *r_vector_xyz(2)/(r**3)
        a(3) = -k *r_vector_xyz(3)/(r**3)



        fv_vector(1) = v_vector_xyz(1) + a(1)*step/real(2)
        fv_vector(2) = v_vector_xyz(2) + a(2)*step/real(2)
        fv_vector(3) = v_vector_xyz(3) + a(3)*step/real(2)

        fr_vector(1) = r_vector_xyz(1) + fv_vector(1)*step 
        fr_vector(2) = r_vector_xyz(2) + fv_vector(2)*step 
        fr_vector(3) = r_vector_xyz(3) + fv_vector(3)*step       

        return
    end subroutine LeapFrog3D 


    subroutine RLeapFrog3D(r_vector_xyz,v_vector_xyz,step,finalt,k,fr_vector,fv_vector)
        implicit none 
        integer , intent(in) :: finalt
        real(kind(1.q0)), intent(in) :: step , k
        real(kind(1.q0)), dimension(3), intent(inout) :: v_vector_xyz  , r_vector_xyz
        real(kind(1.q0)), dimension(3), intent(out) :: fr_vector ,fv_vector
        real(kind(1.q0)), dimension(3) :: a  
        real(kind(1.q0)) :: r 


        call dot_product3D(r_vector_xyz,r_vector_xyz,r)

        r = sqrt(r)

        a(1) = -k * r_vector_xyz(1)/(r**3)
        a(2) = -k *r_vector_xyz(2)/(r**3)
        a(3) = -k *r_vector_xyz(3)/(r**3)



        fv_vector(1) = v_vector_xyz(1) + a(1)*step
        fv_vector(2) = v_vector_xyz(2) + a(2)*step
        fv_vector(3) = v_vector_xyz(3) + a(3)*step

        fr_vector(1) = r_vector_xyz(1) + fv_vector(1)*step 
        fr_vector(2) = r_vector_xyz(2) + fv_vector(2)*step 
        fr_vector(3) = r_vector_xyz(3) + fv_vector(3)*step       

        return
    end subroutine RLeapFrog3D 







    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------


    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! --------------------------------------  LEAP FROG IN 3D (with Pertubation)  ------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------

    subroutine LeapFrogPerturbed3D(r_vector_xyz,v_vector_xyz,step,finalt,k,fr_vector,fv_vector)
        implicit none 
        integer , intent(in) :: finalt
        real(kind(1.q0)), intent(in) :: step , k
        real(kind(1.q0)), dimension(3), intent(inout) :: v_vector_xyz  , r_vector_xyz
        real(kind(1.q0)), dimension(3), intent(out) :: fr_vector ,fv_vector
        real(kind(1.q0)), dimension(3) :: a  
        real(kind(1.q0)) :: r
        real(kind(1.q0)) :: re ,j2,j3

        re =6378136.6  
        j2= 1082.64d-6 
        j3 = -2.5d-6

        call dot_product3D(r_vector_xyz,r_vector_xyz,r)

        r = sqrt(r)

        a(1) = -k *( r_vector_xyz(1)/(r**3) ) * ( 1 - j2*(3.0/2.0)*(re/r)**2 * (5*(r_vector_xyz(3)**2/r**2)-1)  )

        a(2) = -k *r_vector_xyz(2)/(r**3)* ( 1 - j2*(3.0/2.0)*(re/r)**2 * (5*(r_vector_xyz(3)**2/r**2)-1)  )

        a(3) = -k *r_vector_xyz(3)/(r**3) * ( 1 + j2*(3.0/2.0)*(re/r)**3 * (3- 5*(r_vector_xyz(3)**2/r**2)) ) 



        fv_vector(1) = v_vector_xyz(1) + a(1)*step/real(2)
        fv_vector(2) = v_vector_xyz(2) + a(2)*step/real(2)
        fv_vector(3) = v_vector_xyz(3) + a(3)*step/real(2)

        fr_vector(1) = r_vector_xyz(1) + fv_vector(1)*step 
        fr_vector(2) = r_vector_xyz(2) + fv_vector(2)*step 
        fr_vector(3) = r_vector_xyz(3) + fv_vector(3)*step       

        return
    end subroutine LeapFrogPerturbed3D 


    subroutine RLeapFrogPerturbed3D(r_vector_xyz,v_vector_xyz,step,finalt,k,fr_vector,fv_vector)
        implicit none 
        integer , intent(in) :: finalt
        real(kind(1.q0)), intent(in) :: step , k
        real(kind(1.q0)), dimension(3), intent(inout) :: v_vector_xyz  , r_vector_xyz
        real(kind(1.q0)), dimension(3), intent(out) :: fr_vector ,fv_vector
        real(kind(1.q0)), dimension(3) :: a  
        real(kind(1.q0)) :: r
        real(kind(1.q0)) :: re ,j2,j3
        
        re =6378136.6  
        j2= 1082.64d-6 
        j3 = -2.5d-6


        call dot_product3D(r_vector_xyz,r_vector_xyz,r)

        r = sqrt(r)

        a(1) = -k *( r_vector_xyz(1)/(r**3) ) * ( 1 - j2*(3.0/2.0)*(re/r)**2 * (5*(r_vector_xyz(3)**2/r**2)-1)  )

        a(2) = -k *r_vector_xyz(2)/(r**3)* ( 1 - j2*(3.0/2.0)*(re/r)**2 * (5*(r_vector_xyz(3)**2/r**2)-1)  )

        a(3) = -k *r_vector_xyz(3)/(r**3) * ( 1 + j2*(3.0/2.0)*(re/r)**3 * (3- 5*(r_vector_xyz(3)**2/r**2)) ) 



        fv_vector(1) = v_vector_xyz(1) + a(1)*step
        fv_vector(2) = v_vector_xyz(2) + a(2)*step
        fv_vector(3) = v_vector_xyz(3) + a(3)*step

        fr_vector(1) = r_vector_xyz(1) + fv_vector(1)*step 
        fr_vector(2) = r_vector_xyz(2) + fv_vector(2)*step 
        fr_vector(3) = r_vector_xyz(3) + fv_vector(3)*step       

        return
    end subroutine RLeapFrogPerturbed3D


    
    
    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ---------------------------------------- ORBITAL CALCULATION IN 3D ---------------------------------------
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

    subroutine SemiLatusRectum_a_e(a,e,p)
        implicit none 
        real(kind(1.q0)) , intent(in) :: a ,e 
        real(kind(1.q0)) , intent(out) :: p 

        p = a * (1-e**2)

        return 

    end subroutine SemiLatusRectum_a_e


    

         



    ! ----------------------------------------------------------------------------------------------

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


    subroutine SemiLatusRectum_h_k(h,k,p)
        implicit none 
        real(kind(1.q0)), intent(in) :: h , k 
        real(kind(1.q0)), intent(out) :: p

        p = h**2/k 
        
        return 
    end subroutine SemiLatusRectum_h_k

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


    subroutine MeanMotion(ra,rp,k,m)
        implicit none 
        real(kind(1.q0)) , dimension(3) , intent(in) :: ra , rp
        real(kind(1.q0)) , intent(in) :: k 
        real(kind(1.q0)) , intent(out) :: m ! mean motion 
        real(kind(1.q0)) :: a , aa , pp 

        call dot_product3D(ra,ra, aa)
        call dot_product3D(rp,rp,pp)
        a = SemiMajorAxis(aa,pp)
        m = sqrt(k/a**3)
        return 
    end subroutine MeanMotion

    subroutine SpecificEnergyVR(r,v,k,e)
        implicit none 
        real(kind(1.q0)) , dimension(3) , intent(in) :: r ,v 
        real(kind(1.q0)) , intent(in) :: k 
        real(kind(1.q0)) , intent(out) :: e

        e = v**2/2.0 - k / r 

        return 
    end subroutine SpecificEnergyVR
    

    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ---------------------------------------- ORBITAL ELEMENTS - 2D -------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------

    subroutine KeplersEquationSolver(M,e,Ei1)   ! This should be used for 0<e<1  - M should be in radians.  
        implicit none 
     
        real(kind(1.q0)) :: fe , dfe,ddfe, diff1, diff2,Ei 
        real(kind(1.q0)) , intent(in) :: M,e
        real(kind(1.q0)) , intent(out) ::  Ei1
        integer :: i 

        if (e>0.99 .or. M<0.0699) then 
            Ei = M
            fe = Ei - e*sin(Ei) - M
            dfe = 1.0 - e*cos(Ei)
            ddfe = e*sin(Ei)
                
            Ei1 = Ei -  (5.0*fe) / (  dfe - 2.0*sqrt(abs(4.0*dfe**2 - 5.0*fe*ddfe))   )
            diff1 = abs(Ei1-Ei)
            Ei1 = Ei -  (5.0*fe) / (  dfe + 2.0*sqrt(abs(4.0*dfe**2 - 5.0*fe*ddfe))   )
            diff2 = abs(Ei1-Ei)

            if (diff1>diff2) then 
                Ei = M
                do i = 1,35

                    fe = Ei - e*sin(Ei) - M
                    dfe = 1.0 - e*cos(Ei)
                    ddfe = e*sin(Ei)
                    Ei1 = Ei -  (5.0*fe) / (  dfe + 2.0*sqrt(abs(4.0*dfe**2 - 5.0*fe*ddfe))   )
                    Ei = Ei1 
                end do 
            else
                do i = 1,35
                    
                    fe = Ei - e*sin(Ei) - M
                    dfe = 1.0 - e*cos(Ei)
                    ddfe = e*sin(Ei)
                    Ei1 = Ei -  (5.0*fe) / (  dfe - 2.0*sqrt(abs(4.0*dfe**2 - 5.0*fe*ddfe))   )
                    Ei = Ei1 
                end do 
            end if

        else
            Ei = M
            do i =1,35 
               
                fe = Ei - e*sin(Ei) - M
                dfe = 1.0 - e*cos(Ei)
                Ei1 = Ei - (fe/dfe)
                Ei = Ei1
            end do 

        end if 


        return 

    end subroutine KeplersEquationSolver

    subroutine TrueAnomaly_E(E,ecc,v)
        implicit none 
        real(kind(1.q0)) , intent(in) ::   E,ecc
        real(kind(1.q0)), intent(out) :: v 

        v = acos((cos(E)-ecc)/(1-ecc*cos(E)))

        return
    end subroutine TrueAnomaly_E

 
    real function EccentricAnomaly(ta,e)

        real(kind(1.q0)) , intent(in) :: ta , e 
        EccentricAnomaly = 2.0 * atan(tan(ta/2.0)/sqrt((e+1)/(1-e)))
        
        return

    end function EccentricAnomaly
 
    ! ----------------------------------------------------------------------------------------------

    real function MeanAnomaly(ta,e)

        real(kind(1.q0)) , intent(inout) :: ta , e 
        real(kind(1.q0)) :: ecc_anomaly 

        call init_planets()
        ta = ta * pi/180.0

        ecc_anomaly = 2.0 * atan(tan(ta/2.0)/sqrt((e+1)/(1-e)))

        MeanAnomaly = ecc_anomaly - e*sin(ecc_anomaly)

        return
    end function MeanAnomaly


    ! ----------------------------------------------------------------------------------------------

    real function Eccentricity(ra,rp)
        real(kind(1.q0)) :: ra , rp 

        Eccentricity = (ra-rp)/(ra+rp)
    
        return
    end function Eccentricity

    ! ----------------------------------------------------------------------------------------------

    real function SemiMajorAxis(ra,rp)
        
        real(kind(1.q0)) :: ra , rp 

        SemiMajorAxis = (ra+rp)/2.0
        return
    
    end function SemiMajorAxis

    ! ----------------------------------------------------------------------------------------------
    
    real function SpecificAngularMomentum(k,ra,rp)

        real(kind(1.q0)) :: sm,ecc
        real(kind(1.q0)), intent(in) :: k, ra, rp

        sm = SemiMajorAxis(ra,rp)
        ecc = Eccentricity(ra,rp)


        SpecificAngularMomentum = sqrt(k *sm * (1-ecc**2))
        return
    end function SpecificAngularMomentum

    ! ----------------------------------------------------------------------------------------------

    real function SpecificEnergy(k,ra,rp)

        real(kind(1.q0)) :: sm
        real(kind(1.q0)), intent(in) :: ra , rp ,k 

        sm = SemiMajorAxis(ra,rp)

        SpecificEnergy = -k/(2*sm)
        return
    end function SpecificEnergy

    ! ----------------------------------------------------------------------------------------------

    real function SemiLatusRectum(ra,rp)

        real(kind(1.q0)),intent(in) :: ra,rp 
        real(kind(1.q0)) :: sm , ecc

        sm = SemiMajorAxis(ra,rp)
        ecc = Eccentricity(ra,rp)


        SemiLatusRectum = sm * (1-ecc**2)
        return
    end function SemiLatusRectum

    ! ----------------------------------------------------------------------------------------------------------

    real function Period(k,ra,rp)

        real(kind(1.q0)) , intent(in) :: k ,ra , rp 
        real(kind(1.q0)) :: sm 

        sm = SemiMajorAxis(ra,rp)

        Period = (2.0*pi*sm**(1.5))/sqrt(k)
        return
    end function Period







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





    









end module

    

