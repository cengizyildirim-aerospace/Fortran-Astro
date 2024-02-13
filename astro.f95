module astro


    type planet

        character(60) :: name
        real(kind(1.d0)) :: k , radius , mass 

    end type planet

    type (planet) :: Earth , Moon, Sun 



        
    contains

    subroutine init_planets()
        Earth%name = "Earth"
        Earth%k = 3.98600435507d14
        Earth%radius = 6371.0084d3
        Earth%mass = 5.97217d24

        Moon%name = "Moon"
        Moon%k = 4.902800118d12
        Moon%mass = 7.34578924831069d22
        Moon%radius = 1737.4d3

        Sun%name = "Sun"
        Sun%k = 1.32712440041279419d20
        Sun%mass = 1.98840987131d30
        Sun%radius = 695700d3
    end subroutine init_planets





    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! --------------------------------------  LEAP FROG IN 2D  -------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------

   

    subroutine LeapFrog2D(x,vx0,y,vy0,step,finalt,k,fx,fy,fvx,fvy)
        implicit none 
        integer , intent(in) :: finalt
        real(kind(1.d0)) , intent(inout) :: vx0,vy0 ,step , k 
        real(kind(1.d0)) , intent(inout) :: x ,y 
        real(kind(1.d0)) , intent(out) :: fx , fy ,fvx,fvy
        real(kind(1.d0)) :: ax , ay  , r  



        

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
        real(kind(1.d0)) , intent(inout) :: vx0,vy0 ,step , k 
        real(kind(1.d0)) , intent(inout) :: x ,y 
        real(kind(1.d0)) , intent(out) :: fx , fy ,fvx,fvy
        real(kind(1.d0)) :: ax , ay  , r  



        

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
    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------





    
    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ---------------------------------------- ORBITAL CALCULATION IN 3D ---------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------






    subroutine KeplarianElementsToXYZ(a,e,m,i,w,ra,ta,x,y,z)   ! https://ntrs.nasa.gov/api/citations/19650015945/downloads/19650015945.pdf
        implicit none
        real(kind(1.d0)), intent(inout) :: a,e,m,i,w,ra ,ta
        real(kind(1.d0)), intent(out) :: x , y ,z 
        real(kind(1.d0)) :: ecc_anomaly 
        integer :: b 
        
        i = i*3.14159265358979323/180.0
        w = w*3.14159265358979323/180.0
        ra = ra*3.14159265358979323/180.0

        ecc_anomaly = EccentricAnomaly(ta,e)
        

        write(*,"(a , f30.7)") "Eccentric anomoly is : " , ecc_anomaly

        x = a * ((cos(ecc_anomaly)-e)*(cos(w)*cos(ra)-sin(w)*sin(ra)*cos(i))+ &
        sqrt(1-e**2) * sin(ecc_anomaly)* (-sin(w)*cos(ra)-cos(w)*sin(ra)*cos(i)))

        y = a * ((cos(ecc_anomaly)-e)*(cos(w)*sin(ra)+sin(w)*cos(ra)*cos(i))+  &
        sqrt(1-e**2) * sin(ecc_anomaly)* (-sin(w)*sin(ra)+cos(w)*cos(ra)*cos(i)))

        z = a * ((cos(ecc_anomaly)-e)*(sin(w)*sin(i))+  sqrt(1-e**2) * sin(ecc_anomaly)* (cos(w)*sin(i)))

    
        return

    end subroutine KeplarianElementsToXYZ

    ! ----------------------------------------------------------------------------------------------

    subroutine EccentricityVector(r_vector_xyz,v_vector_xyz,k,ecc_vector)
        real(kind(1.d0)) , dimension(3) ,intent(inout) :: r_vector_xyz , v_vector_xyz
        real(kind(1.d0)) , dimension(3) ,intent(out) :: ecc_vector
        real(kind(1.d0)) , dimension(3) :: h , j  ! specific angular momentum and  j= vxh  BATE, page 25 
        real(kind(1.d0)) , intent(in) :: k  ! gravitaonal parameter (GM)
        real(kind(1.d0)) :: v_total , r_total , v, r, s    !s = r * v

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
        real(kind(1.d0)) , dimension(3) ,intent(in) :: h
        real(kind(1.d0)) , dimension(3) ,intent(out) :: node_vector
        real(kind(1.d0)) , dimension(3)  :: k 

        k(1) = 0 
        k(2) = 0
        k(3) = 1

        call cross_product3D(k,h,node_vector)

        return
    end subroutine NodeVector

    subroutine SpecificAngularMomentumVector(r_vector_xyz,v_vector_xyz,h)
        implicit none
        real(kind(1.d0)) , dimension(3) ,intent(in) :: r_vector_xyz , v_vector_xyz
        real(kind(1.d0)) , dimension(3) ,intent(out) :: h 
    
        call cross_product3D(r_vector_xyz,v_vector_xyz,h)

        return
    
    end subroutine SpecificAngularMomentumVector

    ! -------------------------------------------------------------------------------------------------
    ! -----------------------------ORBITAL ELEMENTS (KEPLARIAN) ---------------------------------------
    ! -------------------------------------------------------------------------------------------------


    subroutine inclination(h,i)
        implicit none 
        real(kind(1.d0)) , dimension(3), intent(in) :: h 
        real(kind(1.d0)) , intent(out) :: i 
        real(kind(1.d0)) :: h_total

        call dot_product3D(h,h,h_total)
        h_total = sqrt(h_total)
        i = acos(h(3)/ h_total)
        i = i*180.0/3.14159265358979323
        
        return


    end subroutine inclination


    subroutine LongitudeOfTheAscendingNode(node_vector,ascending_node)
        implicit none 
        real(kind(1.d0)) , dimension(3), intent(in) :: node_vector
        real(kind(1.d0)) , intent(out) :: ascending_node 
        real(kind(1.d0)) :: n_total

        call dot_product3D(node_vector,node_vector,n_total)
        n_total = sqrt(n_total)
        ascending_node = acos(node_vector(1)/n_total)
        ascending_node = ascending_node*180/3.14159265358979323

        if (node_vector(2)<0) ascending_node = 360.0 - ascending_node

        
        return
    end subroutine LongitudeOfTheAscendingNode


    subroutine ArgumentOfPeriapsis(node_vector,EccentricityVector,w)
        real(kind(1.d0)) , dimension(3), intent(in) :: node_vector,EccentricityVector 
        real(kind(1.d0)) , intent(out) :: w
        real(kind(1.d0)) :: n_total , e_total , dot_ne

        call dot_product3D(node_vector,node_vector,n_total)
        call dot_product3D(EccentricityVector,EccentricityVector,e_total)
        call dot_product3D(node_vector,EccentricityVector,dot_ne)


        n_total = sqrt(n_total)
        e_total = sqrt(e_total)

        w = acos(dot_ne/(n_total*e_total))
        w = w*180.0/3.14159265358979323

        if (EccentricityVector(3)<0) w = 360.0 - w 

        return
    
    end subroutine ArgumentOfPeriapsis


    subroutine TrueAnomaly(EccentricityVector,r_vector_xyz,v_vector_xyz,ta)
        implicit none 
        real(kind(1.d0)) , dimension(3), intent(in) :: r_vector_xyz,EccentricityVector ,v_vector_xyz
        real(kind(1.d0)) , intent(out) :: ta
        real(kind(1.d0)) :: r_total , e_total , dot_re, checker

        call dot_product3D(r_vector_xyz, r_vector_xyz ,r_total)
        call dot_product3D(EccentricityVector,EccentricityVector,e_total)

        call dot_product3D(EccentricityVector,r_vector_xyz,dot_re)
        call dot_product3D(r_vector_xyz,v_vector_xyz,checker)

        e_total =sqrt(e_total)
        r_total =sqrt(r_total)

        ta = acos(dot_re/(e_total*r_total))

        ta = ta*180.0/3.14159265358979323

        if (checker<0) ta =360.0 - ta 


        return
    
    end subroutine TrueAnomaly


    subroutine SemiLatusRectum_h_k(h,k,p)
        implicit none 
        real(kind(1.d0)), intent(in) :: h , k 
        real(kind(1.d0)), intent(out) :: p

        p = h**2/k 
        
        return 
    end subroutine SemiLatusRectum_h_k

    subroutine SemiMajorAxis_p_e(p,e,a)
        implicit none 
        real(kind(1.d0)), intent(in) :: p , e 
        real(kind(1.d0)), intent(out) :: a

        a = p / (1-e**2)

        return 

    end subroutine SemiMajorAxis_p_e





    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ---------------------------------------- ORBITAL ELEMENTS - 2D -------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------



 
    real function EccentricAnomaly(ta,e)

        real(kind(1.d0)) , intent(in) :: ta , e 
        EccentricAnomaly = 2.0 * atan(tan(ta/2.0)/sqrt((e+1)/(1-e)))
        
        return

    end function EccentricAnomaly
 
    ! ----------------------------------------------------------------------------------------------

    real function MeanAnomaly(ta,e)

        real(kind(1.d0)) , intent(inout) :: ta , e 
        real(kind(1.d0)) :: ecc_anomaly 

        ta = ta * 3.14159265358979323/180.0

        ecc_anomaly = 2.0 * atan(tan(ta/2.0)/sqrt((e+1)/(1-e)))

        MeanAnomaly = ecc_anomaly - e*sin(ecc_anomaly)

        return
    end function MeanAnomaly


    ! ----------------------------------------------------------------------------------------------

    real function Eccentricity(ra,rp)
        integer, parameter :: dp = kind(1.d0)
        real(dp) :: ra , rp 

        Eccentricity = (ra-rp)/(ra+rp)
    
        return
    end function Eccentricity

    ! ----------------------------------------------------------------------------------------------

    real function SemiMajorAxis(ra,rp)
        integer, parameter :: dp = kind(1.d0)
        real(dp) :: ra , rp 

        SemiMajorAxis = (ra+rp)/2.0
        return
    
    end function SemiMajorAxis

    ! ----------------------------------------------------------------------------------------------
    
    real function SpecificAngularMomentum(k,ra,rp)

        real(kind(1.d0)) :: sm,ecc
        real(kind(1.d0)), intent(in) :: k, ra, rp

        sm = SemiMajorAxis(ra,rp)
        ecc = Eccentricity(ra,rp)


        SpecificAngularMomentum = sqrt(k *sm * (1-ecc**2))
        return
    end function SpecificAngularMomentum

    ! ----------------------------------------------------------------------------------------------

    real function SpecificEnergy(k,ra,rp)

        real(kind(1.d0)) :: sm
        real(kind(1.d0)), intent(in) :: ra , rp ,k 

        sm = SemiMajorAxis(ra,rp)

        SpecificEnergy = -k/(2*sm)
        return
    end function SpecificEnergy

    ! ----------------------------------------------------------------------------------------------

    real function SemiLatusRectum(ra,rp)

        real(kind(1.d0)),intent(in) :: ra,rp 
        real(kind(1.d0)) :: sm , ecc

        sm = SemiMajorAxis(ra,rp)
        ecc = Eccentricity(ra,rp)


        SemiLatusRectum = sm * (1-ecc**2)
        return
    end function SemiLatusRectum

    ! ----------------------------------------------------------------------------------------------------------

    real function Period(k,ra,rp)

        real(kind(1.d0)) , intent(in) :: k ,ra , rp 
        real(kind(1.d0)) :: sm

        sm = SemiMajorAxis(ra,rp)

        Period = (2.0*3.14159265358979323*sm**(1.5))/sqrt(k)
        return
    end function Period







    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ---------------------------------------- LINEAR OPERATIONS  ----------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------------------------------





    subroutine cross_product3D(a,b,c)
        
        implicit none

        real(kind(1.d0)), dimension(3), intent(in) :: a,b  ! the vectors of a x b .. they should be placed in this order!
        real(kind(1.d0)), dimension(3), intent(out) :: c   ! the resulting vector 
        

        ! a(1) is x   , a(2) is y and a(3) is z  in xyz cartesian coordinate system.

        c(1) = a(2)*b(3) - b(2)*a(3)
        c(2) = -(a(1)*b(3) - b(1)*a(3))
        c(3) = a(1)*b(2) - b(1)*a(2)

        return 

    end subroutine cross_product3D


    subroutine dot_product3D(a,b,c)

        implicit none 

        real(kind(1.d0)), dimension(3), intent(in) :: a,b  ! the vectors of a * b .. they should be placed in this order!
        real(kind(1.d0)), intent(out) :: c   !the result 

        c = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

        return 
        
    end subroutine dot_product3D



    subroutine angle_between_vectors(v1,v2,angle)  ! It will print it in rad 
        implicit none
        real(kind(1.d0)) , dimension(3) , intent(in) :: v1 ,v2 
        real(kind(1.d0)) , intent(out) :: angle
        real(kind(1.d0)) :: dt_product , v1_m , v2_m

        call dot_product3D(v1,v2,dt_product)
        call dot_product3D(v1,v1,v1_m)
        call dot_product3D(v2,v2,v2_m)

        angle = acos(dt_product / (  sqrt(v1_m)*sqrt(v2_m)  )) 

        return
    end subroutine angle_between_vectors



end module

    

