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


    subroutine KeplarianElementsToXYZ(a,e,m,i,w,ra,ta,x,y,z)   ! https://ntrs.nasa.gov/api/citations/19650015945/downloads/19650015945.pdf
        implicit none
        real(kind(1.d0)), intent(inout) :: a,e,m,i,w,ra ,ta
        real(kind(1.d0)), intent(out) :: x , y ,z 
        real(kind(1.d0)) :: ecc_anomaly  , dE , ecc
        integer :: b , n = 29 
        
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


    real function EccentricAnomaly(ta,e)
        real(kind(1.d0)) , intent(in) :: ta , e 
        

        EccentricAnomaly = 2.0 * atan(tan(ta/2.0)/sqrt((e+1)/(1-e)))
        
        return
    end function EccentricAnomaly
 

    real function MeanAnomaly(ta,e)
        real(kind(1.d0)) , intent(inout) :: ta , e 
        real(kind(1.d0)) :: ecc_anomaly 

        ta = ta * 3.14159265358979323/180.0

        ecc_anomaly = 2.0 * atan(tan(ta/2.0)/sqrt((e+1)/(1-e)))

        MeanAnomaly = ecc_anomaly - e*sin(ecc_anomaly)

        return
    end function MeanAnomaly

        









          













    real function Eccentricity(ra,rp)
        integer, parameter :: dp = kind(1.d0)
        real(dp) :: ra , rp 

        Eccentricity = (ra-rp)/(ra+rp)
    
        return
    end function Eccentricity

    real function SemiMajorAxis(ra,rp)
        integer, parameter :: dp = kind(1.d0)
        real(dp) :: ra , rp 

        SemiMajorAxis = (ra+rp)/2.0
        return
    
    end function SemiMajorAxis
    
    real function SpecificAngularMomentum(k,ra,rp)
        integer, parameter :: dp = kind(1.d0)
        real(dp) :: sm,ecc
        real(dp), intent(in) :: k, ra, rp

        sm = SemiMajorAxis(ra,rp)
        ecc = Eccentricity(ra,rp)


        SpecificAngularMomentum = sqrt(k *sm * (1-ecc**2))
        return
    end function SpecificAngularMomentum

    real function SpecificEnergy(k,ra,rp)
        integer, parameter :: dp = kind(1.d0)
        real(dp) :: sm
        real(dp), intent(in) :: ra , rp ,k 

        sm = SemiMajorAxis(ra,rp)

        SpecificEnergy = -k/(2*sm)
        return
    end function SpecificEnergy

    real function SemiLatusRectum(ra,rp)
        integer, parameter :: dp = kind(1.d0)
        real(dp),intent(in) :: ra,rp 
        real(dp) :: sm , ecc

        sm = SemiMajorAxis(ra,rp)
        ecc = Eccentricity(ra,rp)


        SemiLatusRectum = sm * (1-ecc**2)
        return
    end function SemiLatusRectum

    real function Period(k,ra,rp)
        integer, parameter :: dp = kind(1.d0)
        real(dp) , intent(in) :: k ,ra , rp 
        real(dp) :: sm , pi=3.14159265358979323_dp

        sm = SemiMajorAxis(ra,rp)

        Period = (2.0*pi*sm**(1.5))/sqrt(k)
        return
    end function Period



        




end module

    

