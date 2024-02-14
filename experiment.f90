program orbit
    use astro
    implicit none 
    
    real(kind(1.d0)) :: ra,rp,spcangmon , kSun , va0,vp0,ecc, T, kEarth, sm ,spe
    
    call init_planets()

    print *, "Earth radius is : " , Earth%radius

    print *, "The gravitational parameter of the Sun : " , Sun%k 

    kSun = 1.32712440041279419d20

    kEarth = 3.98600435507d14



    write(*,"(a)",advance="no") "Enter your ra and rp: "
    read(*,*) ra , rp 

    ecc = Eccentricity(ra,rp)
    write(*,"(a,f30.9,/)") "Here is your Eccentricity:" , ecc

    sm = SemiMajorAxis(ra,rp)
    write(*,"(a,f30.9,/)") "Here is your SemiMajorAxis:" ,sm


    T = Period(kEarth,ra,rp)
    write(*,"(a,f30.9,/)")  "Here is your period in seconds!: " , T
    

    spcangmon = SpecificAngularMomentum(kEarth,ra,rp)
    write(*,"(a,f30.9,/)") "Here is the  Specific Angular momentum:" ,spcangmon

    va0 = spcangmon/ra
    write(*,"(a,f30.9,/)") "Here is the velocity at apogee" ,  va0

    vp0 = spcangmon/rp
    write(*,"(a,f30.9,/)") "Here is the velocity at perigee" ,  vp0

    spe = SpecificEnergy(kEarth,ra,rp)
    write(*,"(a,f30.9,/)") "Here is the specific energy:" ,  spe

    




end program orbit

