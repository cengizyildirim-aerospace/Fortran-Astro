program kepler 
    use astro 
    implicit none 

    real(kind(1.q0)) :: m, ecc ,e ,v, p
    
    

    

    write(*,"(a)",advance="no") " Enter your eccentricity value : " 
    read(*,*) e 
    
    write(*,"(a)",advance="no") " Enter your mean anomaly value in rad: " 
    read(*,*) m


    call KeplersEquationSolver(m,e,ecc)
    write(*,*) " Eccentricity Anomaly is .. " , ecc

    call TrueAnomaly_E(ecc,e,v)

    p = pi
    v = v*180/p

    write(*,*) " True Anomaly is .. " , v



end program kepler 

