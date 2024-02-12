program keptoxyz 
    use astro
    implicit none 
    real(kind(1.d0)) :: e , a , i , w , ra , m , x,y,z , jj ,ta


    e = 0.008999999999
    a = 2d7
    i = 0.0
    w = 360.0
    ra = 0.0
    ta = 100.0
    m = MeanAnomaly(ta,e)
    print*, m 

    call KeplarianElementsToXYZ(a,e,m,i,w,ra,ta,x,y,z)

    write(*,"(f30.7)") x 
    write(*,"(f30.7)") y 
    write(*,"(f30.7)") z 

    write(*,"(f30.7)") sqrt(x**2 + y**2 + z**2)

    


end program 

