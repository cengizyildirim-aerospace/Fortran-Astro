program leap3d

    use astro 
    implicit none 
    character(1) :: gm 
    real(kind(1.d0)), dimension(3) :: r_vector_xyz , v_vector_xyz ,fr_vector,fv_vector
    real(kind(1.d0)) :: gp , step
    integer :: fake_step,wr,i,finalt



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
    else if (gm =="s" ) then
        call init_planets()
        gp = Sun%k
    else if (gm =="m" ) then 
        call init_planets()
        gp = Moon%k 
    else 
        write(*,*) "You need to enter a character and it can only be s for the Sun - e for the Earth and m for the Moon"
    end if 

    open(10, file="leapfrog3d.plt", action="write" , status = "replace",position="rewind" ,iostat = wr )
    if (wr == 0) write(*,*) "Success!"

    write (10, '(a)') "# Orbital Plot Example " 
    write (10, '(a)') "set title ""Molniya Orbit""" 
    write (10,'(a,a)') "splot ""-"" notitle with dots ", "linewidth 2 linecolor 2"


    call LeapFrog3D(r_vector_xyz,v_vector_xyz,step,finalt,gp,fr_vector,fv_vector)

    write(10,"(3f17.7)") fr_vector(1), fr_vector(2) , fr_vector(3)

    r_vector_xyz(1) = fr_vector(1)
    r_vector_xyz(2) = fr_vector(2)
    r_vector_xyz(3) = fr_vector(3)


    v_vector_xyz(1) = fv_vector(1)      
    v_vector_xyz(2) = fv_vector(2)
    v_vector_xyz(3) = fv_vector(3)



    fake_step= finalt/step 

    do i = 1,fake_step
        call RLeapFrog3D(r_vector_xyz,v_vector_xyz,step,finalt,gp,fr_vector,fv_vector)  ! After the first step size, each element will have the same step size. So we use Rleapfrog.

        write(10,"(3f17.7)") fr_vector(1), fr_vector(2) , fr_vector(3)

        r_vector_xyz(1) = fr_vector(1)
        r_vector_xyz(2) = fr_vector(2)
        r_vector_xyz(3) = fr_vector(3)


        v_vector_xyz(1) = fv_vector(1)      
        v_vector_xyz(2) = fv_vector(2)
        v_vector_xyz(3) = fv_vector(3) 

    end do 

    close(10)

end program leap3d