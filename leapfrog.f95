program leapfrog
    use astro
    implicit none
    real(kind(1.d0)) :: x , y  , vx0 , vy0 ,step  , gp,fvx,fvy,fx,fy
    integer :: finalt ,wr ,i,fake_step 
    character(1) :: gm


    write(*,*) " You need to enter your x, y , vx0 , vy0 : " 
    read (*,*) x,y,vx0,vy0


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



    open(10, file="leapfrog.dat", action="write" , status = "replace",position="rewind" ,iostat = wr )
    if (wr == 0) write(*,*) "Success!"


    call LeapFrog2D(x,vx0,y,vy0,step,finalt,gp,fx,fy,fvx,fvy)  !First step of the leapfrog. half of step size for the velocity.

    write(10,"(2f30.7)") fx , fy 

    write(*,*) fx , fy , fvx , fvy

    vx0 = fvx 
    vy0 = fvy 
    x = fx 
    y = fy 

    fake_step= finalt/step 

    do i = 1, fake_step 
        call RLeapFrog2D(x,vx0,y,vy0,step,finalt,gp,fx,fy,fvx,fvy)  ! After the first step size, each element will have the same step size. So we use Rleapfrog.

        write(10,"(2f30.7)") fx , fy
        vx0 = fvx 
        vy0 = fvy 
        x = fx 
        y = fy 
    end do 

    close(10)

end program leapfrog




