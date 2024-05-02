program solar 
    use astro 

    implicit none 
    real(kind(1.q0)) :: rrr , step , vv , rr_orbit , rr_moon , p  , m , a , energy0, totale, test,test1
    real(kind(1.q0)), dimension(3)  :: relative_r  , rr , relative_v , vf , r_moon,v_moon, rf
    
    type planets
        real(kind(1.q0)),dimension(3) ::  p , velocity , acc 
        real(kind(1.q0)) :: gm ,mass
        real(kind(1.q0)) :: energy , v  ! v stand for the total velocity of the object!
        character(15) :: name 
    end type planets

    type (planets),dimension(25) :: planett

    integer(kind=16) :: i , j , t , k , wr , fake_step
    integer(kind=16) :: finalt,day

    ! Lets start by giving the iniital conditions of the planets. Our "almost" intertial reference frame is the Solar System Berycenter
    ! This means the center of the Sun is NOT our center.


    call init_planets()

    
    ! The position and the velocity of every planet on Feb 8, 2024.


    planett(1)%name = "The Sun"
    planett(1)%p(1) = -1.163537052918776q9     ! The sun 
    planett(1)%p(2) = -4.740065576137778q8
    planett(1)%p(3) = 3.118452641906473q7
 
    planett(1)%velocity(1) = 8.893450715842641q0 
    planett(1)%velocity(2) = -1.170047290876101q1
    planett(1)%velocity(3) = -9.076897335120841q-2


    planett(2)%name = "Mercury"
    planett(2)%p(1) = 1.175745966805013q9   ! Mercury
    planett(2)%p(2) = -6.942068692657061q10
    planett(2)%p(3) = -5.817764387641165q9

    planett(2)%velocity(1) = 3.892344716347844q4
    planett(2)%velocity(2) =  4.139599358081010q3
    planett(2)%velocity(3) = -3.23022493584299q3


    planett(3)%name = "Venus"
    planett(3)%p(1) = -4.346854652892029q10 !Venus 
    planett(3)%p(2) =  -1.003969275596239q11
    planett(3)%p(3) = 1.099930209689982q9

    planett(3)%velocity(1) = 3.202060594616277q4  
    planett(3)%velocity(2) =  -1.381876570106473q4
    planett(3)%velocity(3) = -2.036812982080393q3


    planett(4)%name = "The Earth"
    planett(4)%p(1) =  -1.114482120509732q11     !the Earth
    planett(4)%p(2) =  9.755091445282833q10 
    planett(4)%p(3) =  2.677803510279208q7

    planett(4)%velocity(1) = -2.027697054595732q4       
    planett(4)%velocity(2) = -2.239420902641540q4 
    planett(4)%velocity(3) =  1.132016626611332q0


    planett(5)%name = "Mars System"
    planett(5)%p(1) =  3.701864765218908q10  ! Mars system barycenter
    planett(5)%p(2) =  -2.117033921220820q11 
    planett(5)%p(3) =  -5.332148864767030q9

    planett(5)%velocity(1) = 2.476487728790323q4        
    planett(5)%velocity(2) = 6.378058254790010q3 
    planett(5)%velocity(3) = -4.734151917922453q2



    planett(6)%name = "Jupiter System"
    planett(6)%p(1) = 4.894014555935882q11  ! Jupyter Barycenter
    planett(6)%p(2) = 5.625334140898753q11
    planett(6)%p(3) = -1.328302395332387q10

    planett(6)%velocity(1) = -1.000226871469186q4     
    planett(6)%velocity(2) =  9.195529107270879q3
    planett(6)%velocity(3) =  1.856381984935207q2


    planett(7)%name = "Saturn System"
    planett(7)%p(1) = 1.354643742479730q12   ! Saturn Barycenter
    planett(7)%p(2) = -5.269947919948134q11
    planett(7)%p(3) = -4.477177953209031q10  


    planett(7)%velocity(1) = 2.96308837865079q3   
    planett(7)%velocity(2) = 8.983275944771959q3
    planett(7)%velocity(3) = -2.741792164834611q2

    
    planett(8)%name = "Uranus System"
    planett(8)%p(1) = 1.816861039178873q12   ! Uranus barycenter 
    planett(8)%p(2) = 2.301334475229118q12
    planett(8)%p(3) = -1.499057790727472q10

    planett(8)%velocity(1) = -5.395049652933221q3       
    planett(8)%velocity(2) =  3.902465865270330q3
    planett(8)%velocity(3) =  8.43856764548958q1


    planett(9)%name = "Neptune System"
    planett(9)%p(1) =  4.464171781366321q12   ! Neptune barycenter
    planett(9)%p(2) = -2.504295457593737q11
    planett(9)%p(3) =  -9.772425156508410q10

    planett(9)%velocity(1) = 2.685040519735848q2      
    planett(9)%velocity(2) =  5.458978355655157q3
    planett(9)%velocity(3) = -1.186052963155491q2


    planett(10)%name = "The Moon"
    planett(10)%p(1) = -1.113175726324003q11   ! Moon
    planett(10)%p(2) = 9.721142863263611q10
    planett(10)%p(3) =  -5.449223300203681q6

    planett(10)%velocity(1) = -1.929177591167254q4      
    planett(10)%velocity(2) =  -2.195669409085904q4
    planett(10)%velocity(3) = 1.249049169100225q1



    planett(11)%name = "Pluto System"
    planett(11)%p(1) =   2.589355494906734q12   ! Pluto barycenter
    planett(11)%p(2) = -4.534120993663620q12
    planett(11)%p(3) = -2.638172950988719q11

    planett(11)%velocity(1) = 4.861058210676243q3      
    planett(11)%velocity(2) = 1.502368640131439q3
    planett(11)%velocity(3) = -1.566868320916638q3

    







    
    ! -------------------------------------------------

    planett(1)%mass = Sun%mass
    planett(2)%mass = Mercury%mass
    planett(3)%mass = Venus%mass
    planett(4)%mass = Earth%mass
    planett(5)%mass = Marss%mass
    planett(6)%mass = Jupiters%mass
    planett(7)%mass = Saturns%mass
    planett(8)%mass = Uranuss%mass
    planett(9)%mass = Neptunes%mass
    planett(10)%mass = Moon%mass
    planett(11)%mass = Plutos%mass
    
    

    !----------------------------------------------------

    planett(1)%gm = Sun%k
    planett(2)%gm = Mercury%k 
    planett(3)%gm = Venus%k
    planett(4)%gm = Earth%k 
    planett(5)%gm = Marss%k 
    planett(6)%gm = Jupiter%k
    planett(7)%gm = Saturns%k 
    planett(8)%gm = Uranuss%k
    planett(9)%gm = Neptunes%k
    planett(10)%gm = Moon%k
    planett(11)%gm = Plutos%k
  
    

    
    



   
    


    open(10, file="SunPosition.dat", action="write", position="rewind", status="replace", iostat=wr)
    if (wr == 0) write(*, *) "Success!"

    open(11, file="MercuryPosition.dat", action="write", position="rewind", status="replace", iostat=wr)
    if (wr == 0) write(*, *) "Success!"

    open(12, file="VenusPosition.dat", action="write", position="rewind", status="replace", iostat=wr)
    if (wr == 0) write(*, *) "Success!"

    open(13, file="EarthPosition.dat", action="write", position="rewind", status="replace", iostat=wr)
    if (wr == 0) write(*, *) "Success!"

    open(14, file="MarsPosition.dat", action="write", position="rewind", status="replace", iostat=wr)
    if (wr == 0) write(*, *) "Success!"

    open(15, file="JupiterPosition.dat", action="write", position="rewind", status="replace", iostat=wr)
    if (wr == 0) write(*, *) "Success!"

    open(16, file="SaturnPosition.dat", action="write", position="rewind", status="replace", iostat=wr)
    if (wr == 0) write(*, *) "Success!"

    open(17, file="UranusPosition.dat", action="write", position="rewind", status="replace", iostat=wr)
    if (wr == 0) write(*, *) "Success!"

    open(18, file="NeptunePosition.dat", action="write", position="rewind", status="replace", iostat=wr)
    if (wr == 0) write(*, *) "Success!"

    open(19, file="MoonPosition.dat", action="write", position="rewind", status="replace", iostat=wr)
    if (wr == 0) write(*, *) "Success!"

    open(20, file="PlutoPosition.dat", action="write", position="rewind", status="replace", iostat=wr)
    if (wr == 0) write(*, *) "Success!"

    open(22, file="WholeSystemEnergyChange.dat", action="write", position="rewind", status="replace", iostat=wr)
    if (wr == 0) write(*, *) "Success!"

    open(23, file="WholeSystemEnergy.dat", action="write", position="rewind", status="replace", iostat=wr)
    if (wr == 0) write(*, *) "Success!"




   


    !------------------------------------------------- 

    write(*,*) " For how long should we predict the position of the planets? ( seconds-integer!) : " 
    write(*,*) " It should be a multiple of 86400 seconds (1 day)! " // & 
    " Our code writes down the coordinates every 24 hours to save some time. "
    read(*,*) t 

    write(*,*) " Tell us the step size : " 
    read(*,*) step 

    t = int(t/step)

    

    energy0 = 0 
    test = 0 

    do k = 1, t 

        day = int(k*step)
        

        if (day == 0) then 
            day = 1 
        end if 
        

        if (mod(day,86400)== 0 .and. day /= 0) then 
            write(*,*) "day : " , int(day/86400.0) 
        end if 

        do i = 1, 11
            planett(i)%acc = 0 
        end do 


        do i = 1,11
            do j = 1 ,11

                if (i == j )  then 
                    
                    continue 
                
                
                else 
                
                    relative_r = planett(i)%p - planett(j)%p 
                    
                    

                    call dot_product3D(relative_r,relative_r,rrr) 

                    rrr = sqrt(abs(rrr))

                    
                     
                    
                    planett(i)%acc(1) = planett(i)%acc(1) - planett(j)%gm * relative_r(1) / rrr**3
                    planett(i)%acc(2) = planett(i)%acc(2) - planett(j)%gm * relative_r(2) / rrr**3
                    planett(i)%acc(3) = planett(i)%acc(3) - planett(j)%gm * relative_r(3) / rrr**3                   

                    
                end if 

            end do 
        end do 

        

        if (k == 1) then 

            do i = 1,11

                planett(i)%velocity(1) = planett(i)%velocity(1) + planett(i)%acc(1)*step/real(2)
                planett(i)%velocity(2) = planett(i)%velocity(2) + planett(i)%acc(2)*step/real(2)
                planett(i)%velocity(3) = planett(i)%velocity(3) + planett(i)%acc(3)*step/real(2)

                planett(i)%p(1) = planett(i)%p(1) + planett(i)%velocity(1)*step 
                planett(i)%p(2) = planett(i)%p(2) + planett(i)%velocity(2)*step 
                planett(i)%p(3) = planett(i)%p(3) + planett(i)%velocity(3)*step 

            end do

            
        else 
            do i = 1,11
                planett(i)%velocity(1) = planett(i)%velocity(1) + planett(i)%acc(1)*step
                planett(i)%velocity(2) = planett(i)%velocity(2) + planett(i)%acc(2)*step
                planett(i)%velocity(3) = planett(i)%velocity(3) + planett(i)%acc(3)*step

                planett(i)%p(1) = planett(i)%p(1) + planett(i)%velocity(1)*step 
                planett(i)%p(2) = planett(i)%p(2) + planett(i)%velocity(2)*step 
                planett(i)%p(3) = planett(i)%p(3) + planett(i)%velocity(3)*step

                if (mod(day,86400) == 0 .or. k == t) then 
            
                    if (i == 1) then 
                        write(10,"(f30.6,x,f30.6,x,f30.6)") planett(i)%p(1),planett(i)%p(2),planett(i)%p(3)
                    else if ( i == 2) then 
                        write(11,"(f30.6,x,f30.6,x,f30.6)") planett(i)%p(1),planett(i)%p(2),planett(i)%p(3)
                    else if ( i == 3) then 
                        write(12,"(f30.6,x,f30.6,x,f30.6)") planett(i)%p(1),planett(i)%p(2),planett(i)%p(3)
                    else if ( i == 4) then 
                        write(13,"(f30.6,x,f30.6,x,f30.6)") planett(i)%p(1),planett(i)%p(2),planett(i)%p(3)
                    else if ( i == 5) then 
                        write(14,"(f30.6,x,f30.6,x,f30.6)") planett(i)%p(1),planett(i)%p(2),planett(i)%p(3)
                    else if ( i == 6) then 
                        write(15,"(f30.6,x,f30.6,x,f30.6)") planett(i)%p(1),planett(i)%p(2),planett(i)%p(3)
                    else if ( i == 7) then 
                        write(16,"(f30.6,x,f30.6,x,f30.6)") planett(i)%p(1),planett(i)%p(2),planett(i)%p(3)
                    else if ( i == 8) then 
                        write(17,"(f30.6,x,f30.6,x,f30.6)") planett(i)%p(1),planett(i)%p(2),planett(i)%p(3)
                    else if ( i == 9) then 
                        write(18,"(f30.6,x,f30.6,x,f30.6)") planett(i)%p(1),planett(i)%p(2),planett(i)%p(3)
                    else if ( i == 10) then 
                        write(19,"(f30.6,x,f30.6,x,f30.6)") planett(i)%p(1),planett(i)%p(2),planett(i)%p(3)
                    else if ( i == 11) then 
                        write(20,"(f30.6,x,f30.6,x,f30.6)") planett(i)%p(1),planett(i)%p(2),planett(i)%p(3)
                    end if 

                end if 

                


            end do 
        end if 

        

        if (k>1 .and. mod(k,360)== 0) then
            do i = 1, 11 
                planett(i)%energy = 0 
            end do 

            totale = 0 

            do i = 1,11
                do j = 1,11
                    
                    if (i ==j )  then 
                            
                        continue 
                        
                        
                    else 
                        relative_r = planett(i)%p - planett(j)%p 
                        
                        
   
                        call dot_product3D(relative_r,relative_r,rrr) 

                        rrr = sqrt(abs(rrr))
                    
                        planett(i)%energy = planett(i)%energy  - planett(i)%mass*planett(j)%gm /(2.0*rrr )

                        
                        
                            

                            
                    end if 

                end do

                call dot_product3D(planett(i)%velocity,planett(i)%velocity,vv) 
                    
                planett(i)%energy = planett(i)%energy + planett(i)%mass*(vv) /2.0 

                

            end do  


            do i = 1, 11
                totale = totale + planett(i)%energy  
            end do 

            if (abs((totale-energy0)/energy0)>test ) then 
                test = abs((totale-energy0)/energy0)
                write(22,"(i14,x,f60.20)") day, test

            end if 

            test1 = abs((totale-energy0)/energy0) 
            write(23,"(i14,x,f60.20)") day, test1

        else if (k==1) then
            do i = 1, 11 
                planett(i)%energy = 0 
            end do 

            do i = 1,11
                do j = 1,11
                    
                    if (j == i )  then 
                            
                        continue 
                        
                        
                    else 
                        relative_r = planett(i)%p - planett(j)%p 
                        
                        

                        call dot_product3D(relative_r,relative_r,rrr) 

                        rrr = sqrt(abs(rrr))
                    
                        planett(i)%energy = planett(i)%energy  - planett(i)%mass*planett(j)%gm /(2.0*rrr) 

                        
                        
                            

                            
                    end if 

                end do

                call dot_product3D(planett(i)%velocity,planett(i)%velocity,vv) 
                    
                planett(i)%energy = planett(i)%energy + planett(i)%mass*(vv) /2.0 

                

            end do  


            do i = 1,11
                energy0 = energy0 + planett(i)%energy  
            end do  

            

            write(23,"(i14,x,f60.20)") day, real(0)

            
        
            

            

        end if  
            
       
        


        

    end do 



    
    write(*,"(a,f50.6)") " The position of the Earth is x: " , planett(4)%p(1) 
    write(*,"(a,f50.6)") " The position of the Earth is y: " , planett(4)%p(2)
    write(*,"(a,f50.6)") " The position of the Earth is z: " , planett(4)%p(3)  

    write(*,"(a,f50.6)") " Vx Earth: " , planett(4)%velocity(1) 
    write(*,"(a,f50.6)") " Vy Earth: " , planett(4)%velocity(2)
    write(*,"(a,f50.6)") " Vz Earth : " , planett(4)%velocity(3)


    do i = 1, 11 
        call dot_product3D(planett(i)%velocity,planett(i)%velocity,planett(i)%v) 
        planett(i)%v = sqrt(planett(i)%v)
        write(*,"(a,a,f50.6)") " V of " , planett(i)%name, planett(i)%v
    end do 
    
     


    
    

    

    

    

    









    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)
    close(19)
    close(20)
    close(22)
    close(23)





    

   


end program solar 

