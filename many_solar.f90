program many_solar 
    use astro 

    implicit none 
    real(kind(1.q0)) :: rrr , step , vv , rr_orbit , rr_moon , p  , m , a , energy0, totale, test,test1
    real(kind(1.q0)), dimension(3)  :: relative_r  , rr , relative_v , vf , r_moon,v_moon, rf
    
    type planets
        real(kind(1.q0)),dimension(3) ::  p , velocity , acc 
        real(kind(1.q0)) :: gm ,mass
        real(kind(1.q0)) :: energy 
    end type planets

    type (planets),dimension(25) :: planett

    integer(kind=16) :: i , j , t , k , wr , fake_step
    integer(kind=16) :: finalt,day

    ! Lets start by giving the iniital conditions of the planets. Our "almost" intertial reference frame is the Solar System Berycenter
    ! This means the center of the Sun is NOT our center.


    call init_planets()

    
    ! The position and the velocity of every planet on Feb 8, 2024.

    planett(1)%p(1) = -1163537052.918776q0     ! The sun 
    planett(1)%p(2) = -474006557.6137778q0
    planett(1)%p(3) = 31184526.41906473q0
 
    planett(1)%velocity(1) = 8.893450715842641q0 
    planett(1)%velocity(2) = -11.70047290876101q0
    planett(1)%velocity(3) = -0.09076897335120841q0



    planett(2)%p(1) = 1175745966.805013q0   ! Mercury
    planett(2)%p(2) = -69420686926.57061q0
    planett(2)%p(3) = -5817764387.641165q0

    planett(2)%velocity(1) = 38923.44716347844q0 
    planett(2)%velocity(2) =  4139.599358081010q0
    planett(2)%velocity(3) = -3230.224935842994q0



    planett(3)%p(1) = -43468546528.92029q0 !Venus 
    planett(3)%p(2) =  -100396927559.6239q0
    planett(3)%p(3) =  1099930209.689982q0

    planett(3)%velocity(1) = 32020.60594616277q0  
    planett(3)%velocity(2) =  -13818.76570106473q0
    planett(3)%velocity(3) = -2036.812982080393q0



    planett(4)%p(1) =  -111448212050.9732q0     !the Earth
    planett(4)%p(2) =  97550914452.82833q0 
    planett(4)%p(3) =  26778035.10279208q0

    planett(4)%velocity(1) = -20276.97054595732q0       
    planett(4)%velocity(2) = -22394.20902641540q0 
    planett(4)%velocity(3) =  1.132016626611332q0



    planett(5)%p(1) =  37018647652.25172q0  ! Mars
    planett(5)%p(2) = -211703392121.9987q0 
    planett(5)%p(3) = -5332148864.790916q0

    planett(5)%velocity(1) = 24764.87726268154q0        
    planett(5)%velocity(2) =  6378.058270749773q0 
    planett(5)%velocity(3) = -473.4151780243245q0




    planett(6)%p(1) = 489401410294.7523q0   ! Jupyter 
    planett(6)%p(2) = 562533468901.4261q0
    planett(6)%p(3) = -13283023128.02374q0

    planett(6)%velocity(1) = -10001.75006433559q0     
    planett(6)%velocity(2) =  9194.897567773275q0
    planett(6)%velocity(3) =  185.6192334712432q0



    planett(7)%p(1) = 1354644048313.250q0   ! Saturn 
    planett(7)%p(2) = -526994827771.3244q0
    planett(7)%p(3) = -44771791167.80645q0  


    planett(7)%velocity(1) = 2963.258037069641q0  
    planett(7)%velocity(2) = 8984.428483104208q0
    planett(7)%velocity(3) = -274.7875628875702q0

    

    planett(8)%p(1) = 1816861070124.457q0   !Uranus 
    planett(8)%p(2) = 2301334470958.584q0
    planett(8)%p(3) = -14990560295.87567q0

    planett(8)%velocity(1) = -5394.965755528056q0       
    planett(8)%velocity(2) =  3902.419055068386q0
    planett(8)%velocity(3) =  84.18312321275012q0


    planett(9)%p(1) =  4464172106275.985q0   ! Neptune 
    planett(9)%p(2) = -250429368100.1182q0
    planett(9)%p(3) =  -97724357471.04323q0

    planett(9)%velocity(1) = 268.6511540628924q0      
    planett(9)%velocity(2) =  5459.061702264910q0
    planett(9)%velocity(3) = -118.6506672844314q0



    planett(10)%p(1) = -1.113175726324003q11   ! Moon
    planett(10)%p(2) = 9.721142863263611q10
    planett(10)%p(3) =  -5.449223300203681q6

    planett(10)%velocity(1) = -1.929177591167254q4      
    planett(10)%velocity(2) =  -2.195669409085904q4
    planett(10)%velocity(3) = 1.249049169100225q1

    planett(11)%p(1) =  2.589356480215066q12   ! Pluto
    planett(11)%p(2) = -4.534121953364584q12
    planett(11)%p(3) =  -2.638193356204123q11

    planett(11)%velocity(1) = 4.848290480101916q3      
    planett(11)%velocity(2) =  1.483442530065416q3
    planett(11)%velocity(3) = -1.574946376548967q3

    planett(12)%p(1) =  4.895199773547107q11  
    planett(12)%p(2) = 5.629369089534231q11   ! IO
    planett(12)%p(3) =  -1.326710498124310q10

    planett(12)%velocity(1) = -2.665517430567016q4      
    planett(12)%velocity(2) =  1.415563698368067q4
    planett(12)%velocity(3) =  1.161540698989771q2

    planett(13)%p(1) =  4.889355082875217q11     ! EUROPA
    planett(13)%p(2) = 5.620588435455503q11   
    planett(13)%p(3) =  -1.330958536876112q10

    planett(13)%velocity(1) = -1.546706269066096q2      
    planett(13)%velocity(2) = -5.485138841100538q2 
    planett(13)%velocity(3) =  7.497879606478131q1 



    planett(14)%p(1) =  4.903145916266326q11     ! Ganymede
    planett(14)%p(2) = 5.630874979940245q11   
    planett(14)%p(3) =  -1.324903816225049q10

    planett(14)%velocity(1) = -1.564597580880949q4      
    planett(14)%velocity(2) = 1.851415585396170q4 
    planett(14)%velocity(3) =  4.614615739104782q2 


    planett(15)%p(1) =  4.890526414345648q11     !Callisto   
    planett(15)%p(2) = 5.606798236479037q11   
    planett(15)%p(3) =  -1.334576095003527q10

    planett(15)%velocity(1) = -1.949978659331871q3     
    planett(15)%velocity(2) = 7.735093433973734q3
    planett(15)%velocity(3) =  2.474741601007366q2


    
    
    planett(16)%p(1) =  1.35482130055603q12     !Mimas   
    planett(16)%p(2) = -5.269618491871567q11   
    planett(16)%p(3) =  -4.481063551679203q10     

    planett(16)%velocity(1) = -9.975384461854917q2     
    planett(16)%velocity(2) = 2.129107627829789q4
    planett(16)%velocity(3) =  -6.621978811452451q3

    planett(17)%p(1) =  1.354642489401001q12    !Tethys    
    planett(17)%p(2) = -5.267337153093208q11
    planett(17)%p(3) =  -4.49082899937742q10

    planett(17)%velocity(1) = -8.327037745298369q3      
    planett(17)%velocity(2) =  9.482893597245901q3
    planett(17)%velocity(3) = 8.031366314462001q2



    planett(18)%p(1) =  1.354301039006643q12   !Dione 
    planett(18)%p(2) = -5.268424827310864q11
    planett(18)%p(3) =  -4.481821623568347q10

    planett(18)%velocity(1) = -1.173442587952561q3      
    planett(18)%velocity(2) =  1.081806718202701q3
    planett(18)%velocity(3) =  4.263527690000798q3


    planett(19)%p(1) =  1.354603875975271q12   !Rhea 
    planett(19)%p(2) = -5.265282103838131q11
    planett(19)%p(3) =  -4.501408717515120q10

    planett(19)%velocity(1) = -5.465607558692203q3      
    planett(19)%velocity(2) =  8.735642111137274q3
    planett(19)%velocity(3) =  6.275394528207472q2

    planett(20)%p(1) =  1.353398408200447q12   !Titan  
    planett(20)%p(2) = -5.268334562528598q11
    planett(20)%p(3) =  -4.473091702955183q10

    planett(20)%velocity(1) = 2.447694738070457q3      
    planett(20)%velocity(2) =  4.213965630888995q3
    planett(20)%velocity(3) =  2.236824514016694q3

    planett(21)%p(1) =  1.351461596346124q12   !Iapetus    
    planett(21)%p(2) = -5.285043151899385q11
    planett(21)%p(3) =  -4.378158669224045q10

    planett(21)%velocity(1) = 4.427822743486593q3      
    planett(21)%velocity(2) =  6.187387843686237q3
    planett(21)%velocity(3) =  7.658595086595366q1

    planett(22)%p(1) =  4.463969919642986q12   !Triton Neptune!!   
    planett(22)%p(2) = -2.503165359347077q11     
    planett(22)%p(3) =  -9.745539266189398q10

    planett(22)%velocity(1) = 3.063548586374388q3      
    planett(22)%velocity(2) =  8.768002599697214q3
    planett(22)%velocity(3) =  5.890714089677549q2

    planett(23)%p(1) =  1.816435517653067q12   ! Titania - Uranus   
    planett(23)%p(2) = 2.301421392860621q12    
    planett(23)%p(3) =  -1.503629005833972q10

    planett(23)%velocity(1) = -5.673062370925212q3      
    planett(23)%velocity(2) =  4.464290398957321q3
    planett(23)%velocity(3) =  3.670204672685641q3

    planett(24)%p(1) =  1.816491709929089q12   ! Oberon - Uranus   
    planett(24)%p(2) = 2.301350657918233q12    
    planett(24)%p(3) =  -1.544244438785636q10

    planett(24)%velocity(1) = -7.742821176511318q3      
    planett(24)%velocity(2) = 4.690449419677916q3
    planett(24)%velocity(3) =  2.029995590424519q3
    
    planett(25)%p(1) =  1.816601921230485q12   ! Umbriel - Uranus   
    planett(25)%p(2) = 2.301388401562102q12    
    planett(25)%p(3) =  -1.500912503409970q10

    planett(25)%velocity(1) = -5.569775018556443q3      
    planett(25)%velocity(2) = 4.581402604313227q3
    planett(25)%velocity(3) =  4.709885484103392q3






    
    ! -------------------------------------------------

    planett(1)%mass = Sun%mass
    planett(2)%mass = Mercury%mass
    planett(3)%mass = Venus%mass
    planett(4)%mass = Earth%mass
    planett(5)%mass = Mars%mass
    planett(6)%mass = Jupiter%mass
    planett(7)%mass = Saturn%mass
    planett(8)%mass = Uranus%mass
    planett(9)%mass = Neptune%mass
    planett(10)%mass = Moon%mass
    planett(11)%mass = Pluto%mass
    planett(12)%mass = Io%mass
    planett(13)%mass = Europa%mass
    planett(14)%mass = Ganymede%mass
    planett(15)%mass = Callisto%mass
    planett(16)%mass = Mimas%mass
    planett(17)%mass = Tethys%mass
    planett(18)%mass = Dione%mass
    planett(19)%mass = Rhea%mass
    planett(20)%mass = Titan%mass
    planett(21)%mass = Iapetus%mass
    planett(22)%mass = Triton%mass
    planett(23)%mass = Titania%mass
    planett(24)%mass = Oberon%mass
    planett(25)%mass = Umbriel%mass
    

    !----------------------------------------------------

    planett(1)%gm = Sun%k
    planett(2)%gm = Mercury%k 
    planett(3)%gm = Venus%k
    planett(4)%gm = Earth%k 
    planett(5)%gm = Mars%k 
    planett(6)%gm = Jupiter%k
    planett(7)%gm = Saturn%k 
    planett(8)%gm = Uranus%k
    planett(9)%gm = Neptune%k
    planett(11)%gm = Pluto%k
    planett(10)%gm = Moon%k
    planett(12)%gm = Io%k
    planett(13)%gm = Europa%k
    planett(14)%gm = Ganymede%k
    planett(15)%gm = Callisto%k
    planett(16)%gm = Mimas%k
    planett(17)%gm = Tethys%k
    planett(18)%gm = Dione%k
    planett(19)%gm = Rhea%k
    planett(20)%gm = Titan%k
    planett(21)%gm = Iapetus%k
    planett(22)%gm = Triton%k
    planett(23)%gm = Titania%k
    planett(24)%gm = Oberon%k
    planett(25)%gm = Umbriel%k
    
    



   
    


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

    t = int(t/step,kind=16)

    

    energy0 = 0 
    test = 0 

    do k = 1, t 

        day = int(k*step,kind=16)
        

        if (day == 0) then 
            day = 1 
        end if 
        

        if (mod(day,86400)== 0 .and. day /= 0) then 
            write(*,*) "day : " , int(day/86400.0,kind=16) 
        end if 

        do i = 1, 25
            planett(i)%acc = 0 
        end do 


        do i = 1,25
            do j = 1 ,25

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

            do i = 1,25

                planett(i)%velocity(1) = planett(i)%velocity(1) + planett(i)%acc(1)*step/real(2)
                planett(i)%velocity(2) = planett(i)%velocity(2) + planett(i)%acc(2)*step/real(2)
                planett(i)%velocity(3) = planett(i)%velocity(3) + planett(i)%acc(3)*step/real(2)

                planett(i)%p(1) = planett(i)%p(1) + planett(i)%velocity(1)*step 
                planett(i)%p(2) = planett(i)%p(2) + planett(i)%velocity(2)*step 
                planett(i)%p(3) = planett(i)%p(3) + planett(i)%velocity(3)*step 

            end do

            
        else 
            do i = 1,25
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
            do i = 1, 25 
                planett(i)%energy = 0 
            end do 

            totale = 0 

            do i = 1,25
                do j = 1,25
                    
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


            do i = 1, 25
                totale = totale + planett(i)%energy  
            end do 

            if (abs((totale-energy0)/energy0)>test ) then 
                test = abs((totale-energy0)/energy0)
                write(22,"(i14,x,f60.20)") day, test

            end if 

            test1 = abs((totale-energy0)/energy0) 
            write(23,"(i14,x,f60.20)") day, test1

        else if (k==1) then
            do i = 1, 25 
                planett(i)%energy = 0 
            end do 

            do i = 1,25
                do j = 1,25
                    
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


            do i = 1, 25
                energy0 = energy0 + planett(i)%energy  
            end do  

            

            write(23,"(i14,x,f60.20)") day, real(0)

            
        
            

            

        end if  
            
       
        


        

    end do 

    write(*,*) "The end! You may check all the data"

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
    

end program many_solar