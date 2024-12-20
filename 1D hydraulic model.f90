    program inercial

    implicit none

    !************************************************************************
    ! DECLARATION OF VARIABLES
    
    ! CHECK INPUT FILE COTA_RIO!! - Setting the model up / initial conditions - define Z properly
    
    ! Input Variables
    
    real(kind=8)                            :: amp, med, def                ! lateral inflow parameters based on MANICORE
    integer,dimension(:),allocatable        :: CATID, dia, mes, ano         ! reach number and simulation dates
    integer,dimension(:,:),allocatable      :: cotas_fp                     ! ground level for floodplain evaluation
    real(kind=8),dimension(:),allocatable   :: W, L, Depth, L1, L2                  ! river's width, lenght and depth
    real(kind=8),dimension(:),allocatable   :: Area_M, Xcen, Ycen, Zo       ! accumulated area(%), position (lon and lat) and elevation
    real(kind=8),dimension(:),allocatable   :: Y_input, Q_input, Q_input_Wubu, Q_input_Longmen            ! boundary conditions  
    real(kind=8),dimension(:),allocatable   :: Area_fp                      ! flooded area (km2) 
    
    ! Output variables
    
    real(kind=8),dimension(:,:),allocatable :: QTUDO, HTUDO                 ! dayly inflow and depth (NOT IN USE)
    
    
    ! Model variables
    
    character,dimension(6)                  :: lixo_cab                     ! reading and discarding headings
    integer                                 :: iC, iT, iCA, io, i           ! counting variables
    integer                                 :: nT, nC                     ! days to simulate and number of reaches
    integer                                 :: nCA                          ! size of cota_area ('level - flooded area') file
    real(kind=8)                            :: lixo_num, count              ! discarding useless info and counting variable                    
    real(kind=8)                            :: dt, alfa, dtsum, dtmin       ! time step, courant coeficient, sum of dt and minimum dt
    real(kind=8)                            :: Ac, Aprop, phi               ! difference of accum Area and lag for lateral inflow
    real(kind=8)                            :: interp                       ! to interpolate between boundary conditions
    real(kind=8),parameter                  :: g=9.81, day=86400.!, n=0.030  ! gravity acceleration, day to seconds, Manning Coeficient
    real(kind=8),parameter                  :: pi=3.1415926                 ! pi value
    real(kind=8),dimension(:),allocatable   :: Afp                          ! flooded area
    real(kind=8),dimension(:),allocatable   :: At, Rh                       ! x-section area and hydraulic radius
    real(kind=8),dimension(:),allocatable   :: h, Z, Y, hmed                ! height variables
    real(kind=8),dimension(:),allocatable   :: Qlat, Q                      ! inflow variables (main and lateral)
    real(kind=8),dimension(:),allocatable   :: n                            ! n-Manning Coefficient Values

    ! Data Assimilation variables
    
    integer                                 :: nO, nDA, nDAmax              ! count variables
    character(LEN=46)                       :: DAfile                       ! file name to be opened
    real,dimension(6)                       :: lixo                         ! to discard info
    real(kind=8)                            :: dh!, var_H                    ! water level difference between obs and model
    real(kind=8),dimension(:,:),allocatable :: DA, deltaZ                   ! matrix with all altimetry info and difference of Height (Yobs-Y)
    real(kind=8),dimension(:),allocatable   :: var_Z, var_n, K, var_H              ! bottom level variance
    
        
        
    !************************************************************************
    ! OPENING INPUT FILES TO COUNT nT, nC AND nCA
    
    
    ! open files to count nT E nC

    open(4,FILE='input\cota_rio.txt')                             ! river altitude
    open(6,FILE='input\niveis_saida.txt')                         ! levels at the river mouth
    
    ! reading cota_rio.txt to count number of reaches (nC)
    
    io=0
    nC=0
    do while (io==0)              
        read(4,*,IOSTAT=io), lixo_num
        nC=nC+1
    end do
    nC=nC-1
    PRINT*,''
    print*,'n reaches = ', nC
        
    close(4)
    
    ! reading niveis_saida.txt to count period of simulation (nT)
    
    io=0
    nT=0
    do while (io==0)              
        read(6,*,IOSTAT=io), lixo_num
        nT=nT+1
    end do
    nT=nT-1
    print*,''
    print*,'n time = ', nT
    print*,''

    PAUSE
    
    close(6)


    ! reading cota_area_fp.txt to count nCA for floodplain evaluation
    
    open(15,FILE='input\cota_area_rp3.txt')                         ! levels vs flooded area file
    read(15,*), lixo_cab(1:4)    
    
    io=0
    nCA=0
    do while (io==0)              
        read(15,*,IOSTAT=io), lixo(1:4)
        nCA=nCA+1
    end do
    nCA=nCA-1
    
    print*,''
    print*,'n cota_area= ', nCA
    print*,''

    PAUSE
    
    close(15)

    
    !************************************************************************
    ! ALLOCATING VARIABLES
    
    allocate(dia(nT))
    allocate(mes(nT))
    allocate(ano(nT))
    allocate(Y_input(nT))    
    allocate(Q_input(nT))
    allocate(Q_input_Wubu(nT))
    allocate(Q_input_Longmen(nT))
    
    allocate(CatID(nC))
    allocate(W(nC))
    allocate(L(nC))
    allocate(L1(nC))
    allocate(L2(nC))
    allocate(Depth(nC))
    allocate(Area_M(nC))
    allocate(Xcen(nC))
    allocate(Ycen(nC))
    allocate(Zo(nC))
    allocate(At(nC-1))
    allocate(Rh(nC-1))
    allocate(hmed(nC-1))
    allocate(H(nC))
    allocate(Z(nC))
    allocate(Y(nC))
    allocate(Qlat(nC))
    allocate(Q(nC))
    allocate(Afp(nC))
    allocate(var_Z(nC))
    allocate(n(nC))
    allocate(var_n(nC))
    allocate(K(nC))
    allocate(var_H(nC))
    
    allocate(Area_fp(nCA))
    allocate(cotas_fp(nCA,3))
    
    allocate(QTUDO(nC,nT))
    allocate(HTUDO(nC,nT))    

    IF (1==0) then !LQ
    !************************************************************************
    ! DATA ASSIMILATION
    
    ! Opening files

    open(10,FILE='input\alt_obs3.txt')                             ! observation sites (Adrien)
    
    ! reading alt_obs.txt to count number of altimetry observation sites (nO)
    ! also to count number of observations during simulation period (nDA)
    
    io=0
    nO=0
    nDAmax=0
    
    do while (io==0)

        nO=nO+1
        read(10,2,IOSTAT=io), DAfile, lixo(1:2), iC             !
        
        print*,DAfile
        
        i=0
        nDA=0
        open(10+nO,FILE='input/'//DAfile//'.cal')
        do while (i==0)
            read(10+nO,*,IOSTAT=i), lixo(1:6)
            nDA=nDA+1
        end do
        
        close(10+nO)
        
        if (nDA>nDAmax) then
            nDAmax=nDA
        end if
        
    end do
    
    nO=nO-1
    nDA=nDAmax
    
    print*,''
    print*,'n altimetry sites', nO 
    print*,'n max obs', nDA
    print*,''
    
    PAUSE
    close(10)
    
    ! allocating data assimilation variables
    
    allocate(DA(nDA,nO*4))
    allocate(deltaZ(nO,3))                      ! st - position compared to reaches, nd - flag to interpolate, rd - delta Z indeed (Yobs-Y)
    
    end if !LQ
    nDA = 31 !38,14
    nO = 14
    allocate(DA(nDA,5))
    allocate(deltaZ(nO,3)) 
    !************************************************************************
    ! OPENING INPUT FILES AND STORING VARIABLES
    
    ! opening files
    
    open(1,FILE='input\mini_rp2.txt')                              ! reaches info (CatID, Area_M, Lenght, Depth, Xcen and Ycen)
    open(2,FILE='input\coef_hidro_lat1.txt')                       ! lateral inflow info (amp and mean for a sinusoidal function)
    open(3,FILE='input\largura_rio.txt')                          ! river width
    open(4,FILE='input\cota_rio.txt')                             ! river altitude
    open(5,FILE='input\hidrograma_entr.txt')                      ! input flow
    open(6,FILE='input\niveis_saida.txt')                         ! levels at the river mouth
    open(15,FILE='input\cota_area_rp3.txt')                        ! floodplain levels and area  
    
    ! Reading and storing variables
    
    read(1,*), lixo_cab
    read(2,*), lixo_cab(1:2)
    read(15,*), lixo_cab(1:4)
    
    !read(2,*), Amp, Med
    
    do iC=1,nC
        read(1,*), CatID(iC), Area_M(iC), L(iC), Depth(iC), Xcen(iC), Ycen(iC), L1(iC), L2(iC)
        read(3,*), W(iC)
        read(4,*), Zo(iC)
    end do
    
    do iT=1,nT
        read(2,*), Q_input_Wubu(iT), Q_input_Longmen(iT)
        read(5,*), dia(iT), mes(iT), ano(iT), Q_input(iT) !LQ lixo_num, column 4
        read(6,*), Y_input(iT)
    end do
    
    
    do iCA=1,nCA
        read(15,*) cotas_fp(iCA,:), area_fp(iCA)
        area_fp(iCA)=area_fp(iCA)*1000000
    end do
    
    
    close(1)
    close(2)
    close(3)
    close(4)
    close(5)
    close(6)
    close(15)
    !print*, L1
    !print*, L2
    PAUSE
    
    !LQ ! Reading and storing DATA ASSIMILATION variables
    
    IF (1==1) then !LQ
    open(10,FILE='input\alt_obs3.txt')                             ! observation gauges
    
    do io=1,nDA
        read(10,*), DA(io,1:5)
    end do
    
    !DA=-1
    deltaZ=-100
    
    !do io=1,nO
    
        !read(10,2), DAfile, lixo(1:3)
        
        !lixo_num=9999
        !iC=0
        !do i=1,nC
            !count=sqrt((Ycen(i)-lixo(1))**2+(Xcen(i)-lixo(2))**2)
            !if (count<lixo_num) then
                !lixo_num=count
                !iC=i
            !end if
        !end do
        
        !deltaZ(io,1)=iC
        
        !print*, iC
        
        !i=0
        !iC=0
        !open(10+io,FILE='input/'//DAfile//'.cal')
        !do while (i==0)
            !iC=iC+1
            !read(10+io,*,IOSTAT=i),DA(iC,((io-1)*4+1):(4*io)),lixo(1:2)       ! saving date (year month and day) and water level
        !end do
        
    close(10)
        
    !end do
    end if !LQ
    print*, DA
    PAUSE
    
    !************************************************************************
    ! SETTING THE MODEL UP
    
    ! Parameters
    
    alfa=0.7
    def=10
    
    !Initial Conditions
    
    Z=Zo-Depth
    !Z=Zo
    H=Y_input(1)-Z(1)
    Y=Z+H
    Qlat(nC)=0
    var_Z=100
    var_n=0.01*0.01
    var_H=100
    n=0.03
    K=0
    
    do iC=1,nC-1
        hmed(iC)=max(Y(iC),Y(iC+1))-max(Z(iC),Z(iC+1))
        At(iC)=(W(iC)+W(iC+1))*0.667*hmed(iC)/2.
        Rh(iC)=2.*At(iC)/(4.*sqrt(((W(iC)+W(iC+1))/6.0*(W(iC)+W(iC+1))/6.0)+(hmed(iC)*hmed(iC)))+((W(iC)+W(iC+1))*0.333))
        L(iC)=L(iC)!*1000 !LQ His length is in meters
    end do
    
    L(nC)=L(nC)!*1000 !LQ His length is in meters
    
    dt=alfa*L(nC)/SQRT(g*H(1))
    
    count=0.
    dtsum=0.
    
    ! Opening from
    Open(7,FILE='output/QTUDO.txt')
    Open(8,FILE='output/HTUDO.txt')
    Open(9,FILE='output/cota_de_fundo.txt')
    Open(20,FILE='output/var_Z.txt')
    open(22,FILE='output/Kalman.txt')
    Open(23,FILE='output/var_H.txt')
    
    write(20,1),var_Z
    
    do iT=1,nT-1
        
        count=count+1.
        
        ! flooded area
            
        Afp=0.
        i=1
        
        do iCA=1,nCA-1
            if (cotas_fp(iCA,3)<=Y(i)) then
                Afp(i)=Area_fp(iCA)
            end if
            if (cotas_fp(iCA+1,1)/=cotas_fp(iCA,1)) then
                i=i+1
            end if
        end do

        
        ! Calculate whitin a day
        
        do while (dtsum<=(count*day))
            
            ! boundary conditions
            
            interp=(dtsum-iT*day)/(iT*day)
            Q(nC)=Q_input(iT)+interp*(Q_input(iT+1)-Q_input(iT))        
            Y(1)=Y_input(iT)+interp*(Y_input(iT+1)-Y_input(iT))
            H(1)=Y(1)-Z(1)
            !print*, Q(nC),Qlat
            !PAUSE

            do iC=1,nC-1                        ! calcular as coisas
                
                i=nC-iC                         ! changing order to start calculation from upstream
                
                IF (1==1) then !LQ
                ! calculo Qlat
                IF (i>=9) then
                Ac=Area_M(i)-Area_M(i+1)
                Aprop=Ac/Area_M(9)             ! 9 is the WUBU observation spot
                !phi=def*(nC-i)/(nC-91)
                !Qlat(i)=med*Aprop+amp*Aprop*sin(2.*pi*(dtsum+phi*day)/(365*day))   
                Qlat(i) = Aprop*(Q_input_Wubu(iT+1)-Q_input(iT))
                endif
                
                IF (i<9) then
                Ac=Area_M(i)-Area_M(i+1)
                Aprop=Ac/(Area_M(1)-Area_M(9))             ! 1 is the LONGMEN observation spot
                !phi=def*(nC-i)/(nC-91)
                !Qlat(i)=med*Aprop+amp*Aprop*sin(2.*pi*(dtsum+phi*day)/(365*day))   
                Qlat(i) = Aprop*(Q_input_Longmen(iT+1)-Q_input_Wubu(iT))
                endif
                
                !IF (Qlat(i) < 0) then
                !Qlat(i)=1
                !endif
                endif !LQ
                !Qlat=1
                !print*, H(i+1)
                !PAUSE
                ! Calculo Q e H
                
                Q(i)=(Q(i)-dt*g*At(i)*L1(i)*(Y(i)-Y(i+1))/(L(i)+L(i+1)))/(1+dt*g*abs(Q(i))*n(i)*n(i)/(At(i)*Rh(i)**(4./3.)))
                !Q(i)=(Q(i)-dt*g*At(i)*(Y(i)-Y(i+1))/(L2(i)))/(1+dt*g*abs(Q(i))*n(i)*n(i)/(At(i)*Rh(i)**(4./3.)))
                H(i+1)=H(i+1)-dt*(Q(i)-Q(i+1)-Qlat(i+1))/(L(i+1)*W(i+1)+Afp(i+1))
                Y(i+1)=H(i+1)+Z(i+1)
                
                !print*, iT, iC, dt*(Q(i)-Q(i+1)-Qlat(i+1)),Q(i),Q(i+1),Qlat(i+1)
                !PAUSE
            end do
            
            ! cross section (Area, Hydraulic radius and height)
            
            do iC=1,nC-1
                !rect
                !hmed(iC)=max(Y(iC),Y(iC+1))-max(Z(iC),Z(iC+1))
                !At(iC)=(W(iC)+W(iC+1))/2.*hmed(iC)
                !Rh(iC)=2.*At(iC)/(4.*hmed(iC)+W(iC)+W(iC+1))
                !tri
                !IF (iC==19 .or. iC==18 .or. iC==17 .or. iC==16 .or. iC==11 .or. iC==10 .or. iC==9 .or. iC==6 .or. iC==5 .or. iC==4 .or. iC==3 .or. iC==2) then
                !hmed(iC)=max(Y(iC),Y(iC+1))-max(Z(iC),Z(iC+1))
                !At(iC)=(W(iC)+W(iC+1))/2.*hmed(iC)/2.
                !Rh(iC)=2.*At(iC)/(2.*sqrt((W(iC))/2.*((W(iC))/2.)+(hmed(iC)*hmed(iC)))+(2.*sqrt((W(iC+1))/2.*((W(iC+1))/2.)+(hmed(iC)*hmed(iC)))))
                !endif
                !!esp
                !IF (iC==14 .or. iC==13) then
                !hmed(iC)=max(Y(iC),Y(iC+1))-max(Z(iC),Z(iC+1))
                !At(iC)=(W(iC)+W(iC+1))/4.*hmed(iC)*3.1415926/2.
                !Rh(iC)=2.*At(iC)/((3.1415926*hmed(iC)+(2.*((W(iC))/2.-hmed(iC))))+(3.1415926*hmed(iC)+(2.*((W(iC+1))/2.-hmed(iC)))))
                !endif
                !!tri1
                !IF (iC==12 .or. iC==7) then
                !hmed(iC)=max(Y(iC),Y(iC+1))-max(Z(iC),Z(iC+1))
                !At(iC)=W(iC)*hmed(iC)/2.
                !Rh(iC)=(At(iC))/(2.*sqrt((W(iC)/2.*W(iC)/2.)+(hmed(iC)*hmed(iC))))
                !endif
                !!esp
                !IF (iC==15 .or. iC==8 .or. iC==1) then
                !hmed(iC)=max(Y(iC),Y(iC+1))-max(Z(iC),Z(iC+1))
                !At(iC)=W(iC)/2.*hmed(iC)*3.1415926/2.
                !Rh(iC)=(At(iC))/(3.1415926*hmed(iC)+(2.*(W(iC)/2.-hmed(iC))))
                !endif
                
                !rect
                !hmed(iC)=max(Y(iC),Y(iC+1))-max(Z(iC),Z(iC+1))
                !At(iC)=(W(iC)+W(iC+1))/2.*hmed(iC)
                !Rh(iC)=2.*At(iC)/(4.*hmed(iC)+W(iC)+W(iC+1))
                !tri
                !IF (iC==18 .or. iC==17 .or. iC==16 .or. iC==11 .or. iC==10 .or. iC==9 .or. iC==6 .or. iC==5 .or. iC==4 .or. iC==3 .or. iC==2) then
                !hmed(iC)=max(Y(iC),Y(iC+1))-max(Z(iC),Z(iC+1))
                !At(iC)=(W(iC)+W(iC+1))/2.*hmed(iC)/2.
                !Rh(iC)=2.*At(iC)/(2.*sqrt((W(iC))/2.*((W(iC))/2.)+(hmed(iC)*hmed(iC)))+(2.*sqrt((W(iC+1))/2.*((W(iC+1))/2.)+(hmed(iC)*hmed(iC)))))
                !endif
                !!tra
                !IF (iC==14 .or. iC==13) then
                hmed(iC)=max(Y(iC),Y(iC+1))-max(Z(iC),Z(iC+1))
                At(iC)=(W(iC)+W(iC+1))*0.667*hmed(iC)/2.
                Rh(iC)=2.*At(iC)/(4.*sqrt(((W(iC)+W(iC+1))/6.0*(W(iC)+W(iC+1))/6.0)+(hmed(iC)*hmed(iC)))+((W(iC)+W(iC+1))*0.333))
                !endif
                !!tri1
                !IF (iC==12 .or. iC==7) then
                !hmed(iC)=max(Y(iC),Y(iC+1))-max(Z(iC),Z(iC+1))
                !At(iC)= W(iC)*hmed(iC)/2.
                !Rh(iC)=(W(iC)*hmed(iC)/2.)/(2.*sqrt((W(iC)/2.*W(iC)/2.)+(hmed(iC)*hmed(iC))))
                !endif
                !!tra
                !IF (iC==15 .or. iC==8 .or. iC==1) then
                !hmed(iC)=max(Y(iC),Y(iC+1))-max(Z(iC),Z(iC+1))
                !At(iC)=W(iC)*0.667*hmed(iC)
                !Rh(iC)=(W(iC)*0.667*hmed(iC))/(2.*sqrt(((W(iC))/6.*(W(iC))/6.)+(hmed(iC)*hmed(iC)))+((W(iC))*0.333))
                !endif
                
                !rect
                !hmed(iC)=max(Y(iC),Y(iC+1))-max(Z(iC),Z(iC+1))
                !At(iC)=(W(iC)+W(iC+1))/2.*hmed(iC)
                !Rh(iC)=2.*At(iC)/(4.*hmed(iC)+W(iC)+W(iC+1))
                !tri
                !IF (iC==19 .or. iC==18 .or. iC==17 .or. iC==16 .or. iC==11 .or. iC==10 .or. iC==9 .or. iC==6 .or. iC==5 .or. iC==4 .or. iC==3 .or. iC==2) then
                !hmed(iC)=max(Y(iC),Y(iC+1))-max(Z(iC),Z(iC+1))
                !At(iC)=(W(iC)+W(iC+1))/2.*hmed(iC)/2.
                !Rh(iC)=2.*At(iC)/(2.*sqrt((W(iC))/2.*((W(iC))/2.)+(hmed(iC)*hmed(iC)))+(2.*sqrt((W(iC+1))/2.*((W(iC+1))/2.)+(hmed(iC)*hmed(iC)))))
                !endif
                !!esp
                !IF (iC==14 .or. iC==13) then
                !hmed(iC)=max(Y(iC),Y(iC+1))-max(Z(iC),Z(iC+1))
                !At(iC)=(W(iC)+W(iC+1))/4.*hmed(iC)*3.1415926/2.
                !Rh(iC)=2.*At(iC)/((3.1415926*hmed(iC)+(2.*((W(iC))/2.-hmed(iC))))+(3.1415926*hmed(iC)+(2.*((W(iC+1))/2.-hmed(iC)))))
                !endif
                !!tri1
                !IF (iC==12 .or. iC==7) then
                !hmed(iC)=max(Y(iC),Y(iC+1))-max(Z(iC),Z(iC+1))
                !At(iC)=(W(iC)*hmed(iC)/2.+(W(iC+1)/2.*hmed(iC)*3.1415926/2.))/2.
                !!At(iC) = W(iC)*hmed(iC)/2.
                !Rh(iC)=2.*At(iC)/(2.*sqrt((W(iC)/2.*W(iC)/2.)+(hmed(iC)*hmed(iC)))+(3.1415926*hmed(iC)+(2.*(W(iC+1)/2.-hmed(iC)))))
                !endif
                !!esp
                !IF (iC==15 .or. iC==8 .or. iC==1) then
                !hmed(iC)=max(Y(iC),Y(iC+1))-max(Z(iC),Z(iC+1))
                !!At(iC)=W(iC)/2.*hmed(iC)*3.1415926/2.
                !At(iC)=(W(iC)/2.*hmed(iC)*3.1415926/2.+(W(iC+1)*hmed(iC)/2.))/2.
                !Rh(iC)=2.*At(iC)/(3.1415926*hmed(iC)+(2.*(W(iC)/2.-hmed(iC)))+ (2.*sqrt((W(iC+1)/2.*W(iC+1)/2.)+(hmed(iC)*hmed(iC)))))
                !endif


            end do
            
            dt=day
            
            do iC=1,nC                    ! Estimar proximo dt
                dtmin=alfa*L(iC)/sqrt(g*H(iC))
                if (dtmin<=dt) then
                    dt=dtmin
                end if
            end do
    
            dtsum=dtsum+dt
            !print*, dtsum, (count*day)
            !PAUSE
        end do
        
        ! Data assimilation section (CHOOSE THE ASSIMILATION METHOD)
        
        ! interpolation
        
        !do io=1,nDA
        !    do iC=1,nO
        !        if (DA(io,((iC-1)*4+1))==ano(iT) .and. DA(io,((iC-1)*4+2))==mes(iT) .and. DA(io,((iC-1)*4+3))==dia(iT)) then
        !            deltaZ(iC,2)=DA(io,4*iC)-Y(deltaZ(iC,1))
        !            deltaZ(iC,3)=1
        !        end if
        !    end do
        !end do
        !
        !if (sum(deltaZ(:,3))==16) then
        !    call interpol(Z,deltaZ,nC,nO)
        !    deltaZ(:,2:3)=0
        !end if
        
        ! data assimilation indeed
        !LQ ! no DATA ASSIMILATION FOR 1st test
        
        IF (1==1) then !LQ
        do io=1,nDA
            !do iC=nO,1,-1
                if (DA(io,1)==ano(iT) .and. DA(io,2)==mes(iT) .and. DA(io,3)==dia(iT)) then
                    
                    !print*,'DA ',int(deltaZ(iC,1))
                    
                    dh=DA(io,4)-Y(DA(io,5))
                    
                    call data_assimilation(Z,H,Y,n,var_Z,var_n,L,dh,nC,int(DA(io,5)),Q,W,K,var_H)
               
                    write(20,1),var_Z
                    Y=Z+H                    
         
                    write(22,1),K
                    write(23,*),var_H(int(DA(io,5)))
                    
                    !print*, Y(deltaZ(iC,1)),H(deltaZ(iC,1))
                    
                    !pause
        
                end if
            !end do
        end do
        end if !LQ
        
        ! End of data assimilation section
        !print*, H
        write(7,1),Q
        write(8,1),H
        write(9,1),Z
    
        print*,iT
        
    end do
    
    close(7)
    close(8)
    close(9)
    
    ! Writing an altimetry data file that contains all observations together
    
    Open(11,FILE='output/altimetry_data.txt')
    Open(12,FILE='output/position_alt.txt')
    
    do iC=1,nO
        write(12,*),int(deltaZ(iC,1))
    end do
        
    do iC=1,nDA
        write(11,3),DA(iC,:)        
    end do    
        
    print*, 'FINISH'
    
    open(30,FILE='output/manning.txt')
    do iC=1,nC
        write(30,*),n(iC)
    end do
    
    
    PAUSE
    
1   FORMAT(211f10.2)
2   FORMAT(a46,f25.3,f9.3,i8)    
3   FORMAT(64f10.3)
    
    end program inercial
