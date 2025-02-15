! data assimilation

dhdhx2(position)=0  !new
K=0.

if (dh>=0.1) then
    m=int(dh/0.1)
else
    m=1
end if

!do o=1,m
!
!    acum=0.
!    var_H(1)=0.0025
!    
!    do i=2,nC
!        !dhdhx=exp(-3.3*(Y(i)-Y(i-1))/H(i-1))        
!        !fx(i)=1+(H(i-1)/(3.3*(Y(i)-Y(i-1))))*(dhdhx-1)
!    
!        Sf=((n(i)+n(i-1))/2.*(Q(i)+Q(i-1))/2./((H(i)+H(i-1))/2.)**(5./3.)/((W(i)+W(i-1))/2.))**2
!    
!        !Sf=(n(i-1)*Q(i-1)/H(i-1)**(5./3.)/W(i-1))**2
!        !Sf=2*(Y(i)-Y(i-1))/(L(i)+L(i-1))
!        !dhdhx2(i-1)=-3.3*(Y(i)-Y(i-1))/H(i-1)                   !new
!    
!        dhdhx2(i-1)=-3.3*Sf*(L(i-1)+L(i))/(H(i-1)+H(i))                   !newer
!    
!        !fx(i)=1+(H(i-1)/(3.3*(Y(i)-Y(i-1))))*(exp(dhdhx2(i-1))-1)
!    
!        fx(i)=1.+(H(i-1)+H(i))/(3.3*Sf*(L(i-1)+L(i)))*(exp(dhdhx2(i-1))-1.)
!        
!        var_H(i)=(9./25.*H(i)*H(i)*(std_q*std_q+std_w*std_w+var_n(i)/(n(i)*n(i)))+var_Z(i))*fx(i)*fx(i)+exp(dhdhx2(i-1))**2*var_H(i-1)
!        
!        !var_H(i)=(var_Z(i))*fx(i)*fx(i)+exp(dhdhx2(i-1))**2*var_H(i-1)
!        
!        !L2=0.7*H(i)/sf
!        !var_H2=(9./25.*H(i)*H(i)*(std_q*std_q+std_w*std_w+var_n(i)/(n(i)*n(i)))+var_Z(i))!*fx(i)*fx(i)!+exp(dhdhx2(i-1))**2*var_H
!        !print*,var_H2,fx(i),(exp(dhdhx2(i-1))**2)
!        !print*,H(i),var_H
!        !print*,Sf,(Y(i)-Y(i-1)),H(i)
!        !print*,fx(i),L2,sf*1000
!    end do
!        
!        !i=position
!        !Sf=((n(i)+n(i-1))/2.*(Q(i)+Q(i-1))/2./((H(i)+H(i-1))/2.)**(5./3.)/((W(i)+W(i-1))/2.))**2
!        !dhdhx2(i-1)=-3.3*Sf*(L(i-1)+L(i))/(H(i-1)+H(i))                   !newer
!        !fx(i)=1.+(H(i-1)+H(i))/(3.3*Sf*(L(i-1)+L(i)))*(exp(dhdhx2(i-1))-1.)
!        !var_H(i)=9./25.*H(i)*H(i)*(std_q*std_q+std_w*std_w+var_n(i)/(n(i)*n(i)))+(var_Z(i))*fx(i)*fx(i)+exp(dhdhx2(i-1))**2*var_H(i-1)
!        !
!        !dhdhx2(position)=0  !new
!    
!    !do i=position+1,nC
!    !
!    !    Sf=((n(i)+n(i-1))/2.*(Q(i)+Q(i-1))/2./((H(i)+H(i-1))/2.)**(5./3.)/((W(i)+W(i-1))/2.))**2
!    !
!    !    dhdhx3=-3.3*Sf*(L(i-1)+L(i))/(H(i-1)+H(i))                   !newer
!    !
!    !    fx(i)=1.+(H(i-1)+H(i))/(3.3*Sf*(L(i-1)+L(i)))*(exp(dhdhx3)-1.)
!    !    !var_H(i)=(9./25.*H(i)*H(i)*(std_q*std_q+std_w*std_w+var_n(i)/(n(i)*n(i)))+var_Z(i))*fx(i)*fx(i)+exp(dhdhx2(i-1))**2*var_H(i-1)
!    !    
!    !end do
!
!    do i=2,position-1
!        !dhdhx3=exp(-3.3*(Y(position)-Y(i))*2/(H(i)+H(position)))
!        !dhdhx4=exp(-3.3*(Y(position)-Y(i))/(H(i)))
!    
!        dhdhx=exp(sum(dhdhx2(i:position-1)))  !new
!    
!        ! Updating Z
!        r=exp(-sum(L(i:position-1))*3/100000)
!        covhzx=fx(i)*var_Z(i)*dhdhx+sqrt(var_Z(i)*var_Z(position))*fx(position)*r
!        K(i)=covhzx/(var_H(position)+std_obs*std_obs)
!        Z(i)=Z(i)+K(i)*dh/real(m)
!        var_Z(i)=var_Z(i)-K(i)*covhzx/real(m)
!        acum=fx(i)*K(i)*dh/real(m)+acum*exp(dhdhx2(i-1))
!        Y(i)=Y(i)+acum
!        H(i)=Y(i)-Z(i)
!        
!                    
!        
!        !print*,acum, Y(i), Z(i)
!        !print*,fx(i),K(i),var_Z(i)
!    
!        !! Updating n-Manning
!        !dhdnx=3./5.*H(i)/n(i)
!        !covhzx=fx(i)*var_n(i)*dhdhx*dhdnx
!        !K(i)=covhzx/(var_H+std_obs*std_obs)
!        !print*,fx(i),dhdnx
!        !n(i)=n(i)+K(i)*dh
!        !var_n(i)=var_n(i)-K(i)*covhzx
!        !!print*,covhzx,K(i),var_n(i)
!
!    end do
!    
!    !position
!    i=position    
!    
!    dhdhx=1  !new
!    
!    ! Updating Z
!    r=0
!        covhzx=fx(i)*var_Z(i)*dhdhx+sqrt(var_Z(i)*var_Z(position))*fx(position)*r
!        K(i)=covhzx/(var_H(position)+std_obs*std_obs)
!        Z(i)=Z(i)+K(i)*dh/real(m)
!        var_Z(i)=var_Z(i)-K(i)*covhzx/real(m)
!        acum=fx(i)*K(i)*dh/real(m)+acum*exp(dhdhx2(i-1))
!        Y(i)=Y(i)+acum
!        H(i)=Y(i)-Z(i)
!    
!    do i=position+1,nC
!    
!         !Updating Z
!        dhdhx=0
!    
!        ! Updating Z
!        r=exp(-sum(L(position+1:i))*3/100000)
!        covhzx=fx(i)*var_Z(i)*dhdhx+sqrt(var_Z(i)*var_Z(position))*fx(position)*r
!        K(i)=covhzx/(var_H(position)+std_obs*std_obs)
!        Z(i)=Z(i)+K(i)*dh/real(m)
!        var_Z(i)=var_Z(i)-K(i)*covhzx/real(m)
!        acum=fx(i)*K(i)*dh/real(m)+acum*exp(dhdhx2(i-1))
!        Y(i)=Y(i)+acum
!        H(i)=Y(i)-Z(i)
!    
!    end do
!    !
!    print*,K(position-2),Z(position-2),Y(position-2)
!
!end do
!
!do i=2,nC
!    if(Y(i-1)>=Y(i)) then
!        H(i)=H(i)+Y(i-1)-Y(i)
!        Y(i)=Y(i-1)
!    end if
!end do
!Y=H+Z
!!print*,covhzx,K(position),var_H
   
!print*, H(position),Y(position)

do o=1,m

    acum=0.
    
    ! Calculating backwater effect and reduction factor
    
    Sf=(n(1)*Q(1)/H(1)**(5./3.)/W(1))**2.
    !Sf=2*(Y(2)-Y(1))/(L(2)+L(1))
    
    fx(1)=1.+H(1)/(3.3*Sf*L(1))*(exp(-3.3*Sf*L(1)/H(1))-1.)
    !fx2(1)=1.+H(1)/(3.3*Sf2*L(1))*(exp(-3.3*Sf2*L(1)/H(1))-1.)
    
    do i=2,nC
        
        !Sf=2*(Y(i)-Y(i-1))/(L(i-1)+L(i))
        Sf=((n(i)+n(i-1))/2.*Q(i-1)/((H(i)+H(i-1))/2.)**(5./3.)/((W(i)+W(i-1))/2.))**2
        !fx2(i)=1.+H(i)/(3.3*Sf2*L(i))*(exp(-3.3*Sf2*L(i)/H(i))-1.)
        fx(i)=1.+H(i)/(3.3*Sf*L(i))*(exp(-3.3*Sf*L(i)/H(i))-1.)

        dhdhx2(i-1)=-3.3*Sf*(L(i-1)+L(i))/(H(i-1)+H(i))                   !newer

    end do
    
    ! Calculating Var_H
    
    var_H(1)=0.0025
    i=2
    var_H(i)=(9./25.*H(i)*H(i)*(std_q*std_q+std_w*std_w+std_n*std_n)+var_Z(i))*fx(i)*fx(i)+exp(2.*dhdhx2(i-1))*var_H(i-1)
    
    do i=3,position
        
        var_H(i)=(9./25.*H(i)*H(i)*(std_q*std_q+std_w*std_w+std_n*std_n)+var_Z(i))*fx(i)*fx(i)+exp(2.*dhdhx2(i-1))*var_H(i-1)
        
        covQ=0.
        covW=0.
        covn=0.
        covZ=0.
        
        do j=2,i-1
            
            dhdhx=exp(sum(dhdhx2(j:i-1)))
            rZ=exp(-sum(L(j:(i-1)))*2./100000)
            rn=rZ
            rW=rZ
            
            covQ=covQ+dhdhx*(9./25.*H(i)*H(j)*std_q*std_q*fx(i)*fx(j))*1
            covW=covW+dhdhx*(9./25.*H(i)*H(j)*std_w*std_w*fx(i)*fx(j))*rW
            covn=covn+dhdhx*(9./25.*H(i)*H(j)*std_n*std_n*fx(i)*fx(j))*rn
            covZ=covZ+dhdhx*(sqrt(var_Z(i)*var_Z(j))*fx(i)*fx(j))*rZ    
            
        end do
        
        var_H(i)=var_H(i)+2.*(covQ+covW+covn+covZ)
    
    end do
            
    ! Calculating the remaining (covariance, Kalman Gain, etc...)
    
    ! At the river mouth
    do i=1,1
        
        covhzx=0.
        ! calculating covariance between Zx and Y
        do j=2,position
            
            dhdhx=exp(sum(dhdhx2(j:position))-dhdhx2(position))  !new
            
            if (j<=i) then
                rz=exp((-sum(L(j:i))+L(i))*2./100000)
            else
                rz=exp(-sum(L(i:j-1))*2./100000)                
            end if
            
            covhzx=covhzx+dhdhx*fx(j)*sqrt(var_Z(i)*var_Z(j))*rz
            
        end do
        
        K(i)=covhzx/(var_H(position)+std_obs*std_obs)
        Z(i)=Z(i)+K(i)*dh/real(m)
        var_Z(i)=var_Z(i)-K(i)*covhzx/real(m)
        acum=fx(i)*K(i)*dh/real(m)
        Y(i)=Y(i)+acum
        H(i)=Y(i)-Z(i)

    end do
    
    
    !For the whole river
    do i=2,nC
        
        covhzx=0.
        ! calculating covariance between Zx and Y
        do j=2,position
            
            dhdhx=exp(sum(dhdhx2(j:position))-dhdhx2(position))  !new
            
            if (j<=i) then
                rz=exp((-sum(L(j:i))+L(i))*2./100000.)
            else
                rz=exp(-sum(L(i:j-1))*2./100000.)                
            end if
            
            covhzx=covhzx+dhdhx*fx(j)*sqrt(var_Z(i)*var_Z(j))*rz
            
        end do
        
        K(i)=covhzx/(var_H(position)+std_obs*std_obs)
        Z(i)=Z(i)+K(i)*dh/real(m)
        var_Z(i)=var_Z(i)-K(i)*covhzx/real(m)
        acum=fx(i)*K(i)*dh/real(m)+acum*exp(dhdhx2(i-1))
        Y(i)=Y(i)+acum
        H(i)=Y(i)-Z(i)
        if(H(i)<=0) then
            H(i)=0.3
            Y(i)=Y(i)+0.3
        end if

    end do
    
    !print*, fx(position),fx2(position)
    !print*, var_H(position), K(position)
    
end do

do i=2,nC
    if(Y(i-1)>=Y(i)) then
        H(i)=H(i)+Y(i-1)-Y(i)
        Y(i)=Y(i-1)
    end if
end do
Y=H+Z

end subroutine
