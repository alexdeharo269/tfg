program kawasakiv3fortran
    real*8, dimension (:), allocatable :: s, hglauber, hvecchem
    real*8, dimension (:,:), allocatable :: hkawasaki, epsilon, epsilonelec
    integer, dimension (:,:), allocatable :: indiceschem, raster
    real*8 JJ, T, dranu
    integer :: count_attempt, count_accept


! semillas buenas
! 372321
    iseed=32321
    call dranini(iseed)
!!! ficheros
 open (unit = 12, file="./rednuevo.dat", status="unknown")
 open (unit = 10, file="./filenuevo.dat", status="unknown")
 open (unit = 11, file="./fileraster.dat", status="unknown")
!!! parametros
! montecarlo steps
! pp=1 difusion pp=0 reaccion
! para N=400 me salen quimeras para nelec=5, nchem=10 T=0.5, JJ=0.4 y pp=0.3
! 1) si aumentamos nchem debemos favorecer la dinamica de difusion para que aparezcan quimeras
!     por ejemplo nelec 5, nchem 20 pp=0.9
! 2) Tambien ocurre al reves es decir al bajar nechm hay que bajar tb pp para ver quimeras
!    lo que ocurre es que en este caso algunas veces son metaestables y desaparecen despues de un cierto tiempo
!    por ejemplo nelec=10, nchem=10 y pp=0.1, tambien salen muy bien para nelec=10, nchem=8 y p=0.01
! 3) para nchem pequeño puedo hacer que vuelvan a aparecer haciendo mas pequeño nelec para pp=0
!     por ejempl nelec=7, nchem=7 t pp=0. Tambien salen para pp=0.001

    nmax=300
    N=200
    nelec=50
    nchem=0
    nnchem=2*nchem
    T=0.5
    JJ=0.00000
    pp=1.0
    niter=N*nmax

    write(*,*) niter

    allocate(s(N))
    allocate(hglauber(N))
    allocate(hvecchem(N))
    allocate(hkawasaki(N,N))
    allocate(epsilon(N,N))
    allocate(epsilonelec(N,N))
    allocate(indiceschem(N,2*nchem))
    allocate(raster(nmax,N))

! generamos la topologia no local

    do i = 1, N
        do j = 1, nelec
            ivec1 = i + j
            ivec2 = i - j

            ! Apply periodic boundary conditions
            if (ivec1 > N) then
                ivec1 = ivec1 - N
            end if
            if (ivec2 < 1) then
                ivec2 = ivec2 + N
            end if

            ! Assign values to epsilonelec matrix
            epsilonelec(i, ivec1) = 1
            epsilonelec(ivec1, i) = 1
            epsilonelec(i, ivec2) = 1
            epsilonelec(ivec2, i) = 1
        end do
        epsilonelec(i, i) = 0
    end do



    do i = 1, N
        l=0
        do j = 1, nchem
            ivec3 = i + nelec + j
            ivec4 = i - nelec - j

            ! Apply periodic boundary conditions
            if (ivec3 > N) then
                ivec3 = ivec3 - N
            end if
            if (ivec4 < 1) then
                ivec4 = ivec4 + N
            end if
            l=l+1

            indiceschem(i,l)=ivec3
            l=l+1
            indiceschem(i,l)=ivec4

            ! Assign values to epsilon matrix
            epsilon(i, ivec3) = 1
            epsilon(ivec3, i) = 1
            epsilon(i, ivec4) = 1
            epsilon(ivec4, i) = 1
        end do
        
        epsilon(i, i) = 0
        
    end do

! escribimos la red en un fichero
  do i=1,N
     do j=1,N
        
        write(12,*) i , j, epsilonelec(i,j), epsilon(i,j)
    
     end do
  end do



!! condicion inciales de espines

do i=1,N
    s(i)=1
    xr=dranu()

    if (xr>0.75) then
        s(i)=-1
    endif
enddo

!!!!! campos locales iniciales
!!! Glauber

do i=1,N
    hvecchem(i)=0.
    do l=1,nnchem
        is=indiceschem(i,l)
        hvecchem(i)=hvecchem(i)+1.*s(is)
    enddo
    hglauber(i)=2.*JJ*s(i)*hvecchem(i)
enddo



!!! kawasaki

do i=1,N    
    do j=1,N
        sumahkawa=(hvecchem(i)-hvecchem(j)) 
        sumahkawa=sumahkawa-epsilon(i,j)*s(j)+epsilon(j,i)*s(i)
        hkawasaki(i,j)=JJ*(s(i)-s(j))*sumahkawa
        !write(*,*) 'kawa', hkawasaki(i,j)
    enddo
enddo



!!! Simulacion
npinta=N
write(*,*) 'Simulation started...'
flush(6)

count_attempt = 0
count_accept = 0

do it=1,niter


    iran=1+int(dranu()*N)
    jran=1+int(dranu()*N)

    !write(*,*) iran, jran
    isigno=s(iran)*s(jran)  

    xr1=dranu()
    ss=1
    delta=hglauber(iran)
    !prob=min(exp(-delta/T),1.)
    prob=1
    if (delta>0) then
       prob=exp(-delta/T)
    endif  

    if (xr1.lt.pp) then
        count_attempt = count_attempt + 1

        ss=0
        delta=hkawasaki(iran,jran)
        if(delta.ne.0) then
            write(*,*) delta !delta siempre es 0 ahoramismo
        endif
        prob=epsilonelec(iran,jran)*(min(exp(-delta/T),1.))*(1-isigno)*0.5
        !prob=0
        !if (delta<0) then           
        !    prob=1*epsilonelec(iran,jran)*(1-isigno)*0.5
            !write(*,*) delta, prob, isigno
        !endif
    endif

    xr2=dranu()
    
    
    if (prob.gt.xr2) then
        count_accept = count_accept + 1
        iss=s(iran)
        iss2=s(jran)
        s(iran)=iss2*(1-ss)-iss*ss
        s(jran)=ss*s(jran)+(1-ss)*iss
        
    endif




!! nuevos campos locales

    !!! Glauber

    do i=1,N
        hvecchem(i)=0.
        do l=1,nnchem
            is=indiceschem(i,l)
            hvecchem(i)=hvecchem(i)+1.*s(is)
        enddo
        hglauber(i)=2.*JJ*s(i)*hvecchem(i)
    enddo

    !!! kawasaki

    do i=1,N    
        do j=1,N
            sumahkawa=0.
            sumahkawa=hvecchem(i)-hvecchem(j) 
            sumahkawa=sumahkawa-epsilon(i,j)*s(j)+epsilon(j,i)*s(i)
            hkawasaki(i,j)=JJ*(s(i)-s(j))*sumahkawa
        enddo
    enddo
    


! pintamos los datos

    if (it.gt.npinta) then
        itt=it/N
       !write(*,*) itt
        do j=1,N
            
            write(11,*) int(s(j))
            write(10, '(I1, 1X)', advance='no') int((s(j) + 1) / 2)         
        enddo
        write(10, *)  ! Newline
        npinta=npinta+N
        flush(10)
    endif



enddo

write(*,*) 'Final swap acceptance percentage: ', 100.0 * real(count_accept) / real(count_attempt), '%'
flush(6)


end program kawasakiv3fortran
