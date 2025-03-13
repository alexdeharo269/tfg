program kawasakiv3fortran
!MODIFICACIONES:
! He hecho el cálculo de campokawa en el bucle principal, en lugar de hacerlo en un bucle aparte, ahorrando N^2 operaciones
! La dinámica de Glauber ahora comprueba que están conectados según la matrix epsilon
! En el cálculo de campokawa, se usaba epsilon en lugar de epsilonelec. He solucionado el error.
! He cambiado la forma en que se imprime rednuevo y filenuevo, para que sea más fácil de leer.
! Ahora se puede ejecutar en gnuplot como el mío: plot 'filenuevo.dat' matrix with image
! Imprimo en pantalla el porcentaje de aceptación de swaps al final de la simulación
! Inicializo el anillo con un 75% de espines arriba y un 25% de espines abajo
! Por lo demás solo he tocado los parámetros

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

    nmax=1000
    N=1000
    nelec=20
    nchem=10
    nnchem=2*nchem
    T=0.5
    JJ=1.00000
    pp=0.9
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
        
        write(12,'(I1, 1X)', advance='no')  int(epsilon(i,j))
    
     end do
        write(12, *)  ! Newline
  end do



!! condicion inciales de espines

do i=1,N
    s(i)=1
    xr=dranu()

    if (xr>0.5) then
        s(i)=-1
    endif
enddo

!!!!! campos locales iniciales
!!! Glauber







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
    delta=0.
    do j=1,N
        delta=delta+1.*s(j)*epsilon(i,j)
    enddo
    delta=2.*JJ*s(i)*delta
    !prob=min(exp(-delta/T),1.)
    prob=1
    if (delta>0) then
    !me he quedado por aquí
       prob=exp(-delta/T)
    endif  
    !!! aceptamos o no según estén conectados
    prob=prob*epsilon(iran,jran)

    if (xr1.lt.pp) then
        !write(*,*) "error"
        ss=0
        sumahkawa=0.
        do k=1,N
            sumahkawa=sumahkawa-epsilonelec(jran,k)*s(k)+epsilonelec(iran,k)*s(k)
        enddo 
        sumahkawa=sumahkawa-epsilonelec(iran,jran)*s(jran)+epsilonelec(jran,iran)*s(iran)
        delta=JJ*(s(iran)-s(jran))*sumahkawa

        prob=epsilonelec(iran,jran)*(min(exp(-delta/T),1.))*(1-isigno)*0.5
        
    endif

    xr2=dranu()
    if (prob.gt.xr2) then
        !count_accept = count_accept + 1
        iss=s(iran)
        iss2=s(jran)
        s(iran)=iss2*(1-ss)-iss*ss
        s(jran)=ss*s(jran)+(1-ss)*iss
        !if(ss==0) then 
        !aux = s(jran)
        !s(jran)=s(iran)
        !s(iran)=aux
        !endif
        !if(ss==1) then 
        !s(iran)=-s(iran)
        !endif


    endif
    !if (prob.gt.xr2) then 
    !    s(iran)=-s(iran)
    !endif




!! nuevos campos locales

    !!! Glauber

    !do i=1,N
    !    hvecchem(i)=0.
    !    !do l=1,nnchem
    !    !    is=indiceschem(i,l)
    !    !    hvecchem(i)=hvecchem(i)+1.*s(is)
    !    !aenddo
    !    do j=1,N
    !        hvecchem(i)=hvecchem(i)+1.*s(j)*epsilon(i,j)
    !    enddo

    !    hglauber(i)=2.*JJ*s(i)*hvecchem(i)
    !enddo




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
