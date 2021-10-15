!Used Sun, Earth, Jupiter, Mars, Saturn, Uranus, Neptune
PROGRAM solar_sim

    !Declare variables
    IMPLICIT NONE
    !x,y,z position, velocity, mass, initial acceln, new acceln for 7 bodies
    DOUBLEPRECISION :: r(1:3, 1:7), M(1:7), v(1:3, 1:7),  a(1:3, 1:7, -3:1)
    DOUBLEPRECISION :: COM(1:3), COV(1:3), s_sq(1:7, 1:7), s(1:7, 1:7), GPE(1:7), absv_sq(1:7), dist(1:3, 1:7, 1:7)
    DOUBLEPRECISION :: E_0, E_1, Mtot, AU, dt, G, time
    INTEGER :: num_yr, n, i, j, k, z, b, stepno !Define indexing integers
    CHARACTER(LEN = 1) :: lg

!Initialize variables
G = 6.67e-11
n = 7
dt = 1
AU = 1.496e11
M(1) = 1.99e30
M(2) = 5.97e24
M(3) = 1.898e27
M(4) = 6.4171e23
M(5) = 5.683e26
M(6) = 8.681e25
M(7) = 1.024e26

!Empty array variables
num_yr = 0
s_sq = 0
time = 0
GPE = 0
absv_sq = 0
E_0 = 0
E_1 = 0
a = 0
COM = 0
COV = 0

!Initialize starting positions
r(:,1) = 0 !System centred around Sun
r(1,2) = AU * sin(1.) * -1.
r(2,2) = AU * cos(1.)
r(3,2) = 0
r(1,3) = 5.2 * AU * sin(1.) * -1.
r(2,3) = 5.2 * AU * cos(1.)
r(3,3) = 0
r(1,4) = 1.52 * AU * sin(1.) * -1.
r(2,4) = 1.52 * AU * cos(1.)
r(3,4) = 0
r(1,5) = 9.583 * AU * sin(1.) * -1.
r(2,5) = 9.583 * AU * cos(1.)
r(3,5) = 0
r(1,6) = 19.201 * AU * sin(1.) * -1.
r(2,6) = 19.201 * AU * cos(1.)
r(3,6) = 0
r(1,7) = 30.07 * AU * sin(1.) * -1.
r(2,7) = 30.07 * AU * cos(1.)
r(3,7) = 0

!Initialize starting velocities
v(:,1) = 0
v(1,2) = 29780 * -1. * cos(1.)
v(2,2) = 29780 * -1. * sin(1.)
v(3,2) = 0
v(1,3) = 13060 * -1. * cos(1.)
v(2,3) = 13060 * -1. * sin(1.)
v(3,3) = 0
v(1,4) = 24070 * -1. * cos(1.)
v(2,4) = 24070 * -1. * sin(1.)
v(3,4) = 0
v(1,5) = 9680 * -1. * cos(1.)
v(2,5) = 9680 * -1. * sin(1.)
v(3,5) = 0
v(1,6) = 6800 * -1. * cos(1.)
v(2,6) = 6800 * -1. * sin(1.)
v(3,6) = 0
v(1,7) = 5430 * -1. * cos(1.)
v(2,7) = 5430 * -1. * sin(1.)
v(3,7) = 0

DO i = 1, n-1
    COM(:) = COM(:) + r(:,i)*M(i)
    COV(:) = COV(:) + v(:,i)*M(i)
    Mtot = Mtot + M(i)
END DO

!Correct for COM and COV
DO i = 1,n
    r(:,i) = r(:,i) - COM/Mtot
    v(:,i) = v(:,i) - COV/Mtot

    DO k = 1,3
        absv_sq(i) = absv_sq(i) + v(k,i)**2
    END DO
END DO

!Calculate initial conditions and thus initial acceleration
DO i = 1,n
            DO j = 1,n
            !skip loop if i & j same (no force exerted by body on self)
            IF (i == j) CYCLE

                !Find absolute distance squared between bodies i & j
                s_sq(i,j) = ((r(1,j) - r(1,i))**2 + (r(2,j) - r(2,i))**2 + (r(3,j) - r(3,i))**2)

                GPE(i) = GPE(i) - G*M(i)*M(j)/sqrt(s_sq(i,j))


                !Need + so each object does not overwrite previous
                a(:,i,-3) = a(:,i,-3) + (G * M(j) * (r(:,j) - r(:,i)) * (s_sq(i,j)**-1.5))

            END DO
        END DO

GPE = GPE*0.5


DO i = 1,n
    E_0 = E_0 + 0.5*M(i)*absv_sq(i) + GPE(i)
END DO

WRITE(6,*) "Initial Conditions"
WRITE(6,*) "Position Vectors xyz (AU)"
WRITE(6,*) ""
!Print initial coordinates of all bodies
DO i = 1,n
    WRITE(6,*) (r(:,i)/AU)
END DO

WRITE(6,*) ""
WRITE(6,*) "Velocity xyz (ms^-1)"
WRITE(6,*) ""
!Print initial coordinates of all bodies
DO i = 1,n
    WRITE(6,*) v(:,i)
END DO

WRITE(6,*) ""
WRITE(6,*) "Acceleration xyz (ms^-2)"
WRITE(6,*) ""
DO i = 1,n
    WRITE(6,*) a(:,i,-3)
END DO

WRITE(6,*) ""
WRITE(6,*) "Initial Energy of System (J)"
WRITE(6,*) ""
WRITE(6,*) E_0

WRITE(6,*) ""
WRITE(6,*) "Would you like to save this run? (y/n)"
READ *, lg !Store user input to decide whether to log data
WRITE(6,*) ""

WRITE(6,*) "How many years would you like to simulate for?"
READ *, num_yr !Store user input to decide whether to log data
WRITE(6,*) ""

!Open Log Files to write data to if requested
IF (lg == 'y') THEN
    OPEN(2,file = 'Sun_Motion.csv')
    OPEN(3,file = 'Earth_Motion.csv')
    OPEN(4,file = 'Jupiter_Motion.csv')
    OPEN(7,file = 'Mars_Motion.csv')
    OPEN(8,file = 'Saturn_Motion.csv')
    OPEN(9,file = 'Uranus_Motion.csv')
    OPEN(10,file = 'Neptune_Motion.csv')
END IF

!Bootstrap
DO k = 0,3
    a(:,:,-3) = a(:,:,-2)
    a(:,:,-2) = a(:,:,-1)
    a(:,:,-1) = a(:,:,0)
    a(:,:,0) = 0

    s_sq = 0
    dist = 0



    DO i = 1,n
        DO j = 1,n
        !skip loop if i & j same (no force exerted by body on self)
        IF (i == j) CYCLE

            !Find absolute distance squared between bodies i & j
            s_sq(i,j) = ((r(1,j) - r(1,i))**2 + (r(2,j) - r(2,i))**2 + (r(3,j) - r(3,i))**2)


            !Find acceleration on object i in each dimension
            a(:,i,0) = a(:,i,0) + G * M(j) * (r(:,j) - r(:,i))*(s_sq(i,j)**-1.5)


        END DO
    END DO
    !Find new position, r, after time-step
    r = r + v*dt + 0.5*a(:,:,-1)*dt**2

    !Find new velocities after time-step
    v = v + 0.5*(a(:,:,-1) + a(:,:,0))*dt

    time = time + dt
END DO

WRITE(6,*) ""
WRITE(6,*) "Bootstrap Conditions"
WRITE(6,*) "Position Vectors xyz (AU)"
WRITE(6,*) ""
!Print final coordinates of all bodies
DO i = 1,n
    WRITE(6,*) (r(:,i)/AU)
END DO

WRITE(6,*) ""
WRITE(6,*) "Velocity xyz (ms^-1)"
WRITE(6,*) ""
DO i = 1,n
    WRITE(6,*) v(:,i)
END DO

WRITE(6,*) ""
WRITE(6,*) "Acceleration xyz (ms^-2)"
WRITE(6,*) ""
DO i = 1,n
    DO k = -3,0
        WRITE(6,*) a(:,i,k)
    END DO
END DO

dt = 1000
DO
    DO i = 1,n
        DO j = 1,n
        !skip loop if i & j same (no force exerted by body on self)
        IF (i == j) CYCLE
            dist(:, i, j) = r(:, j) - r(:, i)
        END DO
    END DO

    a(:,:,1) = a(:,:,0) + (dt/24)*(-9*a(:,:,-3) + 37*a(:,:,-2) - 59*a(:,:,-1) + 55*a(:,:,0))


    !Find new velocities after time-step
    v = v + 0.5*(a(:,:,0) + a(:,:,1))*dt

    !Find new position, r, after time-step
    r = r + v*dt + 0.5*a(:,:,1)*dt**2

    a(:,:,-3) = a(:,:,-2)
    a(:,:,-2) = a(:,:,-1)
    a(:,:,-1) = a(:,:,0)
    a(:,:,0) = a(:,:,1)

    !Write every 250th step to log files
    IF ((MOD(stepno, 500) == 0) .and. (lg == 'y')) THEN
        WRITE(2,*) time, ',', dist(1,1,1), ',', dist(2,1,1), ',', dist(3,1,1)
        WRITE(3,*) time, ',', dist(1,1,2), ',', dist(2,1,2), ',', dist(3,1,2)
        WRITE(4,*) time, ',', dist(1,1,3), ',', dist(2,1,3), ',', dist(3,1,3)
        WRITE(7,*) time, ',', dist(1,1,4), ',', dist(2,1,4), ',', dist(3,1,4)
        WRITE(8,*) time, ',', dist(1,1,5), ',', dist(2,1,5), ',', dist(3,1,5)
        WRITE(9,*) time, ',', dist(1,1,6), ',', dist(2,1,6), ',', dist(3,1,6)
        WRITE(10,*) time, ',', dist(1,1,7), ',', dist(2,1,7), ',', dist(3,1,7)
    END IF

    !update elapsed time
    time = time + dt
    stepno = stepno + 1
    IF (time > (num_yr * 31536000)) EXIT
END DO


IF (lg == 'y') THEN
    CLOSE(2)
    CLOSE(3)
    CLOSE(4)
    CLOSE(7)
    CLOSE(8)
    CLOSE(9)
    CLOSE(10)
END IF


WRITE(6,*) ""
WRITE(6,*) "Final Conditions"
WRITE(6,*) "Position Vectors xyz (AU)"
WRITE(6,*) ""
!Print final coordinates of all bodies
DO i = 1,n
    WRITE(6,*) (r(:,i)/AU)
END DO

WRITE(6,*) ""
WRITE(6,*) "Acceleration xyz (ms^-2)"
WRITE(6,*) ""
DO i = 1,n
    WRITE(6,*) a(:,i,1)
END DO


!Loop through each body to find GPE of each
GPE = 0
absv_sq = 0

DO i = 1,n
        DO j = 1,n
            IF (i == j) CYCLE
                s_sq(i,j) = ((r(1,j) - r(1,i))**2 + (r(2,j) - r(2,i))**2 + (r(3,j) - r(3,i))**2)

                !Sum calculated GPE
                GPE(i) = GPE(i) - G*M(i)*M(j)/sqrt(s_sq(i,j))
        END DO

        !Sum kinetic energy in each dimension (k)
        DO k = 1,3
            absv_sq(i) = absv_sq(i) + v(k,i)**2
        END DO
END DO

GPE = GPE*0.5

!Sum energies to get final energy
DO i = 1,n
    E_1 = E_1 + 0.5*M(i)*absv_sq(i) + GPE(i)
END DO

!Report final energy and accuracy test
WRITE(6,*) ""
WRITE(6,*) "Final Energy of System (J)"
WRITE(6,*) E_1
WRITE(6,*) "Percentage of Initial Energy Retained"
WRITE(6,*) (E_1/E_0) * 100

END PROGRAM solar_sim
