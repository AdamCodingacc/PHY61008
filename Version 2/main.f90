!Used Sun, Earth, Jupiter
PROGRAM three_body

    !Declare variables
    IMPLICIT NONE
    !x,y,z position, velocity, mass, initial acceln, new acceln for 3 bodies
    DOUBLEPRECISION :: r(1:3, 1:3), v(1:3, 1:3), M(1:3), a_0(1:3, 1:3), a_1(1:3, 1:3)
    DOUBLEPRECISION :: s_sq(1:3, 1:3), dt, G, time, GPE(1:3), absv_sq(1:3), E_0, E_1
    INTEGER :: i, j, n, k

!Initialize variables
G = 1*6.67e-11
n = 3
dt = 1000
M(1) = 1.99e30
M(2) = 1*5.97e24
M(3) = 1*1.898e27

s_sq = 0
time = 0
E_0 = 0
E_1 = 0

!Initialize starting positions
r(:,1) = 0
r(1,2) = 1.496e11 * (1.)
r(2,2) = 1.496e11 * cos(1.)
r(3,2) = 0
r(1,3) = 5.2 * 1.496e11 * sin(1.)
r(2,3) = 5.2 * 1.496e11 * cos(1.)
r(3,3) = 0

!Initialize starting velocities
v(:,1) = 0
v(1,2) = 29780 * cos(1.)
v(2,2) = 29780 * sin(1.)
v(3,2) = 0
v(1,3) = 13060 * cos(1.)
v(2,3) = 13060 * sin(1.)
v(3,3) = 0

!Calculate initial conditions and thus initial acceleration
DO i = 1,n
            DO j = 1,n
            !skip loop if i & j same (no force exerted by body on self)
            IF (i == j) CYCLE

                !Find absolute distance squared between bodies i & j
                s_sq(i,j) = ((r(1,j) - r(1,i))**2 + (r(2,j) - r(2,i))**2 + (r(3,j) - r(3,i))**2)

                a_0(:,i) = G * M(j) * (r(:,j) - r(:,i))*(s_sq(i,j)**-1.5)

                GPE(i) = GPE(i) - G*M(i)*M(j)/sqrt(s_sq(i,j))
            END DO
        END DO

GPE = GPE*0.5
DO i = 1,n
    DO k = 1,3
        absv_sq(i) = absv_sq(i) + v(k,i)**2
    END DO
    E_0 = E_0 + 0.5*M(i)*absv_sq(i) + GPE(i)
END DO

WRITE(6,*) "Initial Conditions"
WRITE(6,*) "Position Vectors xyz (m)"
WRITE(6,*) ""
!Print initial coordinates of all bodies
DO i = 1,n
    WRITE(6,*) r(:,i)
END DO

WRITE(6,*) ""
WRITE(6,*) "Acceleration xyz (ms^-2)"
WRITE(6,*) ""
DO i = 1,n
    WRITE(6,*) a_0(:,i)
END DO

WRITE(6,*) ""
WRITE(6,*) "Initial Energy of System (J)"
WRITE(6,*) ""
WRITE(6,*) E_0


    DO
        DO i = 1,n
            DO j = 1,n
            !skip loop if i & j same (no force exerted by body on self)
            IF (i == j) CYCLE

                !Find absolute distance, s,  squared between bodies i & j
                s_sq(i,j) = ((r(1,j) - r(1,i))**2 + (r(2,j) - r(2,i))**2 + (r(3,j) - r(3,i))**2)

                a_1(:,i) = G * M(j) * (r(:,j) - r(:,i))*(s_sq(i,j)**-1.5)

            END DO
        END DO

        !Find new velocities after time-step
        v = v + 0.5*(a_0 + a_1)*dt
        a_0 = a_1

        !Find new position, r, after time-step
        r = r + v*dt + 0.5*a_0*dt**2

        !updated elapsed time
        time = time + dt
        IF (time > 31536000) EXIT
    END DO

WRITE(6,*) ""
WRITE(6,*) "Final Conditions"
WRITE(6,*) "Position Vectors xyz (m)"
WRITE(6,*) ""
!Print final coordinates of all bodies
DO i = 1,n
    WRITE(6,*) r(:,i)
END DO

WRITE(6,*) ""
WRITE(6,*) "Acceleration xyz (ms^-2)"
WRITE(6,*) ""
DO i = 1,n
    WRITE(6,*) a_0(:,i)
END DO


!Loop through each body to find GPE of each
GPE = 0
absv_sq = 0

DO i = 1,n
        DO j = 1,n
            IF (i == j) CYCLE
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
Write(6,*) ""
write(6,*) "Final Energy of System (J)"
write(6,*) E_1
Write(6,*) "Percentage of Initial Energy Retained"
write(6,*) (E_1/E_0) * 100

END PROGRAM three_body

