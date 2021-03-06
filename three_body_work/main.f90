
!Used Sun, Earth, Jupiter
PROGRAM three_body

    !Declare variables
    IMPLICIT NONE
    !x,y,z position, velocity, mass, initial acceln, new acceln for 3 bodies
    DOUBLEPRECISION :: r(1:3, 1:7), v(1:3, 1:7), M(1:7), a_0(1:3, 1:7), a_1(1:3, 1:7)
    DOUBLEPRECISION :: COM(1:3), COV(1:3), s_sq(1:7, 1:7), dt, G, time, GPE(1:7), absv_sq(1:7), E_0, E_1, Mtot, AU
    INTEGER :: n, i, j, k !Define indexing integers

!Initialize variables
G = 6.67e-11
n = 6
dt = 100
AU = 1.496e11
M(1) = 1.99e30
M(2) = 5.97e24
M(3) = 1.898e27
M(4) = 6.4171e23
M(5) = 5.683e26
M(6) = 8.681e25
M(7) = 1.024e26

!Empty array variables
s_sq = 0
time = 0
GPE = 0
absv_sq = 0
E_0 = 0
E_1 = 0
a_0 = 0
a_1 = 0
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
        DO k = 1,3
                COM(k) = COM(k) + r(k,i)*M(i)
                COV(k) = Cov(k) + v(k,i)*M(i)
        END DO
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

                DO k = 1,3 !Loop over each dimension
                    !Need + so each object does not overwrite previous
                    a_0(k,i) = a_0(k,i) + (G * M(j) * (r(k,j) - r(k,i)) * (s_sq(i,j)**-1.5))
                END DO
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

                DO k = 1,3
                    !Find acceleration on object i in each dimension
                    a_1(k,i) = a_1(k,i) + G * M(j) * (r(k,j) - r(k,i))*(s_sq(i,j)**-1.5)
                END DO

            END DO
        END DO

        !Find new velocities after time-step
        v = v + 0.5*(a_0 + a_1)*dt
        a_0 = a_1
        a_1 = 0

        !Find new position, r, after time-step
        r = r + v*dt + 0.5*a_0*dt**2

        !updated elapsed time
        time = time + dt
        IF (time > 31536000) EXIT
    END DO

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
WRITE(6,*) ""
WRITE(6,*) "Final Energy of System (J)"
WRITE(6,*) E_1
WRITE(6,*) "Percentage of Initial Energy Retained"
WRITE(6,*) (E_1/E_0) * 100

END PROGRAM three_body
