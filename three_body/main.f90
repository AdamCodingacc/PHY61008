!Used Sun, Earth, Jupiter
PROGRAM three_body

    !Declare variables
    IMPLICIT NONE
    !x,y,z position, velocity, mass, initial acceln, new acceln for 3 bodies
    DOUBLEPRECISION :: r(1:3, 1:3), v(1:3, 1:3), M(1:3), a_0(1:3, 1:3), a_1(1:3, 1:3)
    DOUBLEPRECISION :: ssq(1:3, 1:3), dt, G, time
    INTEGER :: i, j, n

!Initialize variables
G = 1*6.67e-11
n = 3
ssq = 0
time = 0
dt = 100
M(1) = 1.99e30
M(2) = 1*5.97e24
M(3) = 1*1.898e27

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
                ssq(i,j) = ((r(1,j) - r(1,i))**2 + (r(2,j) - r(2,i))**2 + (r(3,j) - r(3,i))**2)

                a_0(:,i) = G * M(j) * (r(:,j) - r(:,i))*(ssq(i,j)**-1.5)
                !a_1(i,2) = G * (r(j,2) - r(i,2))**-1.5 * M(j)
                !a_1(i,3) = G * (r(j,3) - r(i,3))**-1.5 * M(j)
            END DO
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

    DO
        DO i = 1,n
            DO j = 1,n
            !skip loop if i & j same (no force exerted by body on self)
            IF (i == j) CYCLE

                !Find absolute distance squared between bodies i & j
                ssq(i,j) = ((r(1,j) - r(1,i))**2 + (r(2,j) - r(2,i))**2 + (r(3,j) - r(3,i))**2)

                a_1(:,i) = G * M(j) * (r(:,j) - r(:,i))*(ssq(i,j)**-1.5)
                !a_1(i,2) = G * (r(j,2) - r(i,2))**-1.5 * M(j)
                !a_1(i,3) = G * (r(j,3) - r(i,3))**-1.5 * M(j)
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
END PROGRAM three_body

