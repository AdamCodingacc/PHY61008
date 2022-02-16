!Used Sun, Earth, Jupiter, Mars, Saturn, Uranus, Neptune
PROGRAM solar_sim

    !Declare variables
    IMPLICIT NONE
    !x,y,z position, velocity, mass, initial acceln, new acceln for 7 bodies
    DOUBLEPRECISION :: r(1:3, 1:7, 0:2), M(1:7), v(1:3, 1:7, -6:2),  a(1:3, 1:7, -6:2), a_int(1:3, 1:7, 1:2), v_int(1:3, 1:7, 1:2)
    DOUBLEPRECISION :: COM(1:3), COV(1:3), s_sq(1:7, 1:7), s(1:7, 1:7), GPE(1:7), absv_sq(1:7), dist(1:3, 1:7, 1:7)
    DOUBLEPRECISION :: E_0, E_1, E_check, Mtot, AU, dt, G, time, num_yr, RelErr, Small, errcoeff, KE, pot
    INTEGER :: n, i, j, k, z, stepno, counter !Define indexing integers
    CHARACTER(LEN = 1) :: lg

!Initialize variables
G = 6.67e-11
n = 7
dt = 1.
AU = 1.496e11
RelErr = 5.e-13
Small = 1.e-7
errcoeff = 19./270. !Saves doing the calculation every loop

!Masses in kg
M(1) = 1.99e30
M(2) = 5.97e24
M(3) = 6.4171e23
M(4) = 1.898e27
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
v = 0
r = 0
COM = 0
COV = 0
stepno = 0
counter = 0

!Initialize starting positions
r(:,1,0) = 0 !System centred around Sun
r(1,2,0) = AU * sin(1.) * -1.
r(2,2,0) = AU * cos(1.)
r(3,2,0) = 0
r(1,3,0) = 1.52 * AU * sin(1.) * -1.
r(2,3,0) = 1.52 * AU * cos(1.)
r(3,3,0) = 0
r(1,4,0) = 5.2 * AU * sin(1.) * -1.
r(2,4,0) = 5.2 * AU * cos(1.)
r(3,4,0) = 0
r(1,5,0) = 9.583 * AU * sin(1.) * -1.
r(2,5,0) = 9.583 * AU * cos(1.)
r(3,5,0) = 0
r(1,6,0) = 19.201 * AU * sin(1.) * -1.
r(2,6,0) = 19.201 * AU * cos(1.)
r(3,6,0) = 0
r(1,7,0) = 30.07 * AU * sin(1.) * -1.
r(2,7,0) = 30.07 * AU * cos(1.)
r(3,7,0) = 0

!Initialize starting velocities
!Negative ensures planets begin orbit counterclockwise
v(:,1,0) = 0
v(1,2,0) = 29780. * -1. * cos(1.)
v(2,2,0) = 29780. * -1. * sin(1.)
v(3,2,0) = 0
v(1,3,0) = 24070. * -1. * cos(1.)
v(2,3,0) = 24070. * -1. * sin(1.)
v(3,3,0) = 0
v(1,4,0) = 13060. * -1. * cos(1.)
v(2,4,0) = 13060. * -1. * sin(1.)
v(3,4,0) = 0
v(1,5,0) = 9680. * -1. * cos(1.)
v(2,5,0) = 9680. * -1. * sin(1.)
v(3,5,0) = 0
v(1,6,0) = 6800. * -1. * cos(1.)
v(2,6,0) = 6800. * -1. * sin(1.)
v(3,6,0) = 0
v(1,7,0) = 5430. * -1. * cos(1.)
v(2,7,0) = 5430. * -1. * sin(1.)
v(3,7,0) = 0

!Find centre of mass and velocity
DO i = 1, n
    COM(:) = COM(:) + r(:,i,0)*M(i)
    COV(:) = COV(:) + v(:,i,0)*M(i)
    Mtot = Mtot + M(i)
END DO

!Correct for COM and COV
DO i = 1,n
    r(:,i,0) = r(:,i,0) - COM/Mtot
    v(:,i,0) = v(:,i,0) - COV/Mtot

    DO k = 1,3
        absv_sq(i) = absv_sq(i) + v(k,i,0)**2
    END DO
END DO

!==================================================================================================================

!Calculate initial conditions and thus initial acceleration
DO i = 1,n
            DO j = 1,n
            !skip loop if i & j same (no force exerted by body on self)
            IF (i == j) CYCLE

                !Find absolute distance squared between bodies i & j
                s_sq(i,j) = ((r(1,j,0) - r(1,i,0))**2 + (r(2,j,0) - r(2,i,0))**2 + (r(3,j,0) - r(3,i,0))**2)

                !GPE subtracted thus negative value
                GPE(i) = GPE(i) - G*M(i)*M(j)/sqrt(s_sq(i,j))


                !Need + so each object does not overwrite previous
                a(:,i,0) = a(:,i,0) + (G * M(j) * (r(:,j,0) - r(:,i,0)) * (s_sq(i,j)**-1.5))

            END DO
        END DO

GPE = GPE*0.5

KE = 0
pot = 0
DO i = 1,n
    !GPE added as components negative above
    E_0 = E_0 + 0.5*M(i)*absv_sq(i) + GPE(i)
    KE = KE + 0.5*M(i)*absv_sq(i)
    pot = pot + GPE(i)
END DO

WRITE(6,*) "KE/GPE ratio"
WRITE(6,*) KE / pot


WRITE(6,*) "Initial Conditions"
WRITE(6,*) "Position Vectors xyz (AU)"
WRITE(6,*) ""
!Print initial coordinates of all bodies
DO i = 1,n
    WRITE(6,*) (r(:,i,0)/AU)
END DO

WRITE(6,*) ""
WRITE(6,*) "Velocity xyz (ms^-1)"
WRITE(6,*) ""
!Print initial coordinates of all bodies
DO i = 1,n
    WRITE(6,*) v(:,i,0)
END DO

WRITE(6,*) ""
WRITE(6,*) "Acceleration xyz (ms^-2)"
WRITE(6,*) ""
DO i = 1,n
    WRITE(6,*) a(:,i,0)
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

OPEN(11, file = 'Energy_err.csv')

!================================================================================================================

!Bootstrap
DO k = 0,5

    !Find new position, r, after time-step
    !Comes at the start of the array to get the new position before other variables updated
    r(:,:,0) = r(:,:,0) + v(:,:,0)*dt + 0.5*a(:,:,0)*dt**2

    !Shift previous values down arrays
    DO z = -3,-1
        a(:,:,z) = a(:,:,z+1)
        v(:,:,z) = v(:,:,z+1)
    END DO

    DO i = 1,n
        DO j = 1,n
        !skip loop if i & j same (no force exerted by body on self)
        IF (i == j) CYCLE

            !Find absolute distance squared between bodies i & j
            s_sq(i,j) = ((r(1,j,0) - r(1,i,0))**2 + (r(2,j,0) - r(2,i,0))**2 + (r(3,j,0) - r(3,i,0))**2)


            !Find acceleration on object i in each dimension
            a(:,i,0) = a(:,i,0) + G * M(j) * (r(:,j,0) - r(:,i,0))*(s_sq(i,j)**-1.5)


        END DO
    END DO

    !Find new velocities after time-step
    v(:,:,0) = v(:,:,0) + 0.5*(a(:,:,-1) + a(:,:,0))*dt

    time = time + dt
END DO

!==================================================================================================================

!ABM Method
DO

    !Print current Earth distance
    !WRITE(6,*) (r(:,2,0) /AU)

    !Clear new acceleration to not factor into calculation
    a(:,:,1) = 0

    !Predict position and velocity using ABM predictor
    r(:,:,1) = r(:,:,0) + ((dt/24.) * (-9.*v(:,:,-3) + 37.*v(:,:,-2) - 59.*v(:,:,-1) + 55.*v(:,:,0)))
    v(:,:,1) = v(:,:,0) + ((dt/24.) * (-9.*a(:,:,-3) + 37.*a(:,:,-2) - 59.*a(:,:,-1) + 55.*a(:,:,0)))

    !Calculate acceleration for predicted position
    DO i = 1,n
        DO j = 1,n
        IF (i == j) CYCLE
            dist(:, i, j) = r(:,j,1) - r(:,i,1)

            s_sq(i,j) = ((r(1,j,1) - r(1,i,1))**2 + (r(2,j,1) - r(2,i,1))**2 + (r(3,j,1) - r(3,i,1))**2)

            a(:,i,1) = a(:,i,1) + G * M(j) * (r(:,j,1) - r(:,i,1))*(s_sq(i,j)**-1.5)
        END DO
    END DO



    !ABM corrector
    r(:,:,2) = r(:,:,0) + ((dt/24) * (v(:,:,-2) - 5.*v(:,:,-1) + 19.*v(:,:,0) + 9.*v(:,:,1)))
    v(:,:,2) = v(:,:,0) + ((dt/24) * (a(:,:,-2) - 5.*a(:,:,-1) + 19.*a(:,:,0) + 9.*a(:,:,1)))


    !Recalculate acceleration for corrected position
    a(:,:,2) = 0
    DO i = 1,n
        DO j = 1,n
        IF (i == j) CYCLE
            dist(:, i, j) = r(:,j,2) - r(:,i,2)

            s_sq(i,j) = ((r(1,j,2) - r(1,i,2))**2 + (r(2,j,2) - r(2,i,2))**2 + (r(3,j,2) - r(3,i,2))**2)

            a(:,i,2) = a(:,i,2) + G * M(j) * (r(:,j,2) - r(:,i,2))*(s_sq(i,j)**-1.5)
        END DO
    END DO

    !Shift previous values down arrays
    DO z = -6,0
        a(:,:,z) = a(:,:,z+1)
        v(:,:,z) = v(:,:,z+1)
    END DO
    !Set corrected values as current position, velocity and acceleration
    r(:,:,0) = r(:,:,2)
    v(:,:,0) = v(:,:,2)
    a(:,:,0) = a(:,:,2)



GPE = 0
absv_sq = 0
DO i = 1,n
            DO j = 1,n
            !skip loop if i & j same (no force exerted by body on self)
            IF (i == j) CYCLE
                GPE(i) = GPE(i) - G*M(i)*M(j)/sqrt(s_sq(i,j))
            END DO
        END DO
GPE = GPE*0.5

DO i = 1,n
    E_check = E_check + 0.5*M(i)*absv_sq(i) + GPE(i)
END DO
!WRITE(6,*) "Current Energy error"
!WRITE(6,*) E_0 - E_check

!Write every 100th step energy error to a log file
IF ((MOD(stepno, 250) == 0)) THEN
    WRITE(11,*) time, ',', E_0 - E_check, ',', E_check / E_0
END IF





    !Ensure enough data points are stored
    IF (counter > 6) THEN
        !Check if error from predictor-corrector is small enough
        IF (errcoeff * MAXVAL(ABS(r(:,:,2) - r(:,:,1)) / (ABS(r(:,:,2)) + Small)) < (RelErr * 0.01)) THEN

            !Double timestep
            dt = dt * 2.

            !Omit alternate points to adjust timestep
            DO k = 1,3
                a(:,:,-k) = a(:,:,-2*k)
                v(:,:,-k) = v(:,:,-2*k)
            END DO
        END IF

        !Check if error from predictor-corrector is too large
        IF (errcoeff * MAXVAL(ABS(r(:,:,2) - r(:,:,1)) / (ABS(r(:,:,2)) + Small)) > RelErr) THEN

            !halve timestep
            dt = dt * 0.5

            !Interpolate previous data points
            a_int(:,:,1) = (-5. * a(:,:,-4) + 28. * a(:,:,-3) - 70. * a(:,:,-2) + 140. * a(:,:,-1) + 35. * a(:,:,0)) / 128.
            a_int(:,:,2) = (3. * a(:,:,-4) - 20. * a(:,:,-3) + 90. * a(:,:,-2) + 60. * a(:,:,-1) -5. * a(:,:,0)) / 128.

            v_int(:,:,1) = (-5. * v(:,:,-4) + 28. * v(:,:,-3) - 70. * v(:,:,-2) + 140. * v(:,:,-1) + 35. * v(:,:,0)) / 128.
            v_int(:,:,2) = (3. * v(:,:,-4) - 20. * v(:,:,-3) + 90. * v(:,:,-2) + 60. * v(:,:,-1) -5. * v(:,:,0)) / 128.

            !Create new mesh
            !-5 and -6 unimportant, not used in calculation & overwritten before next time-step change
            a(:,:,-4) = a(:,:,-2)
            a(:,:,-3) = a_int(:,:,2)
            a(:,:,-2) = a(:,:,-1)
            a(:,:,-1) = a_int(:,:,1)

            v(:,:,-4) = v(:,:,-2)
            v(:,:,-3) = v_int(:,:,2)
            v(:,:,-2) = v(:,:,-1)
            v(:,:,-1) = v_int(:,:,1)


        END IF

        !Reset counter to let old data fill up
        counter = 0
    END IF



    !Write every 500th step to log files
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
    counter = counter + 1
    IF (time > (num_yr * 31536000.)) EXIT
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
    WRITE(6,*) (r(:,i,0)/AU)
END DO

WRITE(6,*) ""
WRITE(6,*) "Velocity xyz (ms^-1)"
WRITE(6,*) ""
DO i = 1,n
    DO k = -3,1
        WRITE(6,*) v(:,i,k)
    END DO
END DO

WRITE(6,*) ""
WRITE(6,*) "Acceleration xyz (ms^-2)"
WRITE(6,*) ""
DO i = 1,n
    DO k = -3,1
        WRITE(6,*) a(:,i,k)
    END DO
END DO

WRITE(6,*) dt

!Loop through each body to find GPE of each
GPE = 0
absv_sq = 0

DO i = 1,n
        DO j = 1,n
            IF (i == j) CYCLE
                s_sq(i,j) = ((r(1,j,0) - r(1,i,0))**2 + (r(2,j,0) - r(2,i,0))**2 + (r(3,j,0) - r(3,i,0))**2)

                !Sum calculated GPE
                GPE(i) = GPE(i) - G*M(i)*M(j)/sqrt(s_sq(i,j))
        END DO

        !Sum kinetic energy in each dimension (k)
        DO k = 1,3
            absv_sq(i) = absv_sq(i) + v(k,i,0)**2
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

WRITE(6,*) "Absolute distance from the Sun (AU)"
DO k = 1,7
    WRITE(6,*) sqrt(s_sq(1,k)) /AU
END DO
END PROGRAM solar_sim
