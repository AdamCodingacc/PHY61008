!Used Sun, Earth, Jupiter, Mars, Saturn, Uranus, Neptune
PROGRAM solar_sim

    !Declare variables
    IMPLICIT NONE
    DOUBLEPRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: r, v, a, a_int, v_int, dist
    DOUBLEPRECISION, ALLOCATABLE, DIMENSION(:,:) :: s_sq, rtemp
    DOUBLEPRECISION, ALLOCATABLE, DIMENSION(:) :: absv_sq, M, GPE, dt, deltat, counter, phi
    DOUBLEPRECISION :: COM(1:3), COV(1:3)
    DOUBLEPRECISION :: E_0, E_1, Mtot, AU, dtmin, G, time, num_yr, RelErr, Small, errcoeff, KE, pot, exittime, pi
    INTEGER :: n, i, j, k, z, stepno !Define indexing integers
    CHARACTER(LEN = 1) :: lg

!Allocate array space for the number of bodies
n = 8
ALLOCATE(r(1:3, 1:n, 0:2))
ALLOCATE(v(1:3, 1:n, -6:2))
ALLOCATE(a(1:3, 1:n, -6:2))
ALLOCATE(a_int(1:3, 1:n, 1:2))
ALLOCATE(v_int(1:3, 1:n, 1:2))
ALLOCATE(dist(1:3, 1:n, 1:n))
ALLOCATE(s_sq(1:n, 1:n))
ALLOCATE(rtemp(1:3, 1:n))
ALLOCATE(absv_sq(1:n))
ALLOCATE(M(1:n))
ALLOCATE(GPE(1:n))
ALLOCATE(dt(1:n))
ALLOCATE(deltat(1:n))
ALLOCATE(counter(1:n))
ALLOCATE(phi(1:n))

!Initialize variables
G = 6.67e-11
dtmin = 1.
AU = 1.496e11
RelErr = 5.e-13
Small = 1.e-7
errcoeff = 19./270. !Saves doing the calculation every loop
pi = 3.14159

!Initialize random angles
!Call init_random_seed()
!Call random_number(phi)
!phi = phi * 2* pi !angle in radians

!Masses in kg
M(1) = 1.99e30
M(2) = 5.97e24
M(3) = 6.4171e23
M(4) = 1.898e27
M(5) = 5.683e26
M(6) = 8.681e25
M(7) = 1.024e26
M(8) = 0.5 * 1.99e30

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
dt = 0
deltat = 0
rtemp = 0
phi = 1

!Initialize starting positions
r(:,1,0) = 0 !System centred around Sun
r(1,2,0) = AU * sin(phi(1)) * (-1.)
r(2,2,0) = AU * cos(phi(1))
r(3,2,0) = 0
r(1,3,0) = 1.52 * AU * sin(phi(2)) * (-1.)
r(2,3,0) = 1.52 * AU * cos(phi(2))
r(3,3,0) = 0
r(1,4,0) = 5.2 * AU * sin(phi(3)) * (-1.)
r(2,4,0) = 5.2 * AU * cos(phi(3))
r(3,4,0) = 0
r(1,5,0) = 9.583 * AU * sin(phi(4)) * (-1.)
r(2,5,0) = 9.583 * AU * cos(phi(4))
r(3,5,0) = 0
r(1,6,0) = 19.201 * AU * sin(phi(5)) * (-1.)
r(2,6,0) = 19.201 * AU * cos(phi(5))
r(3,6,0) = 0
r(1,7,0) = 30.07 * AU * sin(phi(6)) * (-1.)
r(2,7,0) = 30.07 * AU * cos(phi(6))
r(3,7,0) = 0
r(1,8,0) = 500 * AU * sin(phi(7)) * (-1.)
r(2,8,0) = 500 * AU * cos(phi(7))
r(3,8,0) = 0

!Initialize starting velocities
!Negative ensures planets begin orbit counter-clockwise
v(:,1,0) = 0
v(1,2,0) = 29780. * (-1.) * cos(phi(1))
v(2,2,0) = 29780. * (-1.) * sin(phi(1))
v(3,2,0) = 0
v(1,3,0) = 24070. * (-1.) * cos(phi(2))
v(2,3,0) = 24070. * (-1.) * sin(phi(2))
v(3,3,0) = 0
v(1,4,0) = 13060. * (-1.) * cos(phi(3))
v(2,4,0) = 13060. * (-1.) * sin(phi(3))
v(3,4,0) = 0
v(1,5,0) = 9680. * (-1.) * cos(phi(4))
v(2,5,0) = 9680. * (-1.) * sin(phi(4))
v(3,5,0) = 0
v(1,6,0) = 6800. * (-1.) * cos(phi(5))
v(2,6,0) = 6800. * (-1.) * sin(phi(5))
v(3,6,0) = 0
v(1,7,0) = 5430. * (-1.) * cos(phi(6))
v(2,7,0) = 5430. * (-1.) * sin(phi(6))
v(3,7,0) = 0
v(1,7,0) = 5000. * (-1.) * cos(phi(7))
v(2,7,0) = 5000. * (-1.) * sin(phi(7))
v(3,7,0) = 0

!Find centre of mass and velocity
DO i = 1, n-1
    COM(:) = COM(:) + r(:,i,0)*M(i)
    COV(:) = COV(:) + v(:,i,0)*M(i)
    Mtot = Mtot + M(i)
END DO

!Correct for COM and COV
DO i = 1,n
    r(:,i,0) = r(:,i,0) - COM/Mtot
    v(:,i,0) = v(:,i,0) - COV/Mtot

    !Work out absolute velocity for KE calculation
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
                a(:,i,0) = a(:,i,0) + (G * M(j) * (r(:,j,0) - r(:,i,0)) * (s_sq(i,j)**(-1.5)))

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
WRITE(6,*) "Initial Energy of System (J)"
WRITE(6,*) ""
!Total combined initial energy
WRITE(6,*) E_0

WRITE(6,*) ""
WRITE(6,*) "Would you like to save this run? (y/n)"
READ *, lg !Store user input to decide whether to log data
WRITE(6,*) ""

WRITE(6,*) "How many years would you like to simulate for?"
READ *, num_yr !Store user input to calculate time to exit loop
WRITE(6,*) ""

!saves repeat calculation at end of ABM loop
exittime = num_yr * 31536000.

!Open Log Files to write data to if requested
IF (lg == 'y') THEN
    OPEN(2,file = 'Sun_Motion.csv')
    OPEN(3,file = 'Earth_Motion.csv')
    OPEN(4,file = 'Jupiter_Motion.csv')
    OPEN(7,file = 'Mars_Motion.csv')
    OPEN(8,file = 'Saturn_Motion.csv')
    OPEN(9,file = 'Uranus_Motion.csv')
    OPEN(10,file = 'Neptune_Motion.csv')
    OPEN(11,file = 'Wandering_Motion.csv')
END IF


!================================================================================================================

!Bootstrap
DO k = 0,5

    !Find new position, r, after time-step
    !Comes at the start of the array to get the new position before other variables updated
    r(:,:,0) = r(:,:,0) + v(:,:,0)*dtmin + 0.5*a(:,:,0)*dtmin**2

    !Shift previous values down arrays
    DO i = -3,-1
        a(:,:,i) = a(:,:,i+1)
        v(:,:,i) = v(:,:,i+1)
    END DO

    DO i = 1,n
        DO j = 1,n
        !skip loop if i & j same (no force exerted by body on self)
        IF (i == j) CYCLE

            !Find absolute distance squared between bodies i & j
            s_sq(i,j) = ((r(1,j,0) - r(1,i,0))**2 + (r(2,j,0) - r(2,i,0))**2 + (r(3,j,0) - r(3,i,0))**2)


            !Find acceleration on object i in each dimension
            a(:,i,0) = a(:,i,0) + G * M(j) * (r(:,j,0) - r(:,i,0))*(s_sq(i,j)**(-1.5))


        END DO
    END DO

    !Find new velocities after time-step
    v(:,:,0) = v(:,:,0) + 0.5*(a(:,:,-1) + a(:,:,0))*dtmin

    !Update current time
    time = time + dtmin
END DO

!==================================================================================================================
!Ensure that all bodies are done in the first iteration
dt(:) = dtmin
deltat(:) = 0

!ABM Method
DO
    time = time + dtmin
    deltat(:) = deltat(:) + dtmin

    DO z = 1,n
        !skip bodies when their individual timestep has not been reached
        IF (deltat(z) /= dt(z)) CYCLE

        !Update counter for body
        counter(z) = counter(z) + 1

        !clear temporary interpolated values
        rtemp(:,:) = 0
        s_sq(:,:) = 0


        !Predict position and velocity using ABM predictor
        r(:,z,1) = r(:,z,0) + ((dt(z)/24.) * (-9.*v(:,z,-3) + 37.*v(:,z,-2) - 59.*v(:,z,-1) + 55.*v(:,z,0)))
        v(:,z,1) = v(:,z,0) + ((dt(z)/24.) * (-9.*a(:,z,-3) + 37.*a(:,z,-2) - 59.*a(:,z,-1) + 55.*a(:,z,0)))

        DO i = 1,n
            !skip interpolation for body acceleration is being calculated for
            IF (i==z) CYCLE
            !interpolate position from most recently updated position
            rtemp(:,i) = r(:,i,0) + v(:,i,0)*deltat(i) + 0.5*a(:,i,0)*deltat(i)**2
        END DO
        !Use predicted position for body acceleration is being calculated for
        rtemp(:,z) = r(:,z,1)

        !Clear predicted acceleration to not factor into calculation
        a(:,z,1) = 0
        !Calculate acceleration for predicted position
        DO j = 1,n
            IF (z == j) CYCLE
            !calculate square of absolute distance
            s_sq(z,j) = ((rtemp(1,j) - rtemp(1,z))**2 + (rtemp(2,j) - rtemp(2,z))**2 + (rtemp(3,j) - rtemp(3,z))**2)

            !Calculate new acceleration
            a(:,z,1) = a(:,z,1) + G * M(j) * (rtemp(:,j) - rtemp(:,z))*(s_sq(z,j)**(-1.5))
        END DO

        !ABM corrector
        r(:,z,2) = r(:,z,0) + ((dt(z)/24) * (v(:,z,-2) - 5.*v(:,z,-1) + 19.*v(:,z,0) + 9.*v(:,z,1)))
        v(:,z,2) = v(:,z,0) + ((dt(z)/24) * (a(:,z,-2) - 5.*a(:,z,-1) + 19.*a(:,z,0) + 9.*a(:,z,1)))

        !Correct position for body acceleration is being calculated for
        rtemp(:,z) = r(:,z,2)

        !Recalculate acceleration for corrected position
        a(:,z,2) = 0
        DO j = 1,n
            IF (z == j) CYCLE
            s_sq(z,j) = ((rtemp(1,j) - rtemp(1,z))**2 + (rtemp(2,j) - rtemp(2,z))**2 + (rtemp(3,j) - rtemp(3,z))**2)

            !Calculate corrected acceleration
            a(:,z,2) = a(:,z,2) + G * M(j) * (rtemp(:,j) - rtemp(:,z))*(s_sq(z,j)**(-1.5))
        END DO

        !Shift previous values down arrays
        DO i = -6,-1
            a(:,z,i) = a(:,z,i+1)
            v(:,z,i) = v(:,z,i+1)
        END DO

        !Set corrected values as current position, velocity and acceleration
        r(:,z,0) = r(:,z,2)
        v(:,z,0) = v(:,z,2)
        a(:,z,0) = a(:,z,2)


        !Ensure enough data points are stored
        IF (counter(z) > 6) THEN
            !Check if error from predictor-corrector is small enough
            IF (errcoeff * MAXVAL(ABS(a(:,z,2) - a(:,z,1)) / (ABS(a(:,z,2)) + Small)) < (RelErr * 0.01)) THEN

                !Double timestep
                dt(z) = dt(z) * 2.

                !Omit alternate points to adjust timestep
                DO k = 1,3
                    a(:,z,-k) = a(:,z,-2*k)
                    v(:,z,-k) = v(:,z,-2*k)
                END DO
            END IF


            !Check if error from predictor-corrector is too large
            IF (errcoeff * MAXVAL(ABS(a(:,z,2) - a(:,z,1)) / (ABS(a(:,z,2)) + Small)) > RelErr) THEN

                !Halve timestep
                dt(z) = dt(z) * 0.5

                !Interpolate previous data points for predictor-corrector
                a_int(:,z,1) = (-5. * a(:,z,-4) + 28. * a(:,z,-3) - 70. * a(:,z,-2) + 140. * a(:,z,-1) + 35. * a(:,z,0)) / 128.
                a_int(:,z,2) = (3. * a(:,z,-4) - 20. * a(:,z,-3) + 90. * a(:,z,-2) + 60. * a(:,z,-1) -5. * a(:,z,0)) / 128.

                v_int(:,z,1) = (-5. * v(:,z,-4) + 28. * v(:,z,-3) - 70. * v(:,z,-2) + 140. * v(:,z,-1) + 35. * v(:,z,0)) / 128.
                v_int(:,z,2) = (3. * v(:,z,-4) - 20. * v(:,z,-3) + 90. * v(:,z,-2) + 60. * v(:,z,-1) -5. * v(:,z,0)) / 128.

                !Create new mesh
                !-5 and -6 unimportant, not used in calculation & overwritten before next time-step change
                a(:,z,-4) = a(:,z,-2)
                a(:,z,-3) = a_int(:,z,2)
                a(:,z,-2) = a(:,z,-1)
                a(:,z,-1) = a_int(:,z,1)

                v(:,z,-4) = v(:,z,-2)
                v(:,z,-3) = v_int(:,z,2)
                v(:,z,-2) = v(:,z,-1)
                v(:,z,-1) = v_int(:,z,1)

            END IF

            !Reset counter to let old data fill up
            counter(z) = 0
        END IF
        !update that timestep has been reached
        deltat(z) = 0
    END DO

    IF (MINVAL(dt) > dtmin) THEN
        IF (ALL(MOD(deltat,MINVAL(dt)) .EQ. 0)) then
            dtmin = MINVAL(dt)
        END IF
    ELSEIF (MINVAL(dt)<dtmin) then
        dtmin = MINVAL(dt)
    END IF

    !Write every 500th step to log files
    IF (MOD(stepno, 500) == 0) THEN

    !Write every synched step to log files
    !IF (MAXVAL(deltat) == 0) THEN
        IF (lg == 'y') THEN
            !WRITE(6,*) stepno
            DO i = 1,n
                DO j = 1,n
                    IF (i==j) CYCLE
                    !Calculate absolute distance in each direction for logging files
                    !dist(:, i, j) = r(:,j,0) - r(:,i,0)
                    dist(:, i, j) = rtemp(:,j) - rtemp(:,i)
                END DO
            END DO
            WRITE(2,*) time, ',', dist(1,1,1), ',', dist(2,1,1), ',', dist(3,1,1)
            WRITE(3,*) time, ',', dist(1,1,2), ',', dist(2,1,2), ',', dist(3,1,2)
            WRITE(4,*) time, ',', dist(1,1,3), ',', dist(2,1,3), ',', dist(3,1,3)
            WRITE(7,*) time, ',', dist(1,1,4), ',', dist(2,1,4), ',', dist(3,1,4)
            WRITE(8,*) time, ',', dist(1,1,5), ',', dist(2,1,5), ',', dist(3,1,5)
            WRITE(9,*) time, ',', dist(1,1,6), ',', dist(2,1,6), ',', dist(3,1,6)
            WRITE(10,*) time, ',', dist(1,1,7), ',', dist(2,1,7), ',', dist(3,1,7)
            WRITE(11,*) time, ',', dist(1,1,8), ',', dist(2,1,8), ',', dist(3,1,8)
        END IF
    END IF

    !Update elapsed time

    stepno = stepno + 1

    !Second condition ensures bodies synched at end
    IF (time > exittime) THEN
        IF (MAXVAL(deltat) == 0) EXIT
    END IF
END DO

!Shut all logging files
IF (lg == 'y') THEN
    CLOSE(2)
    CLOSE(3)
    CLOSE(4)
    CLOSE(7)
    CLOSE(8)
    CLOSE(9)
    CLOSE(10)
    CLOSE(11)
END IF

WRITE(6,*) "Final time (yrs)", time / 31536000

WRITE(6,*) ""
WRITE(6,*) "Final Conditions"
WRITE(6,*) "Position Vectors xyz (AU)"
WRITE(6,*) ""
!Print final coordinates of all bodies
DO i = 1,n
    WRITE(6,*) (r(:,i,0)/AU)
END DO

WRITE(6,*) ""
WRITE(6,*) "Min timestep", dtmin
DO i = 1,n
    WRITE(6,*) i, "dt", dt(i)
END DO

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
DO k = 1,n
    WRITE(6,*) sqrt(s_sq(1,k)) /AU
END DO

KE = 0
pot = 0
DO i = 1,n
    KE = KE + 0.5*M(i)*absv_sq(i)
    pot = pot + GPE(i)
END DO

WRITE(6,*) "Energy Ratio", KE / pot




    !Subroutine to seed random number generation
    !From https://gcc.gnu.org/onlinedocs/gcc-4.6.1/gfortran/RANDOM_005fSEED.html
    Contains
    subroutine init_random_seed()

      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))

      CALL SYSTEM_CLOCK(COUNT=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)

      DEALLOCATE(seed)
    end


END PROGRAM solar_sim
