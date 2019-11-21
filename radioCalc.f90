!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                          !!!
!!! Calculate radioactive nuclei and abundances for a given tEvents.in file. !!!
!!! Also takes a list of taus for the calculation.                           !!!
!!!                                                                          !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM radioCalc
    USE MPI
    IMPLICIT NONE
    
    ! MPI variables
    INTEGER::rank, nProc, ierror
    
    ! Program variables
    DOUBLE PRECISION, ALLOCATABLE::taus(:), tEvents(:, :), tArray(:)
    DOUBLE PRECISION, ALLOCATABLE::abundArray(:, :), prodFactor(:, :)
    DOUBLE PRECISION::valFactor
    INTEGER::uni, nTau, nTimes, lenEvent, nEvents, ii, jj, kk, redNEvents
    INTEGER::lenFactor, nFactor
    CHARACTER(5)::sRank
    LOGICAL::isMaster
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Initialize MPI
    CALL MPI_INIT(ierror)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nProc, ierror)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
    
    ! Identify master
    isMaster = (rank.EQ.0)
    
    ! Read tau parameters
    uni = 16
    OPEN(UNIT = uni, FILE = "tauList.in")
    
    READ(uni, *) nTau
    ALLOCATE(taus(nTau))
    READ(uni, *) taus
    
    CLOSE(UNIT = uni)
    
    ! Now read the tEvents
    OPEN(UNIT = uni, FILE = "tEvents.in")
    
    READ(uni, *) lenEvent, nEvents
    ! Distribute the events among the processes to conserve memory
    redNEvents = 0
    DO ii = 1, nEvents
        IF (rank.NE.MOD(ii - 1, nProc)) CYCLE
        
        redNEvents = redNEvents + 1
    END DO
    
    ! Allocate tEvents with the reduced size redNEvents
    ALLOCATE(tEvents(lenEvent, redNEvents))
    
    ! Now read the information
    jj = 1
    DO ii = 1, nEvents
        IF (rank.EQ.MOD(ii - 1, nProc)) THEN
            READ(uni, *) tEvents(:, jj)
            jj = jj + 1
        ELSE
            READ(uni, *)
        END IF
    END DO
    
    CLOSE(UNIT = uni)
    
    ! Now read the prodFactor
    OPEN(UNIT = uni, FILE = "prodFactor.in")
    
    READ(uni, *) lenFactor, nFactor
    
    ! If there is only one production factor, then we fill the array with that
    ! value. We know is only one if lenFactor = nFactor = 1
    IF ((lenFactor.NE.1).AND.(lenFactor.NE.lenEvent)) THEN
        PRINT*, "Error! There should be the same number of events"
        PRINT*, " and production factors!"
        STOP
    END IF
    IF ((nFactor.NE.1).AND.(nFactor.NE.nEvents)) THEN
        PRINT*, "Error! There should be the same number of events"
        PRINT*, " and production factors!"
        STOP
    END IF
    IF (lenFactor.EQ.1) THEN
        READ(uni, *) valFactor
    END IF
    
    ! Allocate prodFactor with the reduced size redNEvents
    ALLOCATE(prodFactor(lenEvent, redNEvents))
    
    ! Now read the information
    IF (lenFactor.EQ.1) THEN
        prodFactor = valFactor
    ELSE
        jj = 1
        DO ii = 1, nEvents
            IF (rank.EQ.MOD(ii - 1, nProc)) THEN
                READ(uni, *) prodFactor(:, jj)
                jj = jj + 1
            ELSE
                READ(uni, *)
            END IF
        END DO
    END IF
    
    CLOSE(UNIT = uni)
    
    ! Allocate and read tArray, allocate abundArray
    OPEN(UNIT = uni, FILE = "tArray.in")
    
    READ(uni, *) nTimes
    ALLOCATE(tArray(nTimes), abundArray(nTimes, redNEvents))
    READ(uni, *) tArray
    
    CLOSE(UNIT = uni)
    
    ! Open output files. Write the rank to the filename.
    WRITE(sRank, '(I5)') rank
    OPEN(UNIT = uni, FILE = "Output"//TRIM(ADJUSTL(sRank))//".txt")
    
    ! Write temporal array first and leave a space
    WRITE(uni, *) tArray
    WRITE(uni, *)
    
    ! Do one separate run for each tau
    DO ii = 1, nTau
        WRITE(uni, *) "#", taus(ii)
        
        ! Calculate
        kk = 1
        DO jj = 1, nEvents
            ! Divide events evenly
            IF (rank.NE.MOD(jj - 1, nProc)) CYCLE
            
            CALL decayingAbund(abundArray(:, kk), tEvents(:, kk), tArray, &
                               taus(ii), prodFactor(:, kk))
            
            kk = kk + 1
        END DO
        
        ! Write
        kk = 1
        DO jj = 1, nEvents
            ! Divide events evenly
            IF (rank.NE.MOD(jj - 1, nProc)) CYCLE
            
            WRITE(uni, *) jj, abundArray(:, kk)
            kk = kk + 1
        END DO
    END DO
    
    CLOSE(uni)
    
    DEALLOCATE(taus, tEvents, tArray, abundArray, prodFactor)
    CALL MPI_FINALIZE(ierror)
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This function generates the decaying abundances array.                   !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -abundArray, the empty abundances array.                                 !!!
!!! -tEvents, the array with the polluting events times.                     !!!
!!! -tArray, the array with the sampling temporal points.                    !!!
!!! -tau, the tau value.                                                     !!!
!!! -prodFactor, the array with the production factor.                       !!!
!!!                                                                          !!!
!!! At the output, abundArray will be updated.                               !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE decayingAbund(abundArray, tEvents, tArray, tau, prodFactor)
    IMPLICIT NONE
    
    ! Input
    DOUBLE PRECISION::abundArray(:), tEvents(:), tArray(:), tau, prodFactor(:)
    
    ! Local
    DOUBLE PRECISION::invTau, currT, prevT, dt, val, thisEvent, nextEvent
    DOUBLE PRECISION::minTime
    INTEGER::ii, tLen, iiEvent, lenEvent
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Set up first point
    iiEvent = 1
    IF (tEvents(1).LE.0.D0) THEN
        abundArray(1) = prodFactor(1)
        iiEvent = 2
    ELSE
        abundArray(1) = 0.D0
    END IF
    
    ! Store inverted tau
    invTau = 1.0/tau
    
    ! Store lengths
    tLen = SIZE(tArray)
    lenEvent = SIZE(tEvents)
    
    ! Now calculate all of the others
    prevT = tArray(1)
    DO ii = 2, tLen
        currT = tArray(ii)
        dt = currT - prevT
        val = abundArray(ii - 1)
        
        ! Define event
        IF (iiEvent.LE.lenEvent) thisEvent = tEvents(iiEvent)
        
        ! Decay the dt if there is no influx
        IF (thisEvent.GT.currT) THEN
            val = val*EXP(-dt*invTau)
        ELSE
        ! Decay taking into account the intermediate events
            
            ! First decay up to new event (because there is at least one)
            val = val*EXP(-(thisEvent - prevT)*invTau)
            
            ! Now decay each event up to currT
            DO WHILE (thisEvent.LT.currT)
                
                ! Calculate nextEvent
                IF (iiEvent.LT.lenEvent) THEN
                    nextEvent = tEvents(iiEvent + 1)
                ELSE
                    nextEvent = 10*tArray(tLen)
                END IF
                
                ! Take the minimum
                minTime = MINVAL((/nextEvent, currT/))
                
                ! Make sure that we never have lower than 0 abundance
                IF ((prodFactor(iiEvent) + val).LT.0) THEN
                    val = 0
                ELSE
                    val = (prodFactor(iiEvent) + val)*EXP(-(minTime - thisEvent)*invTau)
                END IF
                
                ! Advance time
                iiEvent = iiEvent + 1
                thisEvent = nextEvent
            END DO
        END IF
        
        ! Add new value
        abundArray(ii) = val
        prevT = currT
    END DO
    
END SUBROUTINE decayingAbund

END PROGRAM radioCalc
