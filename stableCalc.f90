!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                          !!!
!!! Calculate stable nuclei abundances for a given tEvents.in file.          !!!
!!!                                                                          !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM stableCalc
    USE MPI
    IMPLICIT NONE
    
    ! MPI variables
    INTEGER::rank, nProc, ierror
    
    ! Program variables
    DOUBLE PRECISION, ALLOCATABLE::tEvents(:, :), tArray(:)
    DOUBLE PRECISION, ALLOCATABLE::abundArray(:, :), prodFactor(:, :)
    DOUBLE PRECISION::valFactor
    INTEGER::uni, nTimes, lenEvent, nEvents, ii, jj, redNEvents
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
    
    ! Now read the tEvents
    uni = 16
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
    OPEN(UNIT = uni, FILE = "StableOutput"//TRIM(ADJUSTL(sRank))//".txt")
    
    ! Write temporal array first and leave a space
    WRITE(uni, *) tArray
    WRITE(uni, *)
    
    ! Calculate
    jj = 1
    DO ii = 1, nEvents
        ! Divide events evenly
        IF (rank.NE.MOD(ii - 1, nProc)) CYCLE
        
        CALL stableAbund(abundArray(:, jj), tEvents(:, jj), tArray, prodFactor(:, jj))
        jj = jj + 1
    END DO
    
    ! Write
    jj = 1
    DO ii = 1, nEvents
        ! Divide events evenly
        IF (rank.NE.MOD(ii - 1, nProc)) CYCLE
        
        WRITE(uni, *) ii, abundArray(:, jj)
        jj = jj + 1
    END DO
    
    CLOSE(uni)
    
    DEALLOCATE(tEvents, tArray, abundArray, prodFactor)
    CALL MPI_FINALIZE(ierror)
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This function generates the decaying abundances array.                   !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -abundArray, the empty abundances array.                                 !!!
!!! -tEvents, the array with the polluting events times.                     !!!
!!! -tArray, the array with the sampling temporal points.                    !!!
!!! -prodFactor, the array with the production factor.                       !!!
!!!                                                                          !!!
!!! At the output, abundArray will be updated.                               !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE stableAbund(abundArray, tEvents, tArray, prodFactor)
    IMPLICIT NONE
    
    ! Input
    DOUBLE PRECISION::abundArray(:), tEvents(:), tArray(:), prodFactor(:)
    
    ! Local
    DOUBLE PRECISION::currT, prevT, dt, val, thisEvent, nextEvent, minTime
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
        
        IF (thisEvent.LE.currT) THEN
            
            ! Add each event up to currT
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
                    val = (prodFactor(iiEvent) + val)
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
    
END SUBROUTINE stableAbund

END PROGRAM stableCalc
