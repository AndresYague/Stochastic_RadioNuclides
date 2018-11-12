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
    REAL, ALLOCATABLE::tEvents(:, :), tArray(:)
    REAL, ALLOCATABLE::abundArray(:, :)
    INTEGER, ALLOCATABLE::eventsTags(:)
    INTEGER::uni, nTimes, lenEvent, nEvents, ii, jj, redNEvents
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
    
    ! Allocate the reduced sizes
    ALLOCATE(tEvents(lenEvent, redNEvents), eventsTags(redNEvents))
    
    ! Now read the information
    jj = 1
    DO ii = 1, nEvents
        IF (rank.NE.MOD(ii - 1, nProc)) CYCLE
        
        READ(uni, *) tEvents(:, jj)
        eventsTags(jj) = ii
        
        jj = jj + 1
    END DO
    
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
    ii = 1
    DO jj = 1, nEvents
        ! Divide events evenly
        IF (rank.NE.MOD(jj - 1, nProc)) CYCLE
        
        CALL stableAbund(abundArray(:, ii), tEvents(:, ii), tArray)
    END DO
    
    ! Write
    ii = 1
    DO jj = 1, nEvents
        ! Divide events evenly
        IF (rank.NE.MOD(jj - 1, nProc)) CYCLE
        
        WRITE(uni, *) eventsTags(ii), abundArray(:, ii)
        ii = ii + 1
    END DO
    
    DEALLOCATE(tEvents, tArray, abundArray, eventsTags)
    CALL MPI_FINALIZE(ierror)
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This function generates the decaying abundances array.                   !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -abundArray, the empty abundances array.                                 !!!
!!! -tEvents, the array with the polluting events times.                     !!!
!!! -tArray, the array with the sampling temporal points.                    !!!
!!!                                                                          !!!
!!! At the output, abundArray will be updated.                               !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE stableAbund(abundArray, tEvents, tArray)
    IMPLICIT NONE
    
    ! Input
    REAL::abundArray(:), tEvents(:), tArray(:)
    
    ! Local
    REAL::currT, prevT, dt, val, thisEvent, nextEvent
    REAL::minTime
    INTEGER::ii, tLen, iiEvent, lenEvent
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Set up first point
    iiEvent = 1
    IF (tEvents(1).LE.0.D0) THEN
        abundArray(1) = 1.D0
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
                val = val + 1
                
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
