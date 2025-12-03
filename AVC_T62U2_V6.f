C!  *******************************************************************
      PROGRAM AVC_T62U2_V6
C!  -------------------------------------------------------------------
      IMPLICIT NONE
C!    ----- Parameters: Version Date and number
      CHARACTER VDATE*16
      REAL*4 VINPUT  ! Version number for command files YYMMDD.HHMM 
C!                     (date and time  file structure was last changed)
C!    #################################################################
      PARAMETER(VDATE='23-JUL-2025 (HR)',VINPUT=231004.1743d0) ! ############
C!    #################################################################
C!    Revision history:
C!    ------------------------------------------------------------------------------         
C!    *********** Version V6 ******************************************
C!    23-JUL-2025 Fixed bug in this context with the last track. 
C!    17-JUN-2025 Added output of track summaries. 
C!    *********** Version V5 ******************************************
C!    16-JUN-2025 Corrected several bug in the scoring of the totals 
C!                and the output (introduced on 10-MAR-2025). 
C!    11-JUN-2025 Introduced workaround for segmentation fault caused by
C!                the high number of events in PHITS data at low energy. 
C!    10-JUN-2025 Corrected a bug with the newly introduced arrays of 
C!                energy imparted, total number of ions and clusters 
C!                track. Currently not yet output.
C!    10-MAR-2025 Changed scoring of clusters per track and globally to up
C!                to the maximum ICS to enable output of more detailed ICSDs.
C!    *********** Version V4 *******************************************
C!    09-MAR-2025 Added headers for version number in version history. 
C!                Reduced output to the last radial distance with entries.
C!    07-MAR-2025 Fixed handling data files that contain only ionizations,
C!                in which the number of saved events deviates from the number of 
C!                simulated events plus the case of Geant4-DNA Opt4 15 eV which 
C!                contain data from two simulation runs with the same event labelling. 
C!    06-MAR-2025 Added fix for handling of PTra data files that deviate from 
C!                expected structure: They have energy loss in 5th column and 
C!                energy deposit in a 7th column. 
C!    20-FEB-2025 Bug fix: labelling of columns in the header of output file was incorrect. 
C!    ------------------------------------------------------------------      
C!    *********** Version V3 *******************************************
C!    16-MAY-2024 Bug fix: When moving the first point of a new track to the top of the list
C!                after processing the previous track, the flag for ionization was omitted. 
C!                (Bug reported by Pavel.)
C!    ------------------------------------------------------------------   
C!    *********** Version V2 *******************************************
C!    28-JAN-2024 As decided during the last meeting of the task group
C!                on 22-Jan-2024, the code has been modified to read an
C!                additional flag indicating impact ionization (1) or
C!                electronic excitation with or without subsequent auto-
C!                ionization (2). At present this is just a provision
C!                for the final analysis to take place after a common 
C!                sharepoint will be available (again). 
C!                For the distributed analysis by participant, ICSD are 
C!                determined scoring only impact ionizations.
C!    ------------------------------------------------------------------       
C!    *********** Version V1 *******************************************
C!    07-DEC-2023 Bug fix: When moving the first point of a new track to the top of the list
C!                after processing the previous track, the value of E was not moved along. 
C!                (Bug reported by Zine.)
C!    11-NOV-2023 Fixed several bugs based on problems reported by Zine.
C!                #1: TREIMP and TRION were not normalized to NITER (in DO_AVC)
C!                #2: NPTS was not properly set for the last track.
C!                #3: In subroutine SORT62, the case that the input data are 
C!                    already ascending had not been properly handled. 
C!                #4: E was not read for the first point in a track.
C!                #5: Radial scoring of the target center was flawed (the
C!                    (variable RSITE2 was used instead of RSITE, i.e. the
C!                     square of the site radius instead of the radius.)
C!                Further changes: 
C!                - parameter NPMAX has been increased to 10000. 
C!                - NP has been renamed to NPTS.
C!                - In the main loop, the variable NPTS is now used 
C!                  directly instead of I. 
C!    22-OCT-2023 Bug fix: see inline comments starting with the revision 
C!                date in the corresponding code section below, at the end 
C!                of the 'Process_this_track' loop. 
C!    04-OCT-2023 First released version
C!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C!    ----- Parameters: 
      INTEGER*4 MAXFIL, MAXICS,MAXNRB,NPMAX,NTRMAX
      REAL*4 ONE,PI4BY3,TEN,TWO,ZERO
      PARAMETER (MAXFIL=64, MAXICS=20,MAXNRB=1001,NPMAX=100000,
     &           NTRMAX=100000,                                   ! 10-MAR-2025
     &           PI4BY3= 4.1887902,ONE=1.,TEN=10.,TWO=2.,ZERO=0.)
C!    ----- Input parameters for geometry
      REAL*4 DSITE  ! Diameter of spherical target 
      REAL*4 RMIN   ! Radius of innermost sphere for regional scoring
      REAL*4 RMAX   ! Radius of outermost sphere for regional scoring
      REAL*4 DLOGR  ! logarithmic bin size
C!    ----- Input parameters for scoring
      INTEGER*4 NITER  ! number of iterations per track
      INTEGER*4 ISHARE ! [1/0] Consider distribution of site volume over different shells (if applicable)
      INTEGER*4 ICSMAX ! max cluster size (for which F_n instead of P_n is determined)
C!    ----- Input parameters related to data
      INTEGER*4 NHEADL ! Number of header lines 
      INTEGER*4 NFILES ! Number of files to process
      INTEGER*4 ISORT  ! Flag indicating whether data need sorting 
C                      ! 0: Data don't need sorting
C                      ! 1: Data to be sorted and output to scratch file
C                      ! 2: Data to be sorted and output to file 'S_'//FILENM                           
C!    ----- Functions
      INTRINSIC EXP,ALOG
      CHARACTER TSTAMP*24,FLTSTR*16,INTSTR*8,DATEXT*80
      LOGICAL COLION ! 28-JAN-2024
C!    ----- Global variables
C!    Scoring
C!    Variables already declared under input variables
      LOGICAL LSHARE
      COMMON /CONSTS/ ICSMAX,LSHARE,NITER
C!    Geometry
      INTEGER*4 NRBINS ! Actual number of radial bins
      REAL*4 RRBIN2(MAXNRB),RRBHI2(MAXNRB),RRBLO2(MAXNRB),RSITE,RSITE2
      COMMON /GEOMET/ RRBIN2,RRBHI2,RRBLO2,RSITE,RSITE2,NRBINS
C!    Data format
C!     CODE: Character string encoding the meaning of the entries in a 
C!           line of the input file as follows:
C!           'T'     - number of the primary particle track
C!           'X','Y','Z' - x, y, and z coordinates of the transfer point 
C!           'E'     - energy deposit (if applicable)
C!           'I'     - ionization cluster size (if applicable)
C!           '%'     - additional data that are not used 
C!           Example: The string 'TZ%%EXY ' indicates that there are 7
C!                    entries in each line, of which the first is the 
C!                    track number, the second the z coordinate, the 
C!                    fifth the energy, and the sixth and seventh the x
C!                    and y coordinates
C!     IDT: Array of the columns indices of (1-3) x,y,z coordinates 
C!                   (4) energy deposit (if present), (5) track number, 
C!                   (6) number of ionizations in cluster (if present)
C!                   (7-8) are there for future use
      CHARACTER CODE*8
      INTEGER*2 IDT(8),NDT
      COMMON /RFORMT/IDT,NDT,CODE
C!    ----- Local scalars
      CHARACTER EMPTY*2,FILENM*80,OUTPUT*80,PREFIX*40,USERID*80
      INTEGER*4 I,ICS,IDEBUG,IDUMMY,IEV,IFILE,IOLD,IRB,IRMIN,ITR,J,
     &          MRBINS,NEV,NEVMAX,NOTIME,NPTS
      LOGICAL ASKINP, THERE
      LOGICAL IONFLG ! 28-JAN-2024
      LOGICAL ISPTRA ! 06-MAR-2025
      REAL*4 ELOSS   ! 06-MAR-2025
      REAL*4 DLNR,FNORM,RRBIN,TCION,TEIMP,TTCION,TEIMPT,UTCION,UTEIMP, 
     &       VCHECK,VSITE
C!    ----- Local arrays
      CHARACTER*80 FNAMES(MAXFIL)
      LOGICAL DEBUG(4)
      INTEGER ISION(NPMAX) ! 28-JAN-2024
      REAL*4 TTEIMP(NTRMAX),TTIONS(NTRMAX),TTICS(NTRMAX,MAXICS)  ! 10-MAR-2025
      REAL*4 X(NPMAX),Y(NPMAX),Z(NPMAX),E(NPMAX)
      REAL*4 RBIN(MAXNRB), EIMP(MAXNRB),UEIMP(MAXNRB),TREIMP(MAXNRB)
      REAL*4 CION(MAXNRB),UION(MAXNRB),VOLRB(MAXNRB),TRION(MAXNRB)
      REAL*4 PICS(MAXNRB,MAXICS),UICS(MAXNRB,MAXICS),
     &       TRICS(MAXNRB,MAXICS)
      REAL*4 RESU(2*MAXICS+1,2),TPICS(MAXICS),TTPICS(MAXICS),
     &       UTPICS(MAXICS)
      CHARACTER*12 COLHDR(MAXICS,2)
C!<<<<End declarations<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      
      
      PRINT*, 'Program AVC_T62U2_V4 Version: '//VDATE
C!  -------------------------------------------------------------------
      Read_INPUT: DO ! DUMMY block for readability: Read input 
        DEBUG(1)=.FALSE. ! Main program
        DEBUG(2)=.FALSE. ! CLUSTR main sections
        DEBUG(3)=.FALSE. ! CLUSTR IDIR loop
        DEBUG(4)=.FALSE. ! CLUSTR IPOS loop details
        
        PREFIX='AVC_0.0nm_'
        
C!      Check command file version
        PRINT*, 'Enter 0 for manual input or version (YYMMDD.HHMM) '//
     &          'of command file structure ' 
        READ(*,'(a)') FILENM
        PRINT*,FILENM
        READ(FILENM(1:11),*) VCHECK
        ASKINP=(VCHECK.EQ.ZERO)
        IF(.NOT.ASKINP) THEN
          IF(VCHECK.LT.VINPUT) THEN
            PRINT*, 'Command file structure ',VCHECK,' older than '//
     &              'current version ',VINPUT,'=> STOP.'
            STOP
          END IF       
          READ(*,*) FILENM
          IF(FILENM(1:9).NE.'AVC_T62U2') THEN
            PRINT*, 'Command file appears not to be for AVC_T62U2 '//
     &              'but for '//FILENM//'=> STOP.'
            STOP
          END IF    
        END IF

C!      Read debug options
        IF(ASKINP) PRINT*, 'Enter debug options flag or 0 for none' 
        READ(*,*) IDEBUG        
        IF(IDEBUG.GT.0) THEN
          IDUMMY=IDEBUG
          DO I=1,4
            DEBUG(I)=(MOD(IDUMMY,2).EQ.1)
            IDUMMY=IDUMMY/2
          END DO ! I=1,4
        ENDIF
        
C!      Read geometry parameters 
        IF(ASKINP) PRINT*, 'Enter site diameter in nm'
        READ(*,*) DSITE
        
        IF(ASKINP) PRINT*, 'Radius of innermost sphere (in nm) for '//
     &                     'regional scoring'
        READ(*,*) RMIN
        IF(ASKINP) PRINT*, 'Radius of outermost sphere (in nm) for '//
     &                     'regional scoring'
        READ(*,*) RMAX
        IF(ASKINP) PRINT*, 'Logarithmic bin size for radius (base 10)'
        READ(*,*) DLOGR

        IF(ASKINP) PRINT*, 'Enter number of iterations '
        READ(*,*) NITER

        IF(ASKINP) PRINT*, 'Consider distribution of site volume over '
     &                     //'different shells (if applicable)'      
        READ(*,*) ISHARE

        IF(ASKINP) PRINT*, 'Enter maximum ionization cluster size '//
     &                     ' (for output files)'
        READ(*,*) ICSMAX
        IF(ICSMAX.GT.MAXICS) THEN
          PRINT*,'Maximum ICS cannot exceed ',MAXICS,', which is used ',
     &           'as default instead of ',ICSMAX
          ICSMAX=MAXICS
        ENDIF
        
C!      Read user input
        IF(ASKINP) THEN
          PRINT*, 'Enter your name, institute, city, country'//
     &             ' (all in one line)'
          READ(*,'(a)') USERID        

          PRINT*, 'Enter 8 character code for data file '//
     &       'structure where ''T'' indicates the track ID,'//
     &       ' ''X'',''Y'' and ''Z'' the respective coordinates, '//
     &       '''E'' the energy deposit (if present) and ''&'' any '//
     &       'other data'
          READ(*,*) CODE

          PRINT*, 'Sort data first(0/1 - 2 for output to file)'
          READ(*,*) ISORT

          PRINT*, 'Number of header lines'
          READ(*,*) NHEADL

          PRINT*, 'Number of files to process'
          READ(*,*) NFILES
        ELSE
          READ(*,'(a)') EMPTY        
          READ(*,*) CODE   ! 8 character code for input file structure
          READ(*,*) ISORT  ! Sort data first(0/1 - 2 for output to file)

          READ(*,'(a)') EMPTY        
          READ(*,'(a)') USERID        
          READ(*,*) NHEADL ! Number of header lines
          READ(*,*) NFILES ! 'Number of files to process'
        ENDIF

        IF(NFILES.GT.MAXFIL) THEN
          PRINT*,'Maximum number of files is ',MAXFIL
          PRINT*,' ==> only the first ',MAXFIL,' files are processed.'
          PRINT*,'After this execution completes, you must run the '//
     &           'program again with the list of remaining files' 
          NFILES=MAXFIL
        ENDIF
        
        DO IFILE=1,NFILES
          IF(ASKINP) PRINT*, 'Name of file #',I
          READ(*,*) FNAMES(IFILE)
        END DO
        EXIT Read_INPUT
      END DO Read_INPUT ! DUMMY block for readability: Read input 
C!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      Initialize: DO ! DUMMY block for readability: Initialize 
        LSHARE=(ISHARE.EQ.1) ! Logical FLag for sharing sites over radial bins (if applicable)
C
        WRITE(*,'(2a,/)') 'Data file structure is ',CODE
        CALL RFINIT()
        RSITE=DSITE/2. ! radius of spherical site
        RSITE2=RSITE*RSITE ! Square of site radius
        VSITE=PI4BY3*RSITE**3 ! volume of the site
C
C!      Actual number of radial bins
        NRBINS=INT(ALOG10(RMAX/RMIN)/DLOGR) ! Number of radial bins requested
        IF(NRBINS.GT.MAXNRB) THEN
          PRINT*, 'WARNING: Requested number of radial bins exceeds '
     &            //'the maximum possible value!'
          PRINT*, 'Default to using the maximum possible value.'
          PRINT*, 'This means that RMAX is adjusted to ',
     &            RMIN*EXP(DLOGR*ALOG(10.)*(MAXNRB+1))
          NRBINS=MAXNRB
        END IF
C
C!      Radial bins
        DLNR=DLOGR*ALOG(TEN)
        IF(RMIN.LT.RSITE) THEN
          IRMIN=1
          DO WHILE(RMIN*EXP(DLNR*IRMIN).LT.RSITE)
            IRMIN=IRMIN+1
            NRBINS=NRBINS-1
          END DO
        ELSE
          IRMIN=0
        ENDIF
        
        DO IRB=1,NRBINS
          RBIN(IRB)=RMIN*EXP(DLNR*(IRB+IRMIN-1))
          RRBIN=RBIN(IRB)/RSITE ! outer bin radius in units of site radius
          RRBIN2(IRB)=RRBIN*RRBIN ! square of outer bin radius in units of site radius
          RRBLO2(IRB)=RRBIN2(IRB)-TWO*RRBIN+ONE ! square of relative radius offset by -1
          RRBHI2(IRB)=RRBIN2(IRB)+TWO*RRBIN+ONE ! square of relative radius offset by +1
          VOLRB(IRB)=PI4BY3*RBIN(IRB)**3 ! Volume of sphere 
          IF(IRB.GT.1) VOLRB(IRB)=VOLRB(IRB)-PI4BY3*RBIN(IRB-1)**3 ! radial bin volume
C!*          PRINT*, I,J,RBIN(IRB),VOLRB(IRB),VOLRB(IRB)/RBIN(IRB)**3
        END DO
C
C!      Output file prefix
        IF(DSITE.LT.10.0d0) THEN
          WRITE(PREFIX(5:7),'(f3.1)') DSITE
        ELSE
          IF(DSITE.LT.100.0d0) THEN
            WRITE(PREFIX(5:7),'(f3.0)') DSITE
          ELSE
            WRITE(PREFIX(5:7),'(i3)') INT(DSITE)
          END IF
        END IF
        IF(LSHARE) PREFIX='AVCs_'//PREFIX(5:39)
        PREFIX=TRIM(PREFIX)//'nmax='//TRIM(INTSTR(ICSMAX))//'_'
        IF(NITER.GT.1) PREFIX=TRIM(PREFIX)//'iter='
     &                        //TRIM(INTSTR(NITER))//'_'
C
C!      Column headers of output files
        DO ICS=1,ICSMAX
          IF(ICS.LT.ICSMAX) THEN
            COLHDR(ICS,1)='dP'//TRIM(INTSTR(ICS))//'/dV'
          ELSE
            COLHDR(ICS,1)='dF'//TRIM(INTSTR(ICS))//'/dV'
          ENDIF
          COLHDR(ICS,2)='u('//TRIM(COLHDR(ICS,1))//')'
        END DO
        EXIT Initialize
      END DO Initialize  ! DUMMY block for readability: Initialize
C!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      Main_Loop: DO IFILE=1,NFILES ! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
        FILENM=FNAMES(IFILE)
        IF(ISORT.EQ.0) THEN
          INQUIRE(FILE='S_'//FILENM,EXIST=THERE)
          IF(THERE) FILENM='S_'//FILENM
        ENDIF
        
        INQUIRE(FILE=FILENM,EXIST=THERE)
        IF(.NOT.THERE) THEN
          PRINT*,'File '//TRIM(FILENM)//' not found.'
          CYCLE Main_Loop
        ENDIF
        
C!      Open file 
        IF(ISORT.EQ.0) THEN
          OPEN(11,FILE=FILENM,STATUS='OLD')
        ELSE
          CALL SORT62(FILENM,ISORT)
        ENDIF
        WRITE(*,'(3a,$)') 'Process ',TRIM(FILENM),'.'
        
C! 28-JAN-2024 (HR) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
C!      Check file for presence of ionization flag        
        IONFLG=COLION(11,NHEADL)
        ISPTRA=(CODE(5:5).NE.'E') ! 06-MAR-2025
        PRINT*,'<#help ionflg isptra<<<<<<<<<<<<<',IONFLG,ISPTRA  ! 06-MAR-2025
                
C! 28-JAN-2024 (HR) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<       
        
C!      First zero all counters     
C!16-JUN-2025 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        TEIMP=ZERO
        UTEIMP=ZERO ! 16-JUN-2025
        TCION=ZERO
        UTCION=ZERO ! 16-JUN-2025
        DO ICS=1,MAXICS
          TPICS(ICS)=ZERO
          UTPICS(ICS)=ZERO  ! 16-JUN-2025
        END DO
C!16-JUN-2025 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        DO IRB=1,NRBINS
          EIMP(IRB)=ZERO
          UEIMP(IRB)=ZERO
          TREIMP(IRB)=ZERO
          CION(IRB)=ZERO
          UION(IRB)=ZERO
          TRION(IRB)=ZERO
          DO ICS=1,MAXICS
            PICS(IRB,ICS)=ZERO
            UICS(IRB,ICS)=ZERO
            TRICS(IRB,ICS)=ZERO
          END DO
        END DO

C! 10-MAR-2025 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>       
        DO ITR=1,NTRMAX
          TTEIMP(ITR)=ZERO
          TTIONS(ITR)=ZERO
          DO ICS=1,MAXICS
            TTICS(ITR,ICS)=ZERO
          END DO
        END DO
C! 10-MAR-2025 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<     

C!      Initialize scalars
        NPTS=1 ! 11-NOV-2023 replaced I with NPTS
        NEV=1
        NEVMAX=0 ! 07-MAR-2025
        
        IF(ISPTRA) THEN  ! 06-MAR-2025
          READ(11,*) IEV,X(NPTS),Y(NPTS),Z(NPTS),ELOSS,ISION(NPTS),
     &               E(NPTS) ! 06-MAR-2025 (HR)
        ELSE IF(IONFLG) THEN ! 28-JAN-2024 (HR)  ! 06-MAR-2025
          READ(11,*) IEV,X(NPTS),Y(NPTS),Z(NPTS),E(NPTS),ISION(NPTS) ! 28-JAN-2024 (HR)
        ELSE ! 28-JAN-2024 (HR)
          READ(11,*) IEV,X(NPTS),Y(NPTS),Z(NPTS),E(NPTS) ! 11-NOV-2023 replaced I with NPTS; BUG fix: added E(NPTS)
          ISION(NPTS)=1 ! 28-JAN-2024 (HR)
        ENDIF ! 28-JAN-2024 (HR)
        IOLD=IEV
        
        Process_this_track: DO ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          IF(MOD(NEV,100).EQ.0) WRITE(*,'(1H.,a,$)') ''
          DO WHILE(IEV.EQ.IOLD)
            NPTS=NPTS+1
            IF(ISPTRA) THEN  ! 06-MAR-2025
              READ(11,*,END=20) IEV,X(NPTS),Y(NPTS),Z(NPTS),ELOSS,
     &                          ISION(NPTS),E(NPTS) ! 06-MAR-2025 (HR)
            ELSE IF(IONFLG) THEN ! 28-JAN-2024 (HR)
              READ(11,*,END=20) IEV,X(NPTS),Y(NPTS),Z(NPTS),E(NPTS),
     &                          ISION(NPTS) ! 28-JAN-2024 (HR)
            ELSE ! 28-JAN-2024 (HR)
              READ(11,*,END=20) IEV,X(NPTS),Y(NPTS),Z(NPTS),E(NPTS) ! 11-NOV-2023 replaced I with NPTS; BUG fix: added E(NPTS)
              ISION(NPTS)=1 ! 28-JAN-2024 (HR)
            ENDIF ! 28-JAN-2024 (HR)
          END DO
          NPTS=NPTS-1  ! 11-NOV-2023 replaced I with NPTS
          
C!        Tally clusters, energy imparted and ionizations
          CALL DO_AVC(NPTS,X,Y,Z,E,ISION,TRICS,TREIMP,TRION) ! 28-JAN-2024
C! 28-JAN-2024          CALL DO_AVC(NPTS,X,Y,Z,E,TICS,TEIM,TION)

C         Update global counters
C!16-JUN-2025 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          TEIMPT=ZERO
          TTCION=ZERO
          DO ICS=1,MAXICS
            TTPICS(ICS)=ZERO
          END DO
C!16-JUN-2025 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          DO IRB=1,NRBINS
C!16-JUN-2025 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            TEIMPT=TEIMPT+TREIMP(IRB)
            TTCION=TTCION+TRION(IRB)
            DO ICS=1,MAXICS
              TTPICS(ICS)=TTPICS(ICS)+TRICS(IRB,ICS)
            END DO
C!16-JUN-2025 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            EIMP(IRB)=EIMP(IRB)+TREIMP(IRB)
            UEIMP(IRB)=UEIMP(IRB)+TREIMP(IRB)*TREIMP(IRB)
            TREIMP(IRB)=ZERO
            CION(IRB)=CION(IRB)+TRION(IRB)
            UION(IRB)=UION(IRB)+TRION(IRB)*TRION(IRB)
            TRION(IRB)=ZERO
            DO ICS=1,MAXICS
              PICS(IRB,ICS)=PICS(IRB,ICS)+TRICS(IRB,ICS)
              UICS(IRB,ICS)=UICS(IRB,ICS)+TRICS(IRB,ICS)*TRICS(IRB,ICS)
              TRICS(IRB,ICS)=ZERO
            END DO
          END DO
C! 10-MAR-2025 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
          IF(NEV.LE.NTRMAX) THEN 
            TTEIMP(NEV)=TEIMPT
            TTIONS(NEV)=TTCION
            DO ICS=1,MAXICS
              TTICS(NEV,ICS)=TTPICS(ICS)
            END DO
          ENDIF
C! 10-MAR-2025 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C!16-JUN-2025 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          TEIMP=TEIMP+TEIMPT
          UTEIMP=UTEIMP+TEIMPT*TEIMPT
          TCION=TCION+TTCION
          UTCION=UTCION+TTCION*TTCION
          DO ICS=1,MAXICS
            TPICS(ICS)=TPICS(ICS)+TTPICS(ICS)
            UTPICS(ICS)=UTPICS(ICS)+TTPICS(ICS)*TTPICS(ICS)
          END DO
C!16-JUN-2025 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          
C!        HR 22-Oct-2023: There was a bug here previously, which has been 
C!        fixed by using variable IRB in the preceding outer DO loop  
C!        instead of I as before. When I was used as counter in the loop, 
C!        I was not pointing to the first datapoint of the new event in 
C!        the following three lines. Therefore the first datapoint of the 
C!        new event was replaced by a point of the previous track.
          X(1)=X(NPTS+1) ! 11-NOV-2023 replaced I with NPTS+1 (last read data was first in new track)
          Y(1)=Y(NPTS+1) ! 11-NOV-2023 replaced I with NPTS+1 (last read data was first in new track)
          Z(1)=Z(NPTS+1) ! 11-NOV-2023 replaced I with NPTS+1 (last read data was first in new track)
          E(1)=E(NPTS+1) ! 07-DEC-2023 added, since this was missing (according to ZINE)
          ISION(1)=ISION(NPTS+1) ! 16-MAY-2024 added, since this was missing (faulty output reported by Pavel)
          NPTS=1  ! 11-NOV-2023 replaced I with NPTS
          NEV=NEV+1
          IF(IEV.LT.IOLD) NEVMAX=NEVMAX+IOLD  ! 07-MAR-2025
          IOLD=IEV
        
        END DO Process_this_track 
C!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
  20    CONTINUE
        PRINT*,'Done. ',NEV,' sequences of equal eventID processed.'
        CLOSE(11)
C
C!      Process the last track --------------------------
C
C       Tally clusters, energy imparted and ionizations
        NPTS=NPTS-1  ! 11-NOV-2023 added this line for BUG fix: Previously, NPTS was the value from the last but one track
        CALL DO_AVC(NPTS,X,Y,Z,E,ISION,TRICS,TREIMP,TRION) ! 28-JAN-2024
C! 28-JAN-2024        CALL DO_AVC(NPTS,X,Y,Z,E,TICS,TEIM,TION)
C       Final update of global counters 
C!16-JUN-2025 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        TEIMPT=ZERO
        TTCION=ZERO
        DO ICS=1,MAXICS
          TTPICS(ICS)=ZERO
        END DO
C!16-JUN-2025 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        DO IRB=1,NRBINS
C!16-JUN-2025 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          TEIMPT=TEIMPT+TREIMP(IRB)
C!23-JUL-2025          TCION=TCION+TRION(IRB) 
          TTCION=TTCION+TRION(IRB) ! 23-JUL-2025
          DO ICS=1,MAXICS
            TTPICS(ICS)=TTPICS(ICS)+TRICS(IRB,ICS)
          END DO
C!16-JUN-2025 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          EIMP(IRB)=EIMP(IRB)+TREIMP(IRB)
          UEIMP(IRB)=UEIMP(IRB)+TREIMP(IRB)*TREIMP(IRB)
          CION(IRB)=CION(IRB)+TRION(IRB)
          UION(IRB)=UION(IRB)+TRION(IRB)*TRION(IRB)
          DO ICS=1,MAXICS
            PICS(IRB,ICS)=PICS(IRB,ICS)+TRICS(IRB,ICS)
            UICS(IRB,ICS)=UICS(IRB,ICS)+TRICS(IRB,ICS)*TRICS(IRB,ICS)
          END DO
        END DO

C! 10-MAR-2025 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
          IF(NEV.LE.NTRMAX) THEN 
            TTEIMP(NEV)=TEIMPT
            TTIONS(NEV)=TTCION
            DO ICS=1,MAXICS
              TTICS(NEV,ICS)=TTPICS(ICS)
            END DO
          ENDIF
C! 10-MAR-2025 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

C!16-JUN-2025 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        TEIMP=TEIMP+TEIMPT
        UTEIMP=UTEIMP+TEIMPT*TEIMPT
        TCION=TCION+TTCION
        UTCION=UTCION+TTCION*TTCION
        DO ICS=1,MAXICS
          TPICS(ICS)=TPICS(ICS)+TTPICS(ICS)
          UTPICS(ICS)=UTPICS(ICS)+TTPICS(ICS)*TTPICS(ICS)
        END DO
C!16-JUN-2025 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

C!      Calculate final results ###############################      
        Results: DO ! Dummy Loop for output

C!        Get max number of tracks !! 07-MAR_2025
          NEVMAX=NINT(REAL(IEV+NEVMAX)**2/REAL(NEV))  ! 07-MAR-2025

C!        Calculate global scores
*          DO IRB=1,NRBINS
*            TEIMP=TEIMP+EIMP(IRB)
*            UTEIMP=UTEIMP+UEIMP(IRB)  ! 10-MAR-2025
*            TCION=TCION+CION(IRB)
*            UTCION=UTCION+UION(IRB)  ! 10-MAR-2025
*            DO ICS=1,MAXICS
*              TPICS(ICS)=TPICS(ICS)+PICS(IRB,ICS)
*              UTPICS(ICS)=UTPICS(ICS)+UICS(IRB,ICS) ! 10-MAR-2025
*            END DO
*          END DO

C!        Normalize per track 
C! 07-MAR-2025          FNORM=ONE/NEV
          FNORM=ONE/NEVMAX  ! 07-MAR-2025
          TEIMP  = TEIMP*FNORM
          UTEIMP = UTEIMP*FNORM
          TCION  = TCION*FNORM
          UTCION = UTCION*FNORM
          DO ICS=1,MAXICS
            TPICS(ICS)=TPICS(ICS)*FNORM
            UTPICS(ICS)=UTPICS(ICS)*FNORM
          END DO
          
          DO IRB=1,NRBINS
            EIMP(IRB)=EIMP(IRB)*FNORM
            UEIMP(IRB)=UEIMP(IRB)*FNORM
            CION(IRB)=CION(IRB)*FNORM
            UION(IRB)=UION(IRB)*FNORM
            DO ICS=1,MAXICS
              PICS(IRB,ICS)=PICS(IRB,ICS)*FNORM
              UICS(IRB,ICS)=UICS(IRB,ICS)*FNORM
            END DO
          END DO


C!        Calculate standard deviations
          DO IRB=1,NRBINS
            UEIMP(IRB)=SQRT((UEIMP(IRB)-EIMP(IRB)*EIMP(IRB))*FNORM)
            UION(IRB)=SQRT((UION(IRB)-CION(IRB)*CION(IRB))*FNORM)
            DO ICS=1,MAXICS
              UICS(IRB,ICS)=SQRT((UICS(IRB,ICS)-PICS(IRB,ICS)**2)*FNORM)
            END DO
          END DO

C!10-MAR-2025 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          UTEIMP=SQRT((UTEIMP-TEIMP*TEIMP))
          UTCION=SQRT((UTCION-TCION*TCION))
          DO ICS=1,MAXICS
            UTPICS(ICS)=SQRT((UTPICS(ICS)-TPICS(ICS)*TPICS(ICS)))
          END DO
C!10-MAR-2025 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C!        Normalize per radial bin volume
          DO IRB=1,NRBINS
            FNORM=ONE/VOLRB(IRB)
C           Normalize
            EIMP(IRB)=EIMP(IRB)*FNORM
            UEIMP(IRB)=UEIMP(IRB)*FNORM
            CION(IRB)=CION(IRB)*FNORM
            UION(IRB)=UION(IRB)*FNORM
            DO ICS=1,MAXICS
              PICS(IRB,ICS)=PICS(IRB,ICS)*FNORM
              UICS(IRB,ICS)=UICS(IRB,ICS)*FNORM
            END DO
          END DO
          EXIT Results
        END DO Results

C!      Save results to output file ###############################      
        Output: DO ! Dummy Loop for output
          NOTIME=0
          IF(IDEBUG.LT.0) NOTIME=1
          OPEN(11,FILE=DATEXT(PREFIX,FILENM,NOTIME),STATUS='UNKNOWN')
C!        Info on program and user
          WRITE(11,'(a)') '# Output from AVC_T62U2_V6.f - Version '//
     &             VDATE
          WRITE(11,'(a)') '# produced by '//TRIM(USERID)//
     &                    ' on '//TSTAMP()
          WRITE(11,'(a)') '# ******************************************'  
C!        Info on processed data 
          WRITE(11,'(a)') '# Processed input data file: '//
     &                    TRIM(FILENM)//' | '//TRIM(INTSTR(NEV))//
     &                    ' sequences of equal eventID.'
          WRITE(11,'(a)') '# Parameters:'//
     &                    ' | NITER = '//TRIM(INTSTR(NITER))//
     &                    ' | ICSMAX = '//TRIM(INTSTR(ICSMAX))//
     &                    ' | RMIN = '//TRIM(FLTSTR(RMIN,5))//
     &                    ' | RMAX = '//TRIM(FLTSTR(RMAX,5))//
     &                    ' | DLOGR = '//TRIM(FLTSTR(DLOGR,4))//' |'
          WRITE(11,'(a)') '# Results written to file: '//
     &                    DATEXT(PREFIX,FILENM,NOTIME)     
          WRITE(11,'(a)') '# ******************************************'    
C!        Global scores 
          WRITE(11,'(a,f8.3,a)') '# Total energy imparted per track: '
     &                         //TRIM(FLTSTR(TEIMP,2))//' eV'
     &                         //'  |  Total number of ionizations: '
     &                         //TRIM(FLTSTR(TCION,2)) 
     &                         //'  |  Average energy per ionization: '
     &                         //TRIM(FLTSTR(TEIMP/TCION,2))//' eV'
C! 10-MAR-2025 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          WRITE(11,'(a,f8.3,a)') '# Corresponding uncertainties:     '
     &                         //TRIM(FLTSTR(UTEIMP,2))//' eV'
     &                         //'  |                               '
     &                         //TRIM(FLTSTR(UTCION,2)) 
C!     &                         //'  |  Average energy per ionization: '
C!     &                         //TRIM(FLTSTR(TEIMP/TCION,2))//' eV'

*         IF(ICSMAX.LT.MAXICS) THEN
*           UTPICS(ICSMAX)=UTPICS(ICSMAX)**2
*           DO ICS=ICSMAX+1,MAXICS
*             TPICS(ICSMAX)=TPICS(ICSMAX)+TPICS(ICS)
*             UTPICS(ICSMAX)=UTPICS(ICSMAX)+UTPICS(ICS)**2
*           END DO
*           UTPICS(ICSMAX)=SQRT(UTPICS(ICSMAX))
*         ENDIF
C!10-MAR-2025 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

          WRITE(11,'(a)') '# Total number of cluster sites: '
          WRITE(11,'(a,50i11)')   '# n_ion=      ',(ICS,ICS=1,ICSMAX)
          WRITE(11,'(a,$)')       '# N(n_ion)=   '
          DO ICS=1,ICSMAX
            IF(TPICS(ICS).GE.10.) THEN
              WRITE(11,'(f11.2,$)') TPICS(ICS)
            ELSE IF(TPICS(ICS).GE.ONE) THEN
              WRITE(11,'(f11.3,$)') TPICS(ICS)
            ELSE IF(TPICS(ICS).GT.0.1) THEN
              WRITE(11,'(f11.4,$)') TPICS(ICS)
            ELSE
              WRITE(11,'(e11.4,$)') TPICS(ICS)
            ENDIF
          END DO
          
C! 10-MAR-2025 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          WRITE(11,*) ! 16-JUN-2025
          WRITE(11,'(a,$)') '# u[N(n_ion)]='
          DO ICS=1,ICSMAX
            IF(TPICS(ICS).GE.10.) THEN
              WRITE(11,'(f11.2,$)') UTPICS(ICS)
            ELSE IF(TPICS(ICS).GE.ONE) THEN
              WRITE(11,'(f11.3,$)') UTPICS(ICS)
            ELSE IF(TPICS(ICS).GT.0.1) THEN
              WRITE(11,'(f11.4,$)') UTPICS(ICS)
            ELSE
              WRITE(11,'(e11.4,$)') UTPICS(ICS)
            ENDIF
          END DO
C!10-MAR-2025 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          
          WRITE(11,*)
          WRITE(11,'(a)') '# ******************************************'    
C!        Radial densities and cluster probabilities
          WRITE(11,'(a,f8.3,a)') '# Radial densities per track of '
     &        //'energy imparted, ionizations, and ionization clusters '
     &        //' in electron tracks for DSITE=',DSITE,' nm '   
          WRITE(11,'(a)') '# column 1: outer shell radius in nm'
          WRITE(11,'(a)') '# column 2: volume of radial shell in nm^3'
          WRITE(11,'(a)') '# column 3: mean volume density of energy '
     &                    //'imparted (in eV/nm^3)'
          WRITE(11,'(a)') '# column 4: mean volume density of the '
     &                    //'number of ionizations (in nm^-3) '
          WRITE(11,'(a)') '# columns 5...'//TRIM(INTSTR(ICSMAX+3))
     &                    //': mean volume density (in nm^-3) of the '
     &                    //'number of sites with 1...'
     &                    //TRIM(INTSTR(ICSMAX-1))//' ionizations'
          WRITE(11,'(a)') '# column '//TRIM(INTSTR(4+ICSMAX))//': '
     &                    //'mean volume density (in nm^-3) of sites '
     &                    //'with '//TRIM(INTSTR(ICSMAX))//' or more '
     &                    //'ionizations'
C!        Uncertainties of radial densities and cluster probabilities
          WRITE(11,'(a)') '# column '//TRIM(INTSTR(5+ICSMAX))//': '
     &                    //'standard deviation of the volume density '
     &                    //'of energy imparted (in eV/nm^3)'
          WRITE(11,'(a)') '# column '//TRIM(INTSTR(6+ICSMAX))//': '
     &                    //'standard deviation of the volume density '
     &                    //'of the number of ionizations (in nm^-3)'
          WRITE(11,'(a)') '# columns '//TRIM(INTSTR(7+ICSMAX))//'...'
     &                    //TRIM(INTSTR(5+2*ICSMAX))//': standard '
     &                    //'deviation of the volume density (in nm^-3)'
     &                    //' of the number of sites with 1...'
     &                    //TRIM(INTSTR(ICSMAX-1))//' ionizations'
          WRITE(11,'(a)') '# column '//TRIM(INTSTR(6+2*ICSMAX))//': '
     &                    //'Standard deviation of the volume '
     &                    //'density (in nm^-3) of sites '
     &                    //'with '//TRIM(INTSTR(ICSMAX))//' or more '
     &                    //'ionizations'
          WRITE(11,'(a)') '# ******************************************' 
                    
          WRITE(11,'(100a12)') 'Rbin','Vbin','dEimp/dV','dNi/dV',
     &                         (ADJUSTR(COLHDR(ICS,1)),ICS=1,ICSMAX),  
     &                         'u(dEimp/dV)','u(dNi/dV)',
     &                         (ADJUSTR(COLHDR(ICS,2)),ICS=1,ICSMAX)
          

          DO IRB=1,NRBINS                      !09-MAR-2025
            IF(PICS(IRB,1).GT.ZERO) MRBINS=IRB !09-MAR-2025
          END DO                               !09-MAR-2025
          MRBINS=MRBINS+1                      !09-MAR-2025
          
          DO IRB=1,MRBINS !09-MAR-2025
C!09-MAR-2025          DO IRB=1,MRBINS
            DO ICS=1,ICSMAX
              RESU(ICS,1)=PICS(IRB,ICS)
              RESU(ICS,2)=UICS(IRB,ICS)
            END DO
C!10-MAR-2025 Corrected 16-JUN-2025 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            IF(ICSMAX.LT.MAXICS) THEN
              RESU(ICSMAX,2)=RESU(ICSMAX,2)**2
              DO ICS=ICSMAX+1,MAXICS
                RESU(ICSMAX,1)=RESU(ICSMAX,1)+PICS(IRB,ICS)
                RESU(ICSMAX,2)=RESU(ICSMAX,2)+UICS(IRB,ICS)**2
              END DO
            END IF
            RESU(ICSMAX,2)=SQRT(RESU(ICSMAX,2))
C!10-MAR-2025 Corrected 16-JUN-2025 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

C!*            PRINT*,RBIN(IRB),VOLRB(IRB),EIMP(IRB),
C!*     &          UEIMP(IRB),CION(IRB),UION(IRB),J
            WRITE(11,'(f12.3,100e12.4)') RBIN(IRB),VOLRB(IRB),
     &          EIMP(IRB),CION(IRB),(RESU(ICS,1),ICS=1,ICSMAX),
     &          UEIMP(IRB),UION(IRB),(RESU(ICS,2),ICS=1,ICSMAX)            
          END DO
          CLOSE(11)
          WRITE(*,'(2a,/)') 'Results saved to ',
     &                      DATEXT(PREFIX,FILENM,NOTIME)


C! 17-JUN-2025 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
          OPEN(11,FILE='TRT_'//DATEXT(PREFIX,FILENM,NOTIME), 
     &            STATUS='UNKNOWN')
C!        Info on program and user
          WRITE(11,'(a)') '# Output from AVC_T62U2_V6.f - Version '//
     &             VDATE
          WRITE(11,'(a)') '# produced by '//TRIM(USERID)//
     &                    ' on '//TSTAMP()
          WRITE(11,'(a)') '# ******************************************'  
C!        Info on processed data 
          WRITE(11,'(a)') '# Processed input data file: '//
     &                    TRIM(FILENM)//' | '//TRIM(INTSTR(NEV))//
     &                    ' sequences of equal eventID.'
          WRITE(11,'(a)') '# Parameters:'//
     &                    ' | NITER = '//TRIM(INTSTR(NITER))//
     &                    ' | ICSMAX = '//TRIM(INTSTR(ICSMAX))//
     &                    ' | RMIN = '//TRIM(FLTSTR(RMIN,5))//
     &                    ' | RMAX = '//TRIM(FLTSTR(RMAX,5))//
     &                    ' | DLOGR = '//TRIM(FLTSTR(DLOGR,4))//' |'
          WRITE(11,'(a)') '# Results written to file: '//
     &                    'TRT_'//DATEXT(PREFIX,FILENM,NOTIME)     
          WRITE(11,'(a)') '# ******************************************'  
          IF(NEV.GT.NTRMAX) NEV = NTRMAX 
          DO ITR=1,NEV
            DO ICS=1,ICSMAX
              RESU(ICS,1)=TTICS(ITR,ICS)
            END DO
            DO ICS=ICSMAX+1,MAXICS
              RESU(ICSMAX,1)=RESU(ICSMAX,1)+TTICS(ITR,ICS)
            END DO
            WRITE(11,'(i7,100e12.4)') ITR,TTEIMP(ITR),TTIONS(ITR),
     &                                    (RESU(ICS,1),ICS=1,ICSMAX)         
          END DO
          CLOSE(11)
C! 17-JUN-2025 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          
          EXIT Output
        END DO Output
C!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      END DO Main_Loop ! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
      
      END
C!  *******************************************************************

C!  *******************************************************************
      SUBROUTINE DO_AVC(NPTS,X,Y,Z,E,ISION,TRICS,TREIMP,TRION) ! 28-JAN-2024
C! 28-JAN-2024      SUBROUTINE DO_AVC(NPTS,X,Y,Z,E,AVC,TEIM,TION)
C!  -------------------------------------------------------------------
      IMPLICIT NONE
C!    ----- Input Parameters: 
      INTEGER*4 NPMAX
      PARAMETER (NPMAX=100000)
C!    ----- Input variables: 
      INTEGER*4 NPTS,ISION(NPMAX)
      REAL*4 X(NPMAX),Y(NPMAX),Z(NPMAX),E(NPMAX)
C!    ----- Functions: 
      INTEGER*4 IRBIN
      REAL*4 RAND, RELVOL
C!    ----- Global parameters: 
      INTEGER*4 MAXICS,MAXNRB
      PARAMETER (MAXICS=20,MAXNRB=1001)
C!    ----- Global variables: 
C!    Constants
      INTEGER*4 ICSMAX,NITER
      LOGICAL LSHARE
      COMMON /CONSTS/ ICSMAX,LSHARE,NITER
      DATA NITER / 1 /
C!    Geometry
      INTEGER*4 NRBINS ! Actual number of radial bins
      REAL*4 RRBIN2(MAXNRB),RRBHI2(MAXNRB),RRBLO2(MAXNRB),RSITE,RSITE2
      COMMON /GEOMET/ RRBIN2,RRBHI2,RRBLO2,RSITE,RSITE2,NRBINS
      DATA RSITE2,NITER / 2.25, 1 /
C!    ----- Local parameters: 
      REAL*4 ONE,TWOPI,ZERO
      PARAMETER (ONE=1.,TWOPI=6.2831853,ZERO=0.)
C!    ----- Local arrays: 
      REAL*4 TRICS(MAXNRB,MAXICS),TREIMP(MAXNRB),TRION(MAXNRB)
C!    ----- Local scalars: 
      INTEGER*4 I,ICS,IRB,IRBHI,IRBLO,IREP,J
      REAL*4 R,CTH,STH,PHI,DX,DY,DZ,RRC2,WEIGHT,XC,YC,ZC
C
C!    End of declarations section
C!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      

      DO IRB=1,NRBINS
        TRION(IRB)=ZERO
        TREIMP(IRB)=ZERO
        DO ICS=1,MAXICS
          TRICS(IRB,ICS)=ZERO  
        END DO
      END DO
      
      DO I=1,NPTS ! 28-JAN-2024 HR: Reversed sequence of loops
C       Score energy imparted and ionization
        RRC2=(X(I)*X(I)+Y(I)*Y(I)+Z(I)*Z(I))/RSITE2
        IRB=IRBIN(RRC2,0)
        TREIMP(IRB)=TREIMP(IRB)+E(I)
        
        IF(MOD(ISION(I),2).NE.1) CYCLE ! 28-JAN-2024
        TRION(IRB)=TRION(IRB)+1.
          
        DO IREP=1,NITER ! 28-JAN-2024 HR: Reversed sequence of loops 
C         coordinates of centers
          R=RSITE*EXP(ALOG(RAND())/3.) ! 11-NOV-2023 BUG fix: replaced RSITE2 with RSITE
          CTH=2.*RAND()-1. 
          STH=SQRT(1.-CTH*CTH)
          PHI=TWOPI*RAND() 
          XC=X(I)+R*STH*COS(PHI)
          YC=Y(I)+R*STH*SIN(PHI)
          ZC=Z(I)+R*CTH
          RRC2=(XC*XC+YC*YC+ZC*ZC)/RSITE2
          
C         Ionization cluster size
          ICS=0
          DO J=1,NPTS 
            DX=X(J)-XC
            DY=Y(J)-YC
            DZ=Z(J)-ZC
            IF(DX*DX+DY*DY+DZ*DZ.LE.RSITE2) ICS=ICS+1
          END DO

C         Score radial bins
          IF(LSHARE) THEN ! share between bin that contain part of the site volume
            IRBLO=IRBIN(RRC2,-1)
            IRBHI=IRBIN(RRC2,+1)
            DO IRB=IRBLO,IRBHI
              WEIGHT=RELVOL(RRC2,IRB)
              IF(ICS.GE.MAXICS) THEN
                TRICS(IRB,MAXICS)=TRICS(IRB,MAXICS)+WEIGHT/ICS ! Division by ICS must be done here since info is unavailable afterwards
              ELSE
                TRICS(IRB,ICS)=TRICS(IRB,ICS)+WEIGHT  ! Division by ICS can and will be done at the end to enhance calculation efficiency
              ENDIF
            END DO
          ELSE ! The cluster is only counted in the bin where its center is located
            IRB=IRBIN(RRC2,0)
            WEIGHT=1.0
            IF(ICS.GE.MAXICS) THEN
              TRICS(IRB,MAXICS)=TRICS(IRB,MAXICS)+WEIGHT/ICS ! Division by ICS must be done here since info is unavailable afterwards
            ELSE
              TRICS(IRB,ICS)=TRICS(IRB,ICS)+WEIGHT  ! Division by ICS can and will be done at the end to enhance calculation efficiency
            ENDIF
          ENDIF
        END DO
      END DO
C
      DO IRB=1,NRBINS
        DO ICS=1,MAXICS-1                    ! exclude the max. ICS, for which F_k is scored and the correction was already applied
          TRICS(IRB,ICS)=TRICS(IRB,ICS)/(NITER*ICS)  ! normalize to number of iterations and cluster size to remove multiple counts
        END DO
        ICS=MAXICS
        TRICS(IRB,ICS)=TRICS(IRB,ICS)/NITER          ! for the max. ICS F_k is scored and the correction was already applied
C! 28-JAN-2024 Next two lines are obsolete due to change of loops further up
C! 28-JAN-2024        TRION(IRB)=TRION(IRB)/NITER ! 09-NOV-2023 BUG fix; this normalization was missing before
C! 28-JAN-2024        TREIMP(IRB)=TREIMP(IRB)/NITER ! 09-NOV-2023 BUG fix; this normalization was missing before
      END DO
      
      RETURN
      END
C!  *******************************************************************

C!  *******************************************************************
      REAL*4 FUNCTION RELVOL(XR2,IRB)
C!  -------------------------------------------------------------------
      INTEGER*4 IRB ! Index of radial bin
      REAL*4 XR2 ! square of radial position of site center in units of site radius
      
      REAL*4 ONE,TWO,ZERO
      PARAMETER (MAXNRB=1001,ONE=1.0,TWO=2.0,ZERO=0.0)
      REAL*4 RRBIN2(MAXNRB),RRBHI2(MAXNRB),RRBLO2(MAXNRB),RSITE,RSITE2
      COMMON /GEOMET/ RRBIN2,RRBHI2,RRBLO2,RSITE,RSITE2,NRBINS
      
      REAL*4 RR, RR2, XR, F1, F2, F3, F4,WLOW
      !DATA F1,F2,F3,F4 / -0.1875000, 0.5000000,-0.3750000, 0.0625000/
      DATA F1,F2,F3,F4 / -0.1875000,-0.2500000,0.7500000, 0.5000000/
      
      RR2=RRBIN2(IRB) ! square of outer radius of bin in units of site radius
      XR=SQRT(XR2)
      YY=SQRT(RR2)-XR
      YY2=YY*YY
      
      IF(ABS(YY).LT.ONE) THEN
        RELVOL= F1*((YY2-TWO)*YY2+ONE)/XR + (F2*YY2+F3)*YY+F4
      ELSE
        IF(YY.GE.ONE) THEN
          RELVOL=ONE
        ELSE
          RELVOL=ZERO
        ENDIF
      ENDIF
      
      IF(IRB.GT.1) THEN
        RR2=RRBIN2(IRB-1)
        YY=SQRT(RR2)-XR
        YY2=YY*YY
        IF(ABS(YY).LT.ONE) THEN
          WLOW=F1*((YY2-TWO)*YY2+ONE)/XR + (F2*YY2+F3)*YY+F4
          RELVOL=RELVOL-WLOW
        ELSE
          IF(YY.GE.ONE) THEN
            RELVOL=RELVOL-ONE
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END
C!  *******************************************************************
      
C!  *******************************************************************
      INTEGER*4 FUNCTION IRBIN(RRC2,IFLAG)
C!  -------------------------------------------------------------------
      IMPLICIT NONE
C!    ----- Input variables: 
      INTEGER*4 IFLAG
      REAL*4 RRC2
C!    ----- Global parameters: 
      INTEGER*4 MAXNRB
      PARAMETER (MAXNRB=1001)
C!    ----- Global variables: 
      INTEGER*4 NRBINS
      REAL*4 RRBIN2(MAXNRB),RRBHI2(MAXNRB),RRBLO2(MAXNRB),RSITE,RSITE2
      COMMON /GEOMET/ RRBIN2,RRBHI2,RRBLO2,RSITE,RSITE2,NRBINS
C
C!    End of declarations section
C!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
C
      IRBIN=1
      IF(IFLAG.EQ.0) THEN
        DO WHILE(IRBIN.LT.NRBINS.AND.RRC2.GT.RRBIN2(IRBIN))
          IRBIN=IRBIN+1
        END DO
      END IF
      IF(IFLAG.EQ.-1) THEN
        DO WHILE(IRBIN.LT.NRBINS.AND.RRC2.GT.RRBHI2(IRBIN))
          IRBIN=IRBIN+1
        END DO
      END IF
      IF(IFLAG.EQ.1) THEN
        DO WHILE(IRBIN.LT.NRBINS.AND.RRC2.GT.RRBLO2(IRBIN))
          IRBIN=IRBIN+1
        END DO
      END IF
      
      RETURN
      END
C!  *******************************************************************

C!  *******************************************************************
      SUBROUTINE RFINIT()
C!  -------------------------------------------------------------------
C!    Sets the information for the columns of the input file
C!      
C!    ------ Global variables
C!    ------
C!     CODE: Character string encoding the meaning of the entries in a 
C!           line of the input file as follows:
C!           'T'     - number of the primary particle track
C!           'X','Y','Z' - x, y, and z coordinates of the transfer point 
C!           'E'     - energy deposit (if applicable)
C!           'I'     - ionization cluster size (if applicable)
C!           'C'     - GEANT4 code for interaction type
C!           '%'     - additional data that are not used 
C!           Example: The string 'TZ%%EXY ' indicates that there are 7
C!                    entries in each line, of which the first is the 
C!                    track number, the second the z coordinate, the 
C!                    fifth the energy, and the sixth and seventh the x
C!                    and y coordinates
C!     IDT: Array of the columns indices of (1-3) x,y,z coordinates 
C!                   (4) energy deposit (if present), (5) track number, 
C!                   (6) number of ionizations in cluster (if present)
C!                   (7-8) are there for future use
      CHARACTER CODE*8
      INTEGER*2 IDT(8),NDT
      COMMON /RFORMT/IDT,NDT,CODE
C!    Local variables
      INTEGER*1 J
      
      DO J=1,8
        IDT(J)=0
      END DO
      
      DO J=1,8
        IF(CODE(J:J).EQ.'X'.AND.IDT(1).EQ.0) IDT(1)=J
        IF(CODE(J:J).EQ.'Y'.AND.IDT(2).EQ.0) IDT(2)=J
        IF(CODE(J:J).EQ.'Z'.AND.IDT(3).EQ.0) IDT(3)=J
        IF(CODE(J:J).EQ.'E'.AND.IDT(4).EQ.0) IDT(4)=J
        IF(CODE(J:J).EQ.'T'.AND.IDT(5).EQ.0) IDT(5)=J
        IF(CODE(J:J).EQ.'I'.AND.IDT(6).EQ.0) IDT(6)=J
        IF(CODE(J:J).EQ.'C'.AND.IDT(7).EQ.0) IDT(7)=J
        IF(CODE(J:J).NE.' ') NDT=J
      END DO
      
      IF(IDT(1).EQ.0) PRINT*, 'No data column for X'
      IF(IDT(2).EQ.0) PRINT*, 'No data column for Y'
      IF(IDT(3).EQ.0) PRINT*, 'No data column for Z'
      IF(IDT(5).EQ.0) PRINT*, 'No data column for track number'
      IF(IDT(1)*IDT(2)*IDT(3)*IDT(5).EQ.0) STOP
      
      END
C!  *******************************************************************

C!  *******************************************************************
      CHARACTER*16 FUNCTION FLTSTR(X,NDIG)
C!  -------------------------------------------------------------------
      INTEGER*4 NDIG
      REAL*4 X
      CHARACTER FMZ*7
      
      FMZ='(F16.0)'
      WRITE(FMZ(6:6),'(i1)') NDIG
      FLTSTR=''
      WRITE(FLTSTR,FMZ) X
      N=16
      DO WHILE(FLTSTR(N:N).EQ.'0')
        FLTSTR(N:N)=' '
        N=N-1
      END DO
      FLTSTR=ADJUSTL(FLTSTR)
      RETURN
      END
C!  *******************************************************************

C!  *******************************************************************
      CHARACTER*8 FUNCTION INTSTR(N)
C!  -------------------------------------------------------------------
      INTEGER*4 N
      INTSTR=''
      WRITE(INTSTR,'(I8)') N
      INTSTR=ADJUSTL(INTSTR)
      RETURN
      END
C!  *******************************************************************

C!  *******************************************************************
      CHARACTER*24 FUNCTION TSTAMP()
C!  -------------------------------------------------------------------
      CHARACTER DATE*8, TIME*10, ZONE*5, TIMEST*24
      INTEGER*4 IDT(8)
      CALL DATE_AND_TIME(DATE,TIME,ZONE,IDT)
C!              123456789 123456789 1234
      TIMEST = 'DD-MMM-YYYY HH:MM +HH:MM'
      TSTAMP=DATE(7:8)//'-'//DATE(5:6)//'-'//DATE(1:4)//' '//
     &       TIME(1:2)//':'//TIME(3:4)//' '//ZONE(1:3)//':'//ZONE(4:5)
      RETURN
      END
C!  *******************************************************************

C!  *******************************************************************
      CHARACTER*80 FUNCTION DATEXT(PREFIX,FILENM,NOTIME)
C!  -------------------------------------------------------------------
      CHARACTER DATE*8, FILENM*80,PREFIX*40,TIME*10,ZONE*5
      INTEGER*4 IDT(8)
      INTEGER IPOS, I, NOTIME
      LOGICAL HASEXT
C 
      IPOS=0
      DO I=2,LEN_TRIM(FILENM)
        IF(FILENM(I:I).EQ.'.') IPOS=I
      END DO 
C 
      CALL DATE_AND_TIME(DATE,TIME,ZONE,IDT)
C!
      IF(IPOS.GT.0) THEN
        IF(NOTIME.GT.0) THEN
          DATEXT=FILENM(1:IPOS-1)//'_'//DATE//
     &           FILENM(IPOS:LEN_TRIM(FILENM))
        ELSE
          DATEXT=FILENM(1:IPOS-1)//'_'//DATE//'_'//TIME(1:4)//
     &           FILENM(IPOS:LEN_TRIM(FILENM))
        ENDIF
      ELSE
        IF(NOTIME.GT.0) THEN
          DATEXT=TRIM(FILENM)//'_'//DATE
        ELSE
          DATEXT=TRIM(FILENM)//'_'//DATE//'_'//TIME(1:4)
        ENDIF
      ENDIF
      DATEXT=TRIM(PREFIX)//DATEXT
      RETURN
      END
C!  *******************************************************************

C!  *******************************************************************
      SUBROUTINE SORT62(FILENM,IFLAG)
C!  -------------------------------------------------------------------
      IMPLICIT NONE
      
      CHARACTER*80 FILENM
      INTEGER*4 IFLAG
      
      INTEGER*4 MAXELM, MAXINT, MAXTRC
      PARAMETER(MAXELM=10000000,MAXTRC=10000,MAXINT=10000)
      
      INTEGER*4 IEVENT(MAXELM),ISTACK(MAXELM,2),
     &          LINDEX(MAXELM),LSTACK(MAXELM,2),
     &          IALIAS(MAXTRC),NMB(MAXTRC),IFIRST(MAXTRC),ILAST(MAXTRC)
      CHARACTER*9 XYZE(4)
      CHARACTER*50 LINES(MAXELM)
      CHARACTER*80 WRLINE
      CHARACTER*10 OUTEV*8   
      INTEGER*4 I,ICOPY,IEV,IFILE,ISTART,ITAKE,ITOOK,J,NALIAS,NGIVES, 
     &          NEVENT,NFILES,NLAST,NLINES,NTAKES
C!  -------------------------------------------------------------------

      DO I=1,MAXTRC
        IALIAS(I)=-1
      END DO

C!    Open filename 
      OPEN(11,FILE=FILENM,STATUS='OLD')
      WRITE(*,'(a,$)') 'Sorting '//TRIM(FILENM)
        
C!    Read data
      NALIAS=0       
      DO I=1,MAXELM
        READ(11,'(a)',END=10) LINES(I)
        READ(LINES(I),*) NEVENT
        NEVENT=NEVENT+1
        IF(IALIAS(NEVENT).EQ.-1) THEN
          NALIAS=NALIAS+1
          IALIAS(NEVENT)=NALIAS
          NMB(NALIAS)=NEVENT
          IFIRST(NALIAS)=I
          IF(MOD(NALIAS,100).EQ.0) WRITE(*,'(a,$)') '_'
        ENDIF
        IEVENT(I)=IALIAS(NEVENT)
        LINDEX(I)=I
        ILAST(IEVENT(I))=I
        NLINES=I
      END DO 
  10  CONTINUE
      CLOSE(11)
        
C!    Find first 'split' track
      ISTART=0
      DO ISTART=1,NALIAS-1
        IF(IFIRST(ISTART+1).LT.ILAST(ISTART)) EXIT
      END DO
      
      IF(ISTART.LT.NALIAS-1) THEN ! 11-NOV-2023 BUG fix; when the data are already sorted, do nothing, since else the data of the first track are overwritten with some of the last one
C!    1st pass 
        ITAKE=0
        ICOPY=IFIRST(ISTART+1)
        DO J=IFIRST(ISTART+1),ILAST(ISTART)
          IF(IEVENT(J).NE.ISTART) THEN ! other events go to stack
            ITAKE=ITAKE+1
            ISTACK(ITAKE,1)=IEVENT(J)
            LSTACK(ITAKE,1)=LINDEX(J)
          ELSE ! same event data are appended to those at the head
            LINDEX(ICOPY)=LINDEX(J)
            IEVENT(ICOPY)=IEVENT(J)
            ICOPY=ICOPY+1
          ENDIF
        END DO
        NLAST=ILAST(ISTART)
      
        NTAKES=2
        NGIVES=1
        DO IEV=ISTART+1,NALIAS
          IF(MOD(IEV,100).EQ.0) WRITE(*,'(a,$)') '.'
          ITOOK=ITAKE
          ITAKE=0
C!        First take elements from the stack
          DO J=1,ITOOK
            IF(ISTACK(J,NGIVES).NE.IEV) THEN ! other events go to stack
              ITAKE=ITAKE+1
              ISTACK(ITAKE,NTAKES)=ISTACK(J,NGIVES)
              LSTACK(ITAKE,NTAKES)=LSTACK(J,NGIVES)
            ELSE ! same event data are appended to those at the head
              LINDEX(ICOPY)=LSTACK(J,NGIVES)
              IEVENT(ICOPY)=ISTACK(J,NGIVES)
              ICOPY=ICOPY+1
            ENDIF
          END DO
C!        Then those that are in the original array          
          DO J=NLAST+1,ILAST(IEV)
            IF(IEVENT(J).NE.IEV) THEN ! other events go to stack
              ITAKE=ITAKE+1
              ISTACK(ITAKE,NTAKES)=IEVENT(J)
              LSTACK(ITAKE,NTAKES)=LINDEX(J)
            ELSE ! same event data are appended to those at the head
              LINDEX(ICOPY)=LINDEX(J)
              IEVENT(ICOPY)=IEVENT(J)
              ICOPY=ICOPY+1
            ENDIF
          END DO
C!        When done with this event, switch stacks
          NGIVES=3-NGIVES
          NTAKES=3-NTAKES
          IF(ILAST(IEV).GT.NLAST) NLAST=ILAST(IEV)
        END DO
      ENDIF ! 11-NOV-2023 end BUG fix
      
      IF(IFLAG.EQ.1) THEN
        OPEN(11,STATUS='SCRATCH')
      ELSE IF(IFLAG.EQ.2) THEN
        OPEN(11,FILE='S_'//FILENM,STATUS='UNKNOWN')
      ENDIF
      
      DO I=1,NLINES
        WRLINE=LINES(LINDEX(I))
        WRITE(11,'(a)') TRIM(WRLINE)
      END DO
      WRITE(*,'(a)') 'Done.'
      REWIND(11)
      
      END
C!  *******************************************************************

C! 28-JAN-2024 (HR) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
C!  *******************************************************************
      LOGICAL FUNCTION COLION(LUN,NHEADL)
      INTEGER*4 LUN,NHEADL
C!    Local variables
      CHARACTER LINE*80
      INTEGER*4 IEV,ITHERE
      REAL*4 X,Y,Z,E
C!    --------------------------------------------------------------
      COLION=.FALSE.  
      
C!    Read header lines and first data line
      DO I=1,NHEADL+1
        READ(LUN,'(a)') LINE
      END DO  
      
      READ(LINE,*,ERR=10,END=10) IEV,X,Y,Z,E,ITHERE
      COLION=.TRUE.
  10  CONTINUE

C!    Rewind file and read header lines
      REWIND(LUN)
      DO I=1,NHEADL
        READ(LUN,'(a)') LINE
      END DO
      
      RETURN
      END
C!  *******************************************************************
C! 28-JAN-2024 (HR) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<          