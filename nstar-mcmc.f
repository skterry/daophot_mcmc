      SUBROUTINE  NSTAR  (PAR, MAXPAR, PSF, MAXPSF, MAXEXP, DATA, 
     .     NCOL, NROW, MAXCOL, WATCH, FITRAD, E1, E2)
C
C=======================================================================
C
C Photometry for many stars by simultaneous multiple PSF fits.
C
C              OFFICIAL DAO VERSION:  1991 April 18
C
C Currently operates on a picture no larger than 640 X 1024 pixels 
C total, regardless of format, and no more than 60 stars at a time.  
C The latter restriction may be altered by changing the first 
C parameter.
C
C Arguments:
C
C  FITRAD (INPUT) is the fitting radius specified as a user-definable 
C         option.  It governs how many pixels out from the centroid of
C         the star will actually be considered in computing the least-
C         squares profile fits.
C
C   WATCH (INPUT) is the 'watch progress' parameter specified by the
C         user.  If WATCH > 0, information relating to the progress of
C         the reductions will be typed on the terminal during execution.
C
C=======================================================================
C
      IMPLICIT NONE
      INTEGER MAXSTR, MAXEXP, MAXPSF, MAXPAR, NCOL, NROW, MAXCOL
      PARAMETER  (MAXSTR=60)
C
C Parameters:
C
C MAXSTR The maximum number of stars in a single group.  This parameter
C        is determined primarily by the execution time per iteration--
C        at MAXSTR=60, our VAX 11/780 takes around 2.5 CPU minutes per 
C        iteration.  For MAXSTR > 150 or so, the accuracy of inverting 
C        the REAL*4 design matrix would also begin to suffer.
C
C MAXPSF the largest PSF look-up table that can be accomodated.  If

C        MAXPSF = 2*[2*(MAXRAD+1)]+7.
C
      REAL DATA(MAXCOL,*), PSF(MAXPSF,MAXPSF,MAXEXP)
      REAL PAR(MAXPAR)
      REAL C(3*MAXSTR+1,3*MAXSTR+1), V(3*MAXSTR+1)
      REAL XC(MAXSTR+1), YC(MAXSTR+1), MAG(MAXSTR+1), RPIXSQ(MAXSTR)
      REAL SKY(MAXSTR+1), CHI(MAXSTR+1), SUMWT(MAXSTR+1)
      REAL NUMER(MAXSTR+1), DENOM(MAXSTR+1), SHARP(MAXSTR+1)
      REAL MAGERR(MAXSTR+1)
      REAL XERR(MAXSTR+1)
      REAL YERR(MAXSTR+1)
      REAL CLAMP(3*MAXSTR+1), XOLD(3*MAXSTR+1), X(3*MAXSTR+1)
      INTEGER ID(MAXSTR+1), NPIX(MAXSTR), RDPSF
      LOGICAL SKIP(MAXSTR)
C
      REAL AMIN1, AMAX1, ABS, SQRT, USEPSF, PROFIL
      INTEGER MIN0, MAX0
C
      CHARACTER LINE*80
      CHARACTER COOFIL*40, MAGFIL*40, PSFFIL*40, FITFIL*40, GRPFIL*40, 
     .     MCMCFIL*80, CHISQPIXFIL*80, BESTFIL*80, SWITCH*40, EXTEND*40,
     .     MCMCALL*80
      CHARACTER CASE*4


      REAL IX10, IY10, IX20, IY20, IF0
      REAL IX30, IY30, IFB0, IFB, FB
      REAL IX10_IN, IY10_IN, IX20_IN, IY20_IN, IF0_IN
      REAL IX30_IN, IY30_IN, IFB0_IN
      REAL X1, Y1, X1MIN, Y1MIN
      REAL X2, Y2, X2MIN, Y2MIN
      REAL X3, Y3, X3MIN, Y3MIN
      REAL FMIN, FBMIN, EMIN, ZMIN, SSEPMIN
      REAL IX1M, IY1M, IX2M, IY2M, EEM_CHI2
      REAL IX3M, IY3M, IFBM, LSSEPMIN, BSEPMIN
      REAL IFM, FIT_STARS, NLINES
      REAL, ALLOCATABLE :: PPU_MIN(:), EMIN_CHI2(:), EMINPIX(:)
      REAL PROB, PROB_RAND
      REAL*8 EE0, EE0_NORM, RFAC, PSTEP, FSTEP
      REAL*8 PTOT, FTOT, PIXSCALE, ACCEPT_RAT, ACC_ARRAY
      REAL*8 Z0, LSSEP, BSEP
      REAL*8 CHISQ, CHISQ_MIN
      REAL*8 F, F1, F2
      REAL SSEP, SSEP2, ZM, EEM
      REAL, ALLOCATABLE :: FU(:), FU_MIN(:)
      REAL, ALLOCATABLE :: SGU(:), SGL(:), SGU2(:), PPU(:)
      INTEGER IF2, Steps, IT, GRIDSIZE
      INTEGER U, US, Counter
      INTEGER NIMU, NIT
      INTEGER NSTU
      INTEGER PPOS
      REAL, DIMENSION(9) :: RNARRAY
      CHARACTER(LEN=15) A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14
      CHARACTER(LEN=15) A15,A16,A17,A18,A19,A20,A21,A22,A23,A24

      REAL LOBAD, DF, DX, DY, ERR, PSFMAG, BRIGHT, XPSF, YPSF
      REAL SEPCRIT, PSFRAD, RADIUS, THRESH, AP1, PHPADU, RONOIS
      REAL DUM, CUTOFF, RADSQ, SUMSKY, CHIGRP, HIBAD, SEPMIN
      REAL PSFRSQ, SKYBAR, WCRIT, XMIN, XMAX, YMIN, YMAX, SEP
      REAL SUMGRP, GRPWT, D, WT, DELTAX, DELTAY, VAL, DVDXC, DVDYC
      REAL RSQ, DPOS, RELERR, SIGSQ, RHOSQ, DWT, WATCH, FITRAD, E1, E2
      REAL PERERR, PKERR, DFDSIG, FAINT, RDNOIS
      INTEGER I, J, K, IPSTYP, NPSF, NPAR, NEXP, NFRAC, NL, NTOT, NSTR
      INTEGER IX, IY, I3, I3M2, LX, LY, L, NI, IFAINT
      INTEGER NITER, NTERM, IER, IDUM, IXMIN, IXMAX, IYMIN, IYMAX
      INTEGER MAXUNK
      LOGICAL OMIT, REDO, CLIP
      COMMON /FILNAM/ COOFIL, MAGFIL, PSFFIL, FITFIL, GRPFIL, MCMCFIL
      COMMON /FILNAM/ CHISQPIXFIL, BESTFIL, MCMCALL

C
C-----------------------------------------------------------------------
C
C SECTION 1
C
C Get ready, get set, . . .
C
      MAXUNK=MAXSTR*3+1
      CALL TBLANK
C
C Read the entire picture into memory.
C
      LX=1
      LY=1
      CALL RDSECT (1, LY, NROW, DATA, MAXCOL, IER)
      IF (IER .NE. 0) GO TO 9100
C
C Read the point-spread function into memory.
C
  940 CALL GETNAM ('File with the psf:', PSFFIL)
      IF ((PSFFIL .EQ. 'END-OF-FILE') .OR.
     .     (PSFFIL .EQ. 'GIVE UP')) THEN
         PSFFIL = ' '
         RETURN
      END IF
C
      PSFFIL = EXTEND(PSFFIL, CASE('psf'))
      IER = RDPSF (PSFFIL, IPSTYP, PAR, MAXPAR, NPAR,
     .     PSF, MAXPSF, MAXEXP, NPSF, NEXP, NFRAC,
     .     PSFMAG, BRIGHT, XPSF, YPSF)
      IF (IER .NE. 0) THEN
         PSFFIL = 'GIVE UP'
         GO TO 940
      END IF
      PERERR = 0.01 * E1
      PKERR = 0.01 * E2 / (PAR(1)*PAR(2))
C
C Stars will be checked for merger if they are separated by less than
C 1 FWHM of the image core.
C
C     Crit. sep. = 2.355*sigma, where
C          sigma = SQRT [ (sigma(X)**2 + sigma(Y)**2)/2 ]
C
C Use the quadratic mean of PAR(1) and PAR(2) as the HWHM.
C
      SEPCRIT=2.*(PAR(1)**2+PAR(2)**2)
C
C SEPCRIT contains the square of the critical separation.
C
      SEPMIN = 0.14*SEPCRIT
C
      PSFRAD = (REAL(NPSF-1)/2. - 1.)/2.
      PSFRSQ = PSFRAD**2
      RADIUS = AMIN1(FITRAD, PSFRAD)
C
C Ascertain the name of the file with the stellar groups, and open it.
C
      CALL GETNAM ('File with stellar groups:', GRPFIL)
      IF (GRPFIL .EQ. 'END-OF-FILE') THEN
         RETURN
      END IF
C
      CALL INFILE (2, GRPFIL, IER)
      IF (IER .LT. 0) GO TO 9300
C
      CALL RDHEAD (2, NL, IDUM, IDUM, LOBAD, HIBAD, THRESH, AP1, 
     .     PHPADU, RONOIS, DUM)
      IF (NL .NE. 3) GO TO 9200
C
C Inquire the name of the output file, and open it.
C
      FITFIL=SWITCH(GRPFIL, CASE('.nst'))
  980 CALL GETNAM ('File for NSTAR results:', FITFIL)
      IF ((FITFIL .EQ. 'END-OF-FILE') .OR.
     .     (FITFIL .EQ. 'GIVE UP')) THEN
         CALL CLFILE (2)
         FITFIL = ' '
         RETURN
      END IF
C
      CALL OUTFIL (1, FITFIL, IER)
      IF (IER .NE. 0) THEN
         CALL STUPID ('Error opening output file '//FITFIL)
         FITFIL = 'GIVE UP'
         GO TO 980
      END IF
C
C
C Inquire the name of the MCMC output file, and open it.
C      MCMCFIL=SWITCH(GRPFIL, CASE('.mc'))
   24 CALL GETNAM ('File for FULL MCMC chain:', MCMCALL)
      IF ((MCMCFIL .EQ. 'END-OF-FILE') .OR.
     .     (MCMCFIL .EQ. 'GIVE UP')) THEN
         CALL CLFILE (3)
         MCMCFIL = ' '
         RETURN
      END IF
C     
C      CALL OUTFIL (11, MCMCFIL, IER)
C      IF (IER .NE. 0) THEN
C         CALL STUPID ('Error opening output file '//MCMCFIL)
C         MCMCFIL = 'GIVE UP'
C         GO TO 25
C      END IF
C      CLOSE(25)
C
C
C Inquire the name of the MCMC output file, and open it.
C      MCMCFIL=SWITCH(GRPFIL, CASE('.mc'))
   25 CALL GETNAM ('File for ACCEPTED MCMC chains:', MCMCFIL)
      IF ((MCMCFIL .EQ. 'END-OF-FILE') .OR.
     .     (MCMCFIL .EQ. 'GIVE UP')) THEN
         CALL CLFILE (3)
         MCMCFIL = ' '
         RETURN
      END IF
C     
C      CALL OUTFIL (11, MCMCFIL, IER)
C      IF (IER .NE. 0) THEN
C         CALL STUPID ('Error opening output file '//MCMCFIL)
C         MCMCFIL = 'GIVE UP'
C         GO TO 25
C      END IF
C      CLOSE(25)
C
C
C Inquire the name of the chisq_pix output file, and open it.
C      CHISQPIXFIL=SWITCH(MCMCFIL, CASE('.ch2'))
   26 CALL GETNAM ('File for CHISQ_PIXEL results:', CHISQPIXFIL)
      IF ((CHISQPIXFIL .EQ. 'END-OF-FILE') .OR.
     .     (CHISQPIXFIL .EQ. 'GIVE UP')) THEN
         CALL CLFILE (4)
         CHISQPIXFIL = ' '
         RETURN
      END IF
C      
C      CALL OUTFIL (1, CHISQPIXFIL, IER)
C      IF (IER .NE. 0) THEN
C         CALL STUPID ('Error opening output file '//CHISQPIXFIL)
C         CHISQPIXFIL = 'GIVE UP'
C         GO TO 9300
C      END IF     
C
C
C Inquire the name of the BEST_FIT output file, and open it.
C      BESTFIL=SWITCH(CHISQPIXFIL, CASE('.bf'))
   27 CALL GETNAM ('File for BEST_FIT results:', BESTFIL)
      IF ((BESTFIL .EQ. 'END-OF-FILE') .OR.
     .     (BESTFIL .EQ. 'GIVE UP')) THEN
         CALL CLFILE (5)
         BESTFIL = ' '
         RETURN
      END IF
C      
C      CALL OUTFIL (1, BESTFIL, IER)
C      IF (IER .NE. 0) THEN
C         CALL STUPID ('Error opening output file '//BESTFIL)
C         BESTFIL = 'GIVE UP'
C         GO TO 9300
C      END IF
      
      CALL WRHEAD (1, 1, NCOL, NROW, 7, LOBAD, HIBAD, THRESH, AP1, 
     .     PHPADU, RONOIS, RADIUS)
C
C Get ready to go.
C
      RDNOIS=RONOIS**2
      CALL TBLANK
      IF (WATCH .GT. 0.5) THEN
         WRITE (6,610)
  610    FORMAT (/' It = number of iterations for current group',
     .          //' n* = number of stars in current group',
     .          //' N* = number of stars up through current group'
     .         //)
         CALL OVRWRT ('    It   n*   N*  ', 1)
      END IF
      RADSQ=RADIUS**2
      CUTOFF=0.999998*RADSQ
      NTOT=0
C
C-----------------------------------------------------------------------
C
C SECTION 2
C
C GO.
C
C Loop over stellar groups.
C
 2000 IF (WATCH .GT. 0.5) CALL TBLANK
      I=0
      SUMSKY=0.
C
C Read in the next group of stars.
C
 2010 I=I+1
      CALL RDSTAR (2, 3, ID(I), XC(I), YC(I), MAG(I), SKY(I))
      IF (ID(I) .LT. 0) GO TO 2100
      IF (ID(I) .EQ. 0) GO TO 2110
      IF (I .GT. MAXSTR) GO TO 2020

C A single sky brightness value, equal to the arithmetic mean of the
C skies determined for the individual stars, will be used for the
C group as a whole.  (In its present form NSTAR leaves this group sky
C value a constant.  It could in principle be determined as a fitting
C parameter along with the stellar positions and magnitudes.  My
C experiments along this line have been disappointing, but if you want
C to try it, include in the code those lines flagged with 'If sky
C is to be determined'.)
C
      SUMSKY=SUMSKY+SKY(I)
C
C Convert magnitude to brightness, scaled relative to the PSF.
C
      MAG(I)=10.**(-0.4*(MAG(I)))  !-PSFMAG)) Add const. = 21 to MAG(I)

C
C If PHOTOMETRY was unable to obtain a magnitude for the star (Mag. =
C 99.999), NSTAR will give it the old university attempt anyway.
C
      IF (MAG(I) .LE. 1.E-4) MAG(I)=0.01
      MAGERR(I)=0.0
      XERR(I)=0.0
      YERR(I)=0.0
      SHARP(I)=0.0
      GO TO 2010
C
C The group is too large.  Type out a message, keep reading until
C a blank line is encountered, and then go back and start a completely
C new group.
C
 2020 WRITE (6,620) 7, MAXSTR
  620 FORMAT (' Group with more than ', A1, I2, ' stars.')
 2030 CALL RDSTAR (2, 3, I, DUM, DUM, DUM, DUM)
      IF (I .LT. 0) GO TO 9000
      IF (I .NE. 0) GO TO 2030
      GO TO 2000
C
C Either a blank line or the EOF has been encountered.  If at least one
C real star has been read in since the last blank line, reduce the 
C group.  If it is a blank line and no star has been read in, go back 
C and read in a new group (in case in editing the file, the user has 
C clumsily left several blank lines in a row).  If it is the EOF and 
C no star has been read in, return.
C
 2100 IF (I .EQ. 1) GO TO 9000
 2110 IF (I .EQ. 1) GO TO 2000
      NSTR=I-1
      NTOT=NTOT+NSTR
      SKYBAR=SUMSKY/NSTR
C
C Start reducing the group.
C
      NTERM=3*NSTR
C
C If sky is to be determined: NTERM=NTERM+1
C
C Initialize accumulators and constraints on parameter corrections.
C
      CHIGRP=1.0
      NITER=0
      CLIP=.FALSE.
      DO 2120 I=1,NTERM
      XOLD(I)=0.0
 2120 CLAMP(I)=1.0
C
C Update information on screen.
C
      IF (WATCH .GT. 0.5) THEN
         WRITE (LINE,622) NITER, NSTR, NTOT
         CALL OVRWRT (LINE(1:15), 2)
      END IF
C      
C Begin to iterate solution here.
C
 2200 NITER=NITER+1
 2210 IF (WATCH .GT. 0.5) THEN
         WRITE (LINE,622) NITER, NSTR, NTOT
  622    FORMAT (3I8)
         CALL OVRWRT (LINE(1:15), 2)
      END IF
C
C Set up critical error for star rejection.
C
      WCRIT = 1.E15
      IF (NITER .GE. 5) WCRIT = 1.0
      IF (NITER .GE. 10) WCRIT = 0.66667
      IF (NITER .GE. 15) WCRIT = 0.5
C If sky is to be determined: X(NTERM)=-1.0
C
C If there is more than one star, check to see whether any two stars 
C have merged.  Meanwhile, determine the upper and lower limits in x 
C and y of a rectangular box containing the centroids of all stars, and
C initialize a couple of accumulators.
C
      XMIN=NCOL
      XMAX=1.
      YMIN=NROW
      YMAX=1.
C
      DO 2230 I=1,NSTR
      CHI(I)=0.
      SUMWT(I)=0.
      NUMER(I)=0.
      DENOM(I)=0.
      XMIN=AMIN1(XMIN, XC(I))
      XMAX=AMAX1(XMAX, XC(I))
      YMIN=AMIN1(YMIN, YC(I))
      YMAX=AMAX1(YMAX, YC(I))
      IF (NSTR .EQ. 1) GO TO 2230
C
      DO 2220 J=1,I-1
      SEP=(XC(I)-XC(J))**2+(YC(I)-YC(J))**2
      IF (SEP .GT. SEPCRIT) GO TO 2220

C
C Two stars are overlapping.  Identify the fainter of the two.
C
      K=J
      IF (MAG(I) .LT. MAG(J)) K=I
      IF ((SEP .LT. SEPMIN) .OR. 
     .     (MAGERR(K)/MAG(K) .GT. WCRIT)) GO TO 2240
 2220 CONTINUE
C
 2230 CONTINUE
C
C No two stars have merged.
C
      GO TO 2260
 2240 CONTINUE

C
C Now eliminate the fainter of the two.
C
      IF (MAG(I) .LT. MAG(J)) I=J
C
C The K-th star is now the fainter of the two, the I-th, the brighter.
C
      IF (WATCH .GT. 0.5) THEN
         WRITE (LINE,623) NITER, NSTR, NTOT, ID(K), ID(I)
  623    FORMAT (3I5, 5X, 'Star', I8, ' merged with star', I8,
     .     ', so it''s been deleted.')
         CALL OVRWRT (LINE(1:78), 3)
      END IF
C
C Now replace the centroid of the I-th star with the weighted mean of
C the most recent estimates of the centroids of the I-th and K-th
C stars, and the brightness of the I-th with the sum of the brightnesses
C of the I-th and K-th.
C
      XC(I)=XC(I)*MAG(I)+XC(K)*MAG(K)
      YC(I)=YC(I)*MAG(I)+YC(K)*MAG(K)
      MAG(I)=MAG(I)+MAG(K)
      XC(I)=XC(I)/MAG(I)
      YC(I)=YC(I)/MAG(I)
C
C Remove the K-th star from the group.
C
      CALL DAORMV (K, MAXSTR, NSTR, ID, XC, YC, MAG, SKY)
      NTOT=NTOT-1
      IF (WATCH .GT. 0.5) THEN
         WRITE (LINE,622) NITER-1, NSTR, NTOT
         CALL OVRWRT (LINE(1:15), 2)
      END IF
      NTERM=NSTR*3
C     If sky is to be determined: NTERM=NTERM+1
C
C After deleting a star, release all the clamps and back up the 
C iteration counter before doing another iteration.
C
      DO 2250 I=1,NTERM
      XOLD(I)=0.0
 2250 CLAMP(I)=1.0
      CLIP=.FALSE.
      NITER=MAX0(1, NITER-1)
      GO TO 2210
C
C Now... on with the iteration.
C
 2260 IXMIN=MAX0(1, INT(XMIN-RADIUS)+1)
      IXMAX=MIN0(NCOL, INT(XMAX+RADIUS))
      IYMIN=MAX0(1, INT(YMIN-RADIUS)+1)
      IYMAX=MIN0(NROW, INT(YMAX+RADIUS))
C IXMIN, IXMAX, IYMIN, and IYMAX are now the limits of a rectangular 
C containing all pixels within one fitting radius of any star in the
C group.
C
C Zero the normal matrix and the vector of residuals.
C
      U = 1
      DO 2270 J=1,NTERM
      V(J)=0.0
      DO 2270 I=J,NTERM
 2270 C(I,J)=0.0
C
      DO 2280 I=1,NSTR
 2280 NPIX(I)=0

C
C Now deal with the pixels one by one.
C
      SUMGRP=0.
      GRPWT=0.
      DO 2390 IY=IYMIN,IYMAX
      DO 2380 IX=IXMIN,IXMAX
      IF ((DATA(IX,IY) .LT. LOBAD) .OR. (DATA(IX,IY) .GT. HIBAD))
     .     GO TO 2380
C
C If this pixel is within one fitting radius of at least one star,
C include it in the solution.  Otherwise, skip it.  While figuring
C this out, compute the squared distance of this pixel from the
C centroid of each star in the group; make sure that every star
C has at least four valid pixels within one fitting radius.
C
      OMIT=.TRUE.
      DO 2310 I=1,NSTR
      SKIP(I)=.TRUE.
      RPIXSQ(I)=(FLOAT(IX)-XC(I))**2+(FLOAT(IY)-YC(I))**2
      IF (RPIXSQ(I) .GT. CUTOFF) GO TO 2310
      SKIP(I)=.FALSE.
      OMIT=.FALSE.
 2310 CONTINUE
      IF (OMIT) GO TO 2380
C
      D=DATA(IX,IY)-SKYBAR
      WT=0.
C
C Now loop over the stars, one by one, subtracting from this pixel
C the light contribution from each star within one PSF radius.
C
      DO 2320 I=1,NSTR
      DELTAX=(XC(I)-1.)/XPSF-1.
      DELTAY=(YC(I)-1.)/YPSF-1.
C
C If this pixel is within one PSF radius of this star's center, compute 
C the scaled value of the PSF at this point and subtract it.
C
      IF (RPIXSQ(I) .GE. PSFRSQ) GO TO 2320
      VAL=USEPSF(IPSTYP, FLOAT(IX)-XC(I), FLOAT(IY)-YC(I), BRIGHT, PAR, 
     .     PSF, NPSF, NPAR, NEXP, NFRAC, DELTAX, DELTAY, DVDXC, DVDYC)
      D=D-MAG(I)*VAL
C
C The condition equation for pixel (IX,IY) is of the form
C
C data(IX,IY)-sky-summation{scale*psf(IX-Xcenter,IY-Ycenter)}=residual
C
C Then we will jigger the scale's, Xcenter's, and Ycenter's such that
C
C                Summation{weight * residual**2}
C
C is minimized.  'weight' will be a function (1) of the distance of this
C pixel from the center of the nearest star, (2) of the model-predicted
C brightness of the pixel (taking into consideration the readout noise, 
C the photons/ADU, and the interpolation error of the PSF), and (3) of 
C the size of the residual itself.  (1) is necessary to prevent the
C non-linear least-squares solution from oscillating:  oft-times it will
C come to pass that if you include a pixel in the solution, then the
C predicted shift of the centroid will cause that pixel to be excluded 
C in the next iteration, and the new predicted shift of the centroid
C will cause that pixel to be included again.  This could go on ad
C infinitum.  The cure is to have the weight of a pixel go 
C asymptotically to zero as its distance from the stellar centroid
C approaches the fitting radius.  In a case like that just described,
C the solution can then find a real minimum of the sum of the
C weighted squared residuals with that pixel at some low-weight position
C just inside the fitting radius.  (2) is just sensible weighting.
C (3) is just a crude attempt at making the solution more robust against
C bad pixels.
C
      IF (SKIP(I)) GO TO 2320
      RSQ=RPIXSQ(I)/RADSQ
      WT=AMAX1(WT, 5./(5.+RSQ/(1.-RSQ)))
      I3=I*3
      K=I3-2
      X(K)=-VAL
      K=I3-1
      X(K)=-MAG(I)*DVDXC
      X(I3)=-MAG(I)*DVDYC
 2320 CONTINUE
C
C At this point, the vector X contains the first derivative of
C the condition equation for pixel (IX,IY) with respect to each of
C the fitting parameters for all of the stars.
C
C Now these derivatives will be added into the normal matrix and the
C vector of residuals.
C
C The expected random error in the pixel is the quadratic sum of
C the Poisson statistics, plus the readout noise, plus an estimated
C error of 0.75% of the total brightness for the difficulty of flat-
C fielding and bias-correcting the chip, plus an estimated error of 
C some fraction of the fourth derivative at the peak of the profile,
C to account for the difficulty of accurately interpolating within the 
C point-spread function.  The fourth derivative of the PSF is 
C proportional to H/sigma**4 (sigma is the Gaussian width parameter for
C the stellar core); using the geometric mean of sigma(x) and sigma(y), 
C this becomes H/[sigma(x)*sigma(y)]**2.  The ratio of the fitting 
C error to this quantity is estimated from a good-seeing CTIO frame to 
C be approximately 0.027 (see definition of PKERR above.)
C
      DPOS=AMAX1(0., DATA(IX,IY)-D)
C
C DPOS = raw data minus residual = model-predicted value of the 
C intensity at this point (which presumably is non-negative).
C
      IF ((DPOS .GT. HIBAD) .AND. (NITER .GE. 4)) GO TO 2380
      SIGSQ=DPOS/PHPADU+RDNOIS+(PERERR*DPOS)**2+(PKERR*(DPOS-SKYBAR))**2
      RELERR=ABS(D)/SQRT(SIGSQ)
     
C
C Add this residual into the weighted sum of the absolute relative 
C residuals.
C
      SUMGRP=SUMGRP+RELERR*WT
      GRPWT=GRPWT+WT
C
C Add into the accumulating sums of the weighted absolute relative 
C residuals and of the image sharpness parameter for each of the stars.
C
      DO 2330 I=1,NSTR
      IF (SKIP(I)) GO TO 2330
      NPIX(I)=NPIX(I)+1
      CHI(I)=CHI(I)+RELERR*WT
      SUMWT(I)=SUMWT(I)+WT
      RHOSQ=((XC(I)-FLOAT(IX))/PAR(1))**2+
     .     ((YC(I)-FLOAT(IY))/PAR(2))**2
C
C Include in the sharpness index only those pixels within six times
C the HWHM of the centroid of the object.  (This saves time and
C floating underflows by excluding pixels which contribute less than
C about one part in a million to the index.)
C
      IF (RHOSQ .LE. 36.) THEN
        RHOSQ=0.6931472*RHOSQ
        DFDSIG=EXP(-RHOSQ)*(RHOSQ-1.)
        NUMER(I)=NUMER(I)+DFDSIG*D/SIGSQ
        DENOM(I)=DENOM(I)+DFDSIG**2/SIGSQ
      END IF
 2330 CONTINUE
C
C If clipping is in effect, reduce the weight of a bad pixel.  A pixel 
C having a residual of 2.5 sigma gets reduced to half weight; a pixel 
C having a residual of 5. sigma gets weight 1/257.
C
      WT=WT/SIGSQ
      IF (CLIP) WT=WT/(1.+(0.4*RELERR/CHIGRP)**8)
      DWT=D*WT
C     If sky is to be determined: C(NTERM,NTERM)=C(NTERM,NTERM)+WT
C     If sky is to be determined: V(NTERM)=V(NTERM)-DWT
C
C Now work this pixel into the normal matrix.
C
      DO 2370 I=1,NSTR
      IF (SKIP(I)) GO TO 2370
      I3=I*3
      I3M2=I3-2
      DO 2340 K=I3M2,I3
C     If sky is to be determined: C(NTERM,K)=C(NTERM,K)-X(K)*WT
 2340 V(K)=V(K)+X(K)*DWT
      DO 2360 J=1,I
      IF (SKIP(J)) GO TO 2360
      DO 2350 K=I3M2,I3
      DO 2350 L=3*J-2,MIN0(K, 3*J)
 2350 C(K,L)=C(K,L)+X(K)*X(L)*WT
 2360 CONTINUE
 2370 CONTINUE
C
 2380 CONTINUE
 2390 CONTINUE
C
C Make sure that every star in the group has at least four valid pixels
C within one fitting radius.
C
      REDO=.FALSE.
      DO 2400 I=1,NSTR
         IF (NPIX(I) .GE. 4) GO TO 2400
            REDO=.TRUE.
            NI=INT(ALOG10(ID(I)+0.5))+2
            IF (WATCH .GT. 0.5) THEN
               WRITE (LINE,624) NITER, NSTR, NTOT, ID(I)
               CALL OVRWRT (LINE(1:47), 3)
            ELSE
               WRITE (6,625) ID(I)
            END IF
            CALL DAORMV (I, MAXSTR, NSTR, ID, XC, YC, MAG, SKY)
            NTOT=NTOT-1
            IF (NSTR .LE. 0) GO TO 2000
            IF (WATCH .GT. 0.5) CALL TBLANK
            NTERM=NSTR*3
C           If sky is to be determined: NTERM=NTERM+1
 2400 CONTINUE
      IF (REDO) THEN
         GO TO 2210
      END IF
C
C Reflect the normal matrix across the diagonal.
C
      DO 2410 L=2,NTERM
      DO 2410 K=1,L-1
 2410 C(K,L)=C(L,K)
C
C Compute the robust estimate of the standard deviation of the
C residuals for the group as a whole, and for each star.  This 
C estimate is SQRT(PI/2) * Weighted mean absolute relative residual
C (Do you like that "absolute relative" stuff?):
C
C          CHI = 1.2533 * SUM(weight*resid)/(no. of pixels)
C
C This gets corrected for bias by being multiplied by
C
C              SQRT[(no. of pixels)/(no. of pixels - 3)].
C
      IF (GRPWT.GT.3)CHIGRP=1.2533141*SUMGRP*SQRT(1./(GRPWT*(GRPWT-3.)))
C
C But then I drive the value toward unity, depending on exactly how
C many pixels were involved:  if CHIGRP is based on exactly a total 
C weight of 3, then it is extremely poorly determined, and we just
C want to keep CHIGRP = 1.  The larger GRPWT is, the better determined
C CHIGRP is, and the less we want to force it toward unity.  So,
C just take the weighted average of CHIGRP and unity, with weights
C GRPWT-3 and 3, respectively.
C
      IF (GRPWT .GT. 3) CHIGRP=((GRPWT-3.)*CHIGRP+3.)/GRPWT
C
C CHIGRP has been pulled toward its expected value of unity to keep the 
C statistics of a small number of pixels from compeletely dominating 
C the error analysis.  Similarly, the photometric errors for the 
C individual stars will be pulled toward unity now.  Later on, if the
C number of stars in the group is greater than one, CHI will be nudged
C toward the group average.  In order to work optimally, of 
C course, this requires that PHPADU, RDNOIS, and the other noise 
C contributors which I have postulated properly represent the true 
C errors expected in each pixel.
C
C Store a smoothed CHI value for the star in SUMWT.
C
      DO 2420 I=1,NSTR
      IF (SUMWT(I) .GT. 3.) THEN
         CHI(I) = 1.2533141*CHI(I)/SQRT( (SUMWT(I)-3.)*SUMWT(I) )
         SUMWT(I) = ((SUMWT(I)-3.)*CHI(I) + 3.)/SUMWT(I)
      ELSE
         CHI(I) = CHIGRP
         SUMWT(I) = GRPWT
      END IF
 2420 CONTINUE

C
      CALL INVERS (C, MAXUNK, NTERM, IER)
      CALL VMUL (C, MAXUNK, NTERM, V, X)
      REDO=.FALSE.
      IF (NITER .LE. 1) REDO=.TRUE.
C If sky is to be determined: SKYBAR=SKYBAR-X(NTERM)
C If sky is to be determined: IF(ABS(X(NTERM)).GT.0.01)REDO=.TRUE.
C
C In the beginning, the brightness of each star will be permitted to
C change by no more than two magnitudes per iteration, and the x,y 
C coordinates of each centroid will be permitted to change by no more 
C than 0.4 pixel per iteration.  Any time that the parameter
C correction changes sign from one iteration to the next, the maximum 
C permissible change will be reduced by a factor of two.  These
C clamps are released any time a star disappears.
C
      DO 2520 I=1,NSTR
      L=3*I
      K=L-1
      J=L-2
C
C If any correction has changed sign since the last iteration, reduce
C the maximum permissible change by a factor of 2.
C
      IF (XOLD(J)*X(J) .LT. 0.) CLAMP(J)=0.5*CLAMP(J)
      IF (XOLD(K)*X(K) .LT. 0.) CLAMP(K)=0.5*CLAMP(K)
      IF (XOLD(L)*X(L) .LT. 0.) CLAMP(L)=0.5*CLAMP(L)
C
C Note that the sign of the correction is such that it must be 
C SUBTRACTED from the current value of the parameter to get the 
C improved parameter value.  This being the case, if the correction
C to the brightness is negative (the least-squares thinks that the
C star should be brighter) a change of 1 magnitude is a change of a
C factor of 2.5; if the brightness correction is positive (the star
C should be fainter) a change of 1 magnitude is a change of 60%.
C
      MAG(I)=MAG(I)-X(J)/
     .  (1.+AMAX1(X(J)/(0.84*MAG(I)),-X(J)/(5.25*MAG(I)))/CLAMP(J))
      XC(I)=XC(I)-X(K)/(1.+ABS(X(K))/(CLAMP(K)*0.4))
      YC(I)=YC(I)-X(L)/(1.+ABS(X(L))/(CLAMP(L)*0.4))
      XOLD(J)=X(J)
      XOLD(K)=X(K)
      XOLD(L)=X(L)
      MAGERR(I) = SUMWT(I)*SQRT(C(J,J))
C
C There are two milestones in the convergence process:  the fits
C proceed normally until each star's magnitude changes by less than its
C standard error or 0.005 magnitudes, whichever is greater, and its
C x- and y-centroids change by less than 0.02 pixel.  At this point
C the least-squares begins to apply the down-weighting of pixels
C with large residuals, as described above.  The fits then continue
C until each star's magnitude changes by less than MAX(0.1*standard
C error, 0.0005 magnitude), and its centroids change by less than 0.002
C pixel.
C
C If you already know that the solution hasn't converged, don't bother
C to keep checking.
C
      IF (REDO) GO TO 2510
      IF (CLIP) THEN
         IF (ABS(X(J)) .GT.
     .        AMAX1( 0.1*MAGERR(I), 0.0005*MAG(I) )) THEN
            REDO=.TRUE.
         ELSE
            DF = (0.1*SUMWT(I))**2
            IF (X(K)**2 .GT. AMAX1(DF*C(K,K), 4.E-6)) THEN
               REDO=.TRUE.
            ELSE IF (X(L)**2 .GT. AMAX1(DF*C(L,L), 4.E-6)) THEN
               REDO=.TRUE.
            END IF
         END IF
      ELSE
         IF (ABS(X(J)) .GT. 
     .        AMAX1( MAGERR(I), 0.005*MAG(I) )) THEN
            REDO=.TRUE.
         ELSE
            DF = SUMWT(I)**2
            IF (X(K)**2 .GT. AMAX1(DF*C(K,K), 4.E-4)) THEN
               REDO=.TRUE.
            ELSE IF (X(L)**2 .GT. AMAX1(DF*C(L,L), 4.E-4)) THEN
               REDO=.TRUE.
            END IF
         END IF
      END IF
 2510 CONTINUE
 2520 CONTINUE
C
C Check whether the estimated centroid of any star has moved so far out 
C of the limits of the picture that it has fewer than four or five 
C pixels within one fitting radius.
C
      I=0
 2525 I=I+1
      IF (I .GT. NSTR) GO TO 2528
      DX=AMAX1( 1.-XC(I), XC(I)-NCOL, 0.)
      DY=AMAX1( 1.-YC(I), YC(I)-NROW, 0.)
C
C If the centroid of the star is outside the picture in x or y, then
C DX or DY is its distance from the center of the edge pixel; otherwise 
C DX and DY are zero.
C
      IF ((DX .LE. 0.001) .AND. (DY .LE. 0.001)) GO TO 2525
      IF ( (DX+1.)**2+(DY+1.)**2 .LT. RADSQ) GO TO 2525
      NI=INT(ALOG10(ID(I)+0.5))+2
      IF (WATCH .GT. 0.5) THEN
         WRITE (LINE,624) NITER, NSTR, NTOT, ID(I)
         CALL OVRWRT (LINE(1:47), 3)
      ELSE
         WRITE (6,625) ID(I)
      END IF
      CALL DAORMV (I, MAXSTR, NSTR, ID, XC, YC, MAG, SKY)
      NTOT=NTOT-1
      IF (NSTR .LE. 0) GO TO 2000
      REDO=.TRUE.
C
C Update display on terminal.
C
      IF (WATCH .GT. 0.5) THEN
         WRITE (LINE,622) NITER, NSTR, NTOT
         CALL OVRWRT (LINE(1:15), 4)
      END IF
      IF (I .LT. NSTR) GO TO 2525
C
C End of loop to check that centroids aren't too far from the edge
C of the frame.
C
C Update matrix dimensions.
C
 2528 NTERM=NSTR*3
C     If sky is to be determined: NTERM=NTERM+1

C
C Now check whether any of the stars is too faint (more than 12.5
C magnitudes fainter than the PSF star).  If several stars are
C too faint, delete the faintest one, and set the brightnesses of
C the other faint ones exactly to 12.5 mag below the PSF star.
C That way on the next iteration we will see whether these stars
C want to grow or to disappear.
C
      FAINT=1.0
      IFAINT=0
C
      DO 2540 I=1,NSTR
      IF (MAG(I) .GT. 1.E-5) GO TO 2540
      IF (MAG(I) .GT. FAINT) GO TO 2530
      FAINT=MAG(I)
      IFAINT=I
 2530 MAG(I)=1.E-5
 2540 CONTINUE
C
C If at least one star is more than 12.5 mag. fainter than the
C PSF, then  I  is the index of the faintest of them, and FAINT
C is the relative brightness of the faintest of them.
C
      IF (IFAINT .GT. 0) GO TO 2560
C
C If the solution has not converged and if the number of iterations
C is less than 4, perform another iteration no questions asked.
C
      IF (REDO .AND. (NITER .LT. 4)) GO TO 2200
C
C If the solution doesn't think it has converged, after the fourth
C iteration delete the least certain star if it is less than a one-sigma
C detection; after the eighth iteration delete the least certain star if
C it is less than a 1.50 sigma detection; after the twelfth iteration
C OR if the solution thinks it has converged, delete the least certain 
C star if it is less than a two-sigma detection.
C
      FAINT=0.
      IFAINT=0
C
      DO 2550 I=1,NSTR
      WT=MAGERR(I)/MAG(I)
      IF (WT .LT. FAINT) GO TO 2550
      FAINT=WT
      IFAINT=I
 2550 CONTINUE
C
C If the solution has not converged, and the least certain star still
C has N/S less than WCRIT, do another iteration.
C
      IF (REDO .AND. (NITER .LT. 50) .AND. (FAINT .LT. WCRIT)) 
     .     GO TO 2200
C
C Either the solution has converged, or we have hit 50 iterations,
C or the poorest star has N/S greater than WCRIT.  If the star has N/S
C LESS than half, then either the solution has converged or we have
C hit 50 iterations.  In either case we're done.  Otherwise, delete
C the poorest star and do another iteration as long as there are some
C stars left.
C
      IF (FAINT .LT. 0.5) GO TO 2900
 2560 IF (WATCH .GT. 0.5) THEN
         WRITE (LINE,624) NITER, NSTR, NTOT, ID(IFAINT)
  624    FORMAT (3I6, 5X, 'Star', I8, ' has disappeared.')
         CALL OVRWRT (LINE(1:51), 3)
      ELSE
         WRITE (6,625) ID(IFAINT)
  625    FORMAT(1X, 'Star', I8, ' has disappeared.')
      END IF
      CALL DAORMV (IFAINT, MAXSTR, NSTR, ID, XC, YC, MAG, SKY)
      NTOT=NTOT-1
      IF (NSTR .LE. 0) GO TO 2000
      IF (WATCH .GT. 0.5) THEN
         WRITE (LINE,622) NITER, NSTR, NTOT
         CALL OVRWRT (LINE(1:18), 2)
      END IF
      NTERM=NSTR*3
C     If sky is to be determined: NTERM=NTERM+1
C
C After deleting a star, release all the clamps, back the iteration
C counter up by one, and do the next iteration without incrementing 
C the counter.  That way the second most uncertain star will have
C two chances to get its act together before it comes up for tenure
C review.
C
      DO 2570 I=1,NTERM
      XOLD(I)=0.0
 2570 CLAMP(I)=1.0
      CLIP=.FALSE.
      NITER=MAX0(1, NITER-1)
      GO TO 2210
C
 2900 CONTINUE
C
C Solution has either converged or gone to 50 iterations.
C
      IF ((NITER .LT. 50) .AND. (.NOT. CLIP)) THEN
C
C The first convergence milestone has been reached.  Turn on the 
C clipper, loosen the clamps, and keep iterating.
C
         CLIP=.TRUE.
         DO I=1,NTERM
            XOLD(I)=0.0
            CLAMP(I)=AMAX1(CLAMP(I), 0.25)
         END DO
         GO TO 2200
      END IF
C
C Either there have been 50 iterations or real convergence has been
C achieved.  Write out the results and go on to the next group.
C
      DO 2910 I=1,NSTR
      SHARP(I)=1.4427*PAR(1)*PAR(2)*NUMER(I)/(MAG(I)*BRIGHT*DENOM(I))
      SHARP(I)=AMIN1(99.999,AMAX1(SHARP(I),-99.999))
      ERR=1.085736*MAGERR(I)/MAG(I)
      MAG(I)=PSFMAG-1.085736*ALOG(MAG(I))
      IF (SKY(I) .LT. 9999.999) THEN
         WRITE (1,321) ID(I), XC(I), YC(I), MAG(I), ERR, SKY(I), 
     .        FLOAT(NITER), CHI(I), SHARP(I) !MCMC
  321    FORMAT (I7, 3F9.3, F9.4, F9.3, F9.0, F9.2, F9.3)
      ELSE
         WRITE (1,322) ID(I), XC(I), YC(I), MAG(I), ERR, SKY(I), 
     .     FLOAT(NITER), CHI(I), SHARP(I) !MCMC
  322    FORMAT (I7, 3F9.3, F9.4, F9.2, F9.0, F9.2, F9.3)
      END IF
 2910 CONTINUE
      WRITE (LINE, 622) NITER, NSTR, NTOT
      IF (WATCH .GT. 0.5) CALL OVRWRT (LINE(1:15), 3)
C
      WRITE (1,321)
C      GO TO 2000
 9400 CONTINUE 
C===================================================================================
C--------------------
C             MCMC Version 1.5.0 - 2023 April 4
C             S.K. Terry
C
C Markov chain Monte Carlo routine to fit stellar
C profiles (one, two, or three). Fitting parameters 
C are star centroids (x_i,y_i), flux
C ratio (f), and total flux (z).
C
C--------------------
C===================================================================================
C In order to manually adjust the size of the fitting box, uncomment and change
C the values below this line.
C      IXMIN = 1084
C      IXMAX = 1101
C      IYMIN = 1210
C      IYMAX = 1226

      GRIDSIZE = ((IXMAX-IXMIN+1)*(IYMAX-IYMIN+1))
      ALLOCATE (PPU_MIN(GRIDSIZE),FU(GRIDSIZE),FU_MIN(GRIDSIZE),
     . SGU(GRIDSIZE),SGL(GRIDSIZE),SGU2(GRIDSIZE),PPU(GRIDSIZE),
     . EMIN_CHI2(GRIDSIZE),EMINPIX(GRIDSIZE))

      A1 = "#X1_CENTER"
      A2 = "Y1_CENTER"
      A3 = "X2_CENTER"
      A4 = "Y2_CENTER"
      A5 = "SEPARATION"
      A6 = "S1_F_CONTRIB"
      A7 = "FTOTAL"
      A8 = "CHISQ"

      A9 = "X3_CENTER"
      A10 = "Y3_CENTER"
      A11 = "S1-2_SEP"
      A12 = "S1-3_SEP"
      A13 = "S1_F_CONTRIB"
      A14 = "S3_F_CONTRIB"
      A15 = "F1"
      A16 = "F2"
      A17 = "F3"

      A18 = "#X"
      A19 = "Y"
      A20 = "CHISQ"
      A21 = "INTENSITY"
      OPEN(24,FILE=MCMCALL,STATUS='UNKNOWN')
      OPEN(25,FILE=MCMCFIL,STATUS='UNKNOWN')
      OPEN(26,FILE=CHISQPIXFIL,STATUS='UNKNOWN')
      OPEN(27,FILE=BESTFIL,STATUS='UNKNOWN') 

      IX10_IN = 0.0
      IY10_IN = 0.0
      IX20_IN = 0.0
      IY20_IN = 0.0
      IF0_IN = 0.0
      ACC_ARRAY = 0.0


      WRITE(*,*)  "1, 2, or 3 star fit:"
      READ(*,*) FIT_STARS

      IF (FIT_STARS == 1) THEN
       GO TO 9500
      ELSEIF (FIT_STARS == 2) THEN
       GO TO 9600
      ELSEIF (FIT_STARS == 3) THEN
       GO TO 9700
      ELSE 
       GO TO 9900
      ENDIF
C-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=       
C==============================================================
C              1-STAR-FIT
C==============================================================
C-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=       
 9500 WRITE(*,"(a,I4.2)") "Degrees of Freedom = ",
     . ((IXMAX-IXMIN+1)*(IYMAX-IYMIN+1))-3
      WRITE(*,*) "Renormalization Factor (type '1' if first time run):"
      READ(*,*) RFAC
      WRITE(*,*)  "Number of MCMC Iterations:"
      READ(*,*) Steps
      WRITE(*,*)  "Step size pixel (goal: accept. rate -> 0.1 - 0.5):"
      READ(*,*) PSTEP
      WRITE(*,*)  "Step size flux (goal: accept. rate -> 0.1 - 0.5):"
      READ(*,*) FSTEP
      WRITE(*,*) "Star x,y:"
      READ(*,*) IX10_IN,IY10_IN
      WRITE(*,*) "Flux Contribution (0.0 - 1.0):"
      READ(*,*) IF0_IN
      WRITE(24,*) A1,A2,A7,A8
      WRITE(25,*) A1,A2,A7,A8
      WRITE(26,*) A18,A19,A20,A21
      WRITE(27,*) A1,A2,A7,A8
      WRITE(*,*) "Fit box:", IXMIN,IXMAX,IYMIN,IYMAX
      WRITE(*,*) "Sky Avg:", SKYBAR
      WRITE(*,*) "Zeropoint Mag:", PSFMAG

      IX10 = IX10_IN !initial source x
      IY10 = IY10_IN !initial source y
      IF0 = 1000*IF0_IN !initial flux ratio
      U = 1

        X1 = 1.0*IX10
        Y1 = 1.0*IY10
         F = 0.001*IF0
      DO IX=IXMIN,IXMAX
      DO IY=IYMIN,IYMAX
        FU(U)=(F )*USEPSF(IPSTYP, IX-X1, IY-Y1, BRIGHT,
     .  PAR, PSF, NPSF, NPAR, NEXP, NFRAC, DELTAX, DELTAY, DVDXC, DVDYC)
        PPU(U) = (DATA(IX,IY) - SKYBAR) !raw pixel value - sky
        U = U + 1
         ENDDO
         ENDDO
        PTOT = 0
        FTOT = 0
        U = 1
        DO IX=IXMIN,IXMAX
        DO IY=IYMIN,IYMAX
         PTOT = PTOT + PPU(U)
         FTOT = FTOT + FU(U)
         U = U + 1
         ENDDO
         ENDDO
           Z0 = (PTOT/FTOT)
           EE0 = 0.0D0
        U = 1
        DO IX=IXMIN,IXMAX
        DO IY=IYMIN,IYMAX
C---------Stetson Pixel Noise Calculations--------------------
C          D=0.795*FU(U)
C          DPOS=DATA(IX,IY)-D
C          SGU(U)=DPOS/PHPADU+RDNOIS+(PERERR*DPOS)**2+
C     .    (PKERR*(DPOS-SKYBAR))**2
C          SGL(U)=ABS(D)/SQRT(SGU(U))
C-------------------------------------------------------------
C---------Simple Pixel Noise Model----------------------------
          SGU2(U)=SQRT(16.0+MAX(DATA(IX,IY),0.)
     .    +(0.01*MAX(DATA(IX,IY),0.))**2)
C-------------------------------------------------------------
        EE0=EE0+((((DATA(IX,IY)-SKYBAR)-FU(U)*Z0)/SGU2(U))**2)/RFAC
         U = U + 1
         ENDDO
         ENDDO
        U = 1 
        DO IX=IXMIN,IXMAX
        DO IY=IYMIN,IYMAX
          FU_MIN(U) = FU(U)
          PPU_MIN(U) = (DATA(IX,IY)-SKYBAR)
            U = U + 1
            ENDDO
            ENDDO
         X1MIN = X1
         Y1MIN = Y1
         FMIN = F
         EMIN = EE0
         ZMIN = Z0
      WRITE (25,*) X1MIN, Y1MIN, ZMIN, EMIN
C--------------------------------------------------------------
C Begin MCMC LOOP
      DO IT = 1, Steps
      CALL RANDOM_NUMBER(RNARRAY)
      IX1M = IX10 + (200*RNARRAY(1)-100)*(PSTEP*0.01)
      IY1M = IY10 + (200*RNARRAY(2)-100)*(PSTEP*0.01)
22      IFM = IF0 + (200*RNARRAY(5)-100)*(FSTEP*0.1)

C      IF ((IFM .GE. 1000) .OR. (IFM .LE. 0)) THEN
C       GO TO 22
C      ELSE
C      ENDIF
C--------------------------------------------------------------
        X1 = 1.0*IX1M
        Y1 = 1.0*IY1M

          F = 0.001*IFM
          U = 1
           DO IX=IXMIN,IXMAX
           DO IY=IYMIN,IYMAX
              FU(U)=(F )*USEPSF(IPSTYP, IX-X1, IY-Y1, BRIGHT,
     .        PAR,PSF,NPSF,NPAR,NEXP,NFRAC,DELTAX,DELTAY,DVDXC,DVDYC)
              PPU(U) = (DATA(IX,IY) - SKYBAR) !raw pixel value - sky
              U = U + 1
              ENDDO
              ENDDO
           PTOT = 0.
           FTOT = 0.
           U = 1
           DO IX=IXMIN,IXMAX
           DO IY=IYMIN,IYMAX 
              PTOT = PTOT + PPU(U)
              FTOT = FTOT + FU(U)
              U = U + 1
              ENDDO
              ENDDO
           U = 1
           ZM = (PTOT/FTOT)
           EEM = 0.0D0
           DO IX=IXMIN,IXMAX
           DO IY=IYMIN,IYMAX
             EEM_CHI2=(((DATA(IX,IY)-SKYBAR)-FU(U)*ZM)/SGU2(U))**2
              EEM=EEM+(((DATA(IX,IY)-SKYBAR)-FU(U)*ZM)/SGU2(U))**2
               U = U + 1
           ENDDO
           ENDDO
         EEM = EEM/RFAC
        IF (EEM .LT. EMIN) THEN
              X1MIN = X1
              Y1MIN = Y1
              FMIN = IFM * 0.001
              EMIN = EEM
              ZMIN = ZM
         U = 1
         DO IX=IXMIN,IXMAX
         DO IY=IYMIN,IYMAX
        FU_MIN(U) = FU(U)
        PPU_MIN(U) = (DATA(IX,IY)-SKYBAR)
        EEM_CHI2=(((DATA(IX,IY)-SKYBAR)-FU(U)*ZM)/SGU2(U))**2
         U = U + 1
         WRITE(26,*) IX, IY, EEM_CHI2, DATA(IX,IY)
          ENDDO
          ENDDO
         CLOSE(26)
        ELSE
        ENDIF
      IF (EEM .LE. EE0) THEN
        IX10 = IX1M
        IY10 = IY1M
        IF0 = IFM
        Z0 = ZM
        EE0 = EEM
        ACC_ARRAY = ACC_ARRAY + 1
        WRITE(25,*) X1,Y1,ZM,EEM
         U = 1
         DO IX=IXMIN,IXMAX
         DO IY=IYMIN,IYMAX
         FU_MIN(U) = FU(U)
         PPU_MIN(U) = (DATA(IX,IY)-SKYBAR)
         U = U + 1
          ENDDO
          ENDDO
      ELSE
        PROB_RAND = RNARRAY(6)
        PROB = EXP(-(EEM-EE0)/2.00)
        IF (PROB .GT. PROB_RAND) THEN
        IX10 = IX1M
        IY10 = IY1M
        IF0 = IFM
        Z0 = ZM
        EE0 = EEM
        ACC_ARRAY = ACC_ARRAY + 1
        WRITE(25,*) X1,Y1,ZM,EEM
      ENDIF
      ENDIF
      WRITE(24,*) IX10,IY10,Z0,EE0
      IF (MOD(IT,50000)==0) WRITE(*,"(I6.2)") IT
      ENDDO !End main MCMC iteration
       WRITE(*,*) "Acceptance Rate = ", ACC_ARRAY/Steps
       IF (ACC_ARRAY/Steps .LE. 0.10) WRITE(*,*)
     . "WARNING!! Acceptance rate less than 0.10!"
       IF (ACC_ARRAY/Steps .GE. 0.50) WRITE(*,*)
     . "WARNING!! Acceptance rate larger than 0.50!"
       WRITE(*,*) "      X1               Y1             F_TOT        
     .   CHI2"
       WRITE(*,*) X1MIN,Y1MIN,
     .            ZMIN,EMIN
       WRITE(27,*) X1MIN,Y1MIN,
     .            ZMIN,EMIN
       GO TO 9000
C-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=       
C==============================================================
C              2-STAR-FIT
C==============================================================
C-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 
 9600 WRITE(*,"(a,I4.2)") "Degrees of Freedom =",
     .((IXMAX-IXMIN+1)*(IYMAX-IYMIN+1))-6
      WRITE(*,*) "Pixel Scale (mas/pix):"
      READ(*,*) PIXSCALE
      WRITE(*,*) "Renormalization Factor (type '1' if first time run):"
      READ(*,*) RFAC
      WRITE(*,*)  "Number of MCMC Iterations:"
      READ(*,*) Steps
      WRITE(*,*)  "Step size pixel (goal: accept. rate -> 0.1 - 0.5):"
      READ(*,*) PSTEP
      WRITE(*,*)  "Step size flux (goal: accept. rate -> 0.1 - 0.5):"
      READ(*,*) FSTEP
      WRITE(*,*) "Star 1 x,y:"
      READ(*,*) IX10_IN,IY10_IN
      WRITE(*,*) "Star 2 x,y:"
      READ(*,*) IX20_IN,IY20_IN
      WRITE(*,*) "Star 1 Flux Contribution (star 1/(star1 + star2)):"
      READ(*,*) IF0_IN
      WRITE(24,*) A1,A2,A3,A4,A5,A6,A7,A15,A16,A8
      WRITE(25,*) A1,A2,A3,A4,A5,A6,A7,A15,A16,A8
      WRITE(26,*) A18,A19,A20,A21
      WRITE(27,*) A1,A2,A3,A4,A5,A6,A7,A15,A16,A8
      WRITE(*,*) "Fit box:", IXMIN,IXMAX,IYMIN,IYMAX
      WRITE(*,*) "Sky Avg:", SKYBAR
      WRITE(*,*) "Zeropoint Mag:", PSFMAG

      IX10 = IX10_IN !initial source x
      IY10 = IY10_IN !initial source y
      IX20 = IX20_IN !initial lens x
      IY20 = IY20_IN !initial lens y
      IF0 = 1000*IF0_IN !initial flux ratio
      U = 1

        X1 = 1.0*IX10
        Y1 = 1.0*IY10
        X2 = 1.0*IX20
        Y2 = 1.0*IY20
         F = 0.001*IF0
      DO IX=IXMIN,IXMAX
      DO IY=IYMIN,IYMAX
        FU(U)=(F )*USEPSF(IPSTYP, IX-X1, IY-Y1, BRIGHT,
     .  PAR, PSF, NPSF, NPAR, NEXP, NFRAC, DELTAX, DELTAY, DVDXC, DVDYC)
     .  + (1-F)*USEPSF(IPSTYP, IX-X2, IY-Y2, BRIGHT,
     .  PAR, PSF, NPSF, NPAR, NEXP, NFRAC, DELTAX, DELTAY, DVDXC, DVDYC)
        PPU(U) = (DATA(IX,IY) - SKYBAR) !raw pixel value - sky
        U = U + 1
         ENDDO
         ENDDO
        PTOT = 0
        FTOT = 0
        U = 1
        DO IX=IXMIN,IXMAX
        DO IY=IYMIN,IYMAX
         PTOT = PTOT + PPU(U)
         FTOT = FTOT + FU(U)
         U = U + 1
         ENDDO
         ENDDO
           Z0 = (PTOT/FTOT)
           EE0 = 0.0D0
        U = 1
        DO IX=IXMIN,IXMAX
        DO IY=IYMIN,IYMAX
C---------Stetson Pixel Noise Calculations--------------------
C          D=0.795*FU(U)
C          DPOS=DATA(IX,IY)-D
C          SGU(U)=DPOS/PHPADU+RDNOIS+(PERERR*DPOS)**2+
C     .    (PKERR*(DPOS-SKYBAR))**2
C          SGL(U)=ABS(D)/SQRT(SGU(U))
C-------------------------------------------------------------
C---------Simple Pixel Noise Model----------------------------
          SGU2(U)=SQRT(16.0+MAX(DATA(IX,IY),0.)
     .    +(0.01*MAX(DATA(IX,IY),0.))**2)
C-------------------------------------------------------------
        EE0=EE0+((((DATA(IX,IY)-SKYBAR)-FU(U)*Z0)/SGU2(U))**2)/RFAC
         U = U + 1
         ENDDO
         ENDDO
C         EE0 = EE0 + EXP(((FU(U)*Z0-6284)/231.565)**2)!Flux Constraint here
C      SSEP = PIXSCALE*SQRT((X1-X2)**2+(Y1-Y2)**2)
C        EE0 = EE0 + EXP(((SSEP-75.0)/(5.73))**2)!SSEP Constraint Here
        U = 1
        DO IX=IXMIN,IXMAX
        DO IY=IYMIN,IYMAX
          FU_MIN(U) = FU(U)
          PPU_MIN(U) = (DATA(IX,IY)-SKYBAR)
            U = U + 1
            ENDDO
            ENDDO
         X1MIN = X1
         Y1MIN = Y1
         X2MIN = X2
         Y2MIN = Y2
         FMIN = F
         EMIN = EE0
         ZMIN = Z0
      WRITE (25,*) X1MIN, Y1MIN, X2MIN, Y2MIN,
     . 9.942*SQRT((X1MIN-X2MIN)**2+(Y1MIN-Y2MIN)**2), FMIN, ZMIN, 
     . (FMIN*ZMIN), ((1-FMIN)*(ZMIN)), EMIN
C--------------------------------------------------------------
C Begin MCMC LOOP
      DO IT = 1, Steps
      CALL RANDOM_NUMBER(RNARRAY)
      IX1M = IX10 + (200*RNARRAY(1)-100)*(PSTEP*0.01)
      IY1M = IY10 + (200*RNARRAY(2)-100)*(PSTEP*0.01)
      IX2M = IX20 + (200*RNARRAY(3)-100)*(PSTEP*0.01)
13      IY2M = IY20 + (200*RNARRAY(4)-100)*(PSTEP*0.01)
12      IFM = IF0 + (200*RNARRAY(5)-100)*(FSTEP*0.1)

C      IF ((IFM .GE. 1000) .OR. (IFM .LE. 0)) THEN
C       GO TO 12
C      ELSE
C      ENDIF
C--------------------------------------------------------------
        X1 = 1.0*IX1M
        Y1 = 1.0*IY1M
        X2 = 1.0*IX2M
        Y2 = 1.0*IY2M

          F = 0.001*IFM
          U = 1
           DO IX=IXMIN,IXMAX
           DO IY=IYMIN,IYMAX
              FU(U)=(F )*USEPSF(IPSTYP, IX-X1, IY-Y1, BRIGHT,
     .        PAR,PSF,NPSF,NPAR,NEXP,NFRAC,DELTAX,DELTAY,DVDXC,DVDYC)
     .        + (1-F)*USEPSF(IPSTYP, IX-X2, IY-Y2, BRIGHT,
     .        PAR,PSF,NPSF,NPAR,NEXP,NFRAC,DELTAX,DELTAY,DVDXC,DVDYC)
              PPU(U) = (DATA(IX,IY) - SKYBAR) !raw pixel value - sky
              U = U + 1
              ENDDO
              ENDDO
           PTOT = 0.
           FTOT = 0.
           U = 1
           DO IX=IXMIN,IXMAX
           DO IY=IYMIN,IYMAX 
              PTOT = PTOT + PPU(U)
              FTOT = FTOT + FU(U)
              U = U + 1
              ENDDO
              ENDDO
           U = 1
           ZM = (PTOT/FTOT)
           EEM = 0.0D0
           DO IX=IXMIN,IXMAX
           DO IY=IYMIN,IYMAX
             EEM_CHI2=(((DATA(IX,IY)-SKYBAR)-FU(U)*ZM)/SGU2(U))**2
              EEM=EEM+(((DATA(IX,IY)-SKYBAR)-FU(U)*ZM)/SGU2(U))**2
               U = U + 1
           ENDDO
           ENDDO
C         EEM = EEM + EXP(((FU(U)*Z0-6284)/231.565)**2)!Flux Constraint here
C      SSEP = PIXSCALE*SQRT((X1-X2)**2+(Y1-Y2)**2)
C        EEM = EEM + EXP(((SSEP-75.0)/(5.73))**2)!SSEP Constraint Here
         EEM = EEM/RFAC
        IF (EEM .LT. EMIN) THEN
              X1MIN = X1
              Y1MIN = Y1
              X2MIN = X2
              Y2MIN = Y2
              FMIN = IFM * 0.001
              ZMIN = ZM
              EMIN = EEM
              SSEPMIN = PIXSCALE*SQRT((X1MIN-X2MIN)**2+(Y1MIN-Y2MIN)**2)
         U = 1
         DO IX=IXMIN,IXMAX
         DO IY=IYMIN,IYMAX
        FU_MIN(U) = FU(U)
        PPU_MIN(U) = (DATA(IX,IY)-SKYBAR)
        EEM_CHI2=(((DATA(IX,IY)-SKYBAR)-FU(U)*ZM)/SGU2(U))**2
         U = U + 1
         WRITE(26,*) IX, IY, EEM_CHI2, DATA(IX,IY)
          ENDDO
          ENDDO
         CLOSE(26)
        ELSE
        ENDIF
      IF (EEM .LE. EE0) THEN
        IX10 = IX1M
        IX20 = IX2M
        IY10 = IY1M
        IY20 = IY2M
        IF0 = IFM
        Z0 = ZM
        EE0 = EEM
        SSEP = PIXSCALE*SQRT((X1-X2)**2+(Y1-Y2)**2)
        ACC_ARRAY = ACC_ARRAY + 1
        WRITE(25,*) X1,Y1,X2,Y2,SSEP,IFM*0.001,ZM,
     .  ((IFM*.001)*(ZM)), ((1-IFM*.001)*(ZM)),EEM
         U = 1
         DO IX=IXMIN,IXMAX
         DO IY=IYMIN,IYMAX
         FU_MIN(U) = FU(U)
         PPU_MIN(U) = (DATA(IX,IY)-SKYBAR)
         U = U + 1
          ENDDO
          ENDDO
      ELSE
        PROB_RAND = RNARRAY(6)
        PROB = EXP(-(EEM-EE0)/2.00)
        IF (PROB .GT. PROB_RAND) THEN
        IX10 = IX1M
        IX20 = IX2M
        IY10 = IY1M
        IY20 = IY2M
        IF0 = IFM
        Z0 = ZM
        EE0 = EEM
        SSEP = PIXSCALE*SQRT((X1-X2)**2+(Y1-Y2)**2)
        ACC_ARRAY = ACC_ARRAY + 1
        WRITE(25,*) X1,Y1,X2,Y2,SSEP,IFM*0.001,ZM,
     .  ((IFM*.001)*(ZM)), ((1-IFM*.001)*(ZM)),EEM
        ENDIF
      ENDIF
      WRITE(24,*) IX10,IY10,IX20,IY20,SSEP,IF0*0.001,Z0,
     .  ((IF0*.001)*(Z0)), ((1-IF0*.001)*(Z0)),EE0
      IF (MOD(IT,50000)==0) WRITE(*,"(I8.2)") IT
      ENDDO !End main MCMC iteration
       WRITE(*,*) "Acceptance Rate = ", ACC_ARRAY/Steps
       IF (ACC_ARRAY/Steps .LE. 0.10) WRITE(*,*)
     . "WARNING!! Acceptance rate less than 0.10!"
       IF (ACC_ARRAY/Steps .GE. 0.50) WRITE(*,*)
     . "WARNING!! Acceptance rate larger than 0.50!"
       WRITE(*,*) "      X1               Y1               X2        
     .    Y2            SEP          S1_F_CONTRIB        F_TOTAL
     .           F1               F2               CHI2"
       WRITE(*,*) X1MIN,Y1MIN,X2MIN,Y2MIN,
     .            PIXSCALE*SQRT((X1MIN-X2MIN)**2+(Y1MIN-Y2MIN)**2),
     .            FMIN,ZMIN,(FMIN*ZMIN),((1-FMIN)*(ZMIN)),EMIN
       WRITE(27,*) X1MIN,Y1MIN,X2MIN,Y2MIN,
     .            PIXSCALE*SQRT((X1MIN-X2MIN)**2+(Y1MIN-Y2MIN)**2),
     .            FMIN,ZMIN,(FMIN*ZMIN),((1-FMIN)*(ZMIN)),EMIN
       GO TO 9000
C-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=       
C==============================================================
C              3-STAR-FIT
C==============================================================
C-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=       
 9700 WRITE(*,"(a,I4.2)") "Degrees of Freedom = ",
     . ((IXMAX-IXMIN+1)*(IYMAX-IYMIN+1))-9
      WRITE(*,*) "Pixel Scale (mas/pix):"
      READ(*,*) PIXSCALE
      WRITE(*,*) "Renormalization Factor (type '1' if first time run):"
      READ(*,*) RFAC
      WRITE(*,*)  "Number of MCMC Iterations:"
      READ(*,*) Steps
      WRITE(*,*)  "Step size pixel (goal: accept. rate -> 0.2 - 0.5):"
      READ(*,*) PSTEP
      WRITE(*,*)  "Step size flux (goal: accept. rate -> 0.2 - 0.5):"
      READ(*,*) FSTEP
      WRITE(*,*) "Star 1 x,y:"
      READ(*,*) IX10_IN,IY10_IN
      WRITE(*,*) "Star 2 x,y:"
      READ(*,*) IX20_IN,IY20_IN
      WRITE(*,*) "Star 3 x,y:"
      READ(*,*) IX30_IN,IY30_IN
      WRITE(*,*) "Star1 Flux Contrib. (star 1/(star1 + star2 + star3)):"
      READ(*,*) IF0_IN
      WRITE(*,*) "Star3 Flux Contrib. (star 3/(star1 + star2 + star3)):"
      WRITE(*,*) "This should be the FAINTEST star!"
      READ(*,*) IFB0_IN
      WRITE(24,*) A1,A2,A3,A4,A9,A10,A11,A12,A13,A14,A7,A15,A16,A17,A8
      WRITE(25,*) A1,A2,A3,A4,A9,A10,A11,A12,A13,A14,A7,A15,A16,A17,A8
      WRITE(26,*) A18,A19,A20,A21
      WRITE(27,*) A1,A2,A3,A4,A9,A10,A11,A12,A13,A14,A7,A15,A16,A17,A8
      WRITE(*,*) "Fit box:", IXMIN,IXMAX,IYMIN,IYMAX
      WRITE(*,*) "Sky Avg:", SKYBAR
      WRITE(*,*) "Zeropoint Mag:", PSFMAG

      IX10 = IX10_IN !initial source x
      IY10 = IY10_IN !initial source y
      IX20 = IX20_IN !initial lens x
      IY20 = IY20_IN !initial lens y
      IX30 = IX30_IN !initial blend x
      IY30 = IY30_IN !initial blend y
      IF0 = 1000*IF0_IN !initial flux ratio
      IFB0 = 1000*IFB0_IN !initial flux ratio
      U = 1

        X1 = 1.0*IX10
        Y1 = 1.0*IY10
        X2 = 1.0*IX20
        Y2 = 1.0*IY20
        X3 = 1.0*IX30
        Y3 = 1.0*IY30
         F = 0.001*IF0
         FB = 0.001*IFB0
      DO IX=IXMIN,IXMAX
      DO IY=IYMIN,IYMAX
        FU(U)=(F )*USEPSF(IPSTYP, IX-X1, IY-Y1, BRIGHT,
     .  PAR, PSF, NPSF, NPAR, NEXP, NFRAC, DELTAX, DELTAY, DVDXC, DVDYC)
     .  + (1-F-FB)*USEPSF(IPSTYP, IX-X2, IY-Y2, BRIGHT,
     .  PAR, PSF, NPSF, NPAR, NEXP, NFRAC, DELTAX, DELTAY, DVDXC, DVDYC)
     .  + (FB )*USEPSF(IPSTYP, IX-X3, IY-Y3, BRIGHT,
     .  PAR, PSF, NPSF, NPAR, NEXP, NFRAC, DELTAX, DELTAY, DVDXC, DVDYC)
        PPU(U) = (DATA(IX,IY) - SKYBAR) !raw pixel value - sky
        U = U + 1
         ENDDO
         ENDDO
        PTOT = 0
        FTOT = 0
        U = 1
        DO IX=IXMIN,IXMAX
        DO IY=IYMIN,IYMAX
         PTOT = PTOT + PPU(U)
         FTOT = FTOT + FU(U)
         U = U + 1
         ENDDO
         ENDDO
           Z0 = (PTOT/FTOT)
           EE0 = 0.0D0
        U = 1
        DO IX=IXMIN,IXMAX
        DO IY=IYMIN,IYMAX
C---------Stetson Pixel Noise Calculations--------------------
C          D=0.795*FU(U)
C          DPOS=DATA(IX,IY)-D
C          SGU(U)=DPOS/PHPADU+RDNOIS+(PERERR*DPOS)**2+
C     .    (PKERR*(DPOS-SKYBAR))**2
C          SGL(U)=ABS(D)/SQRT(SGU(U))
C-------------------------------------------------------------
C---------Simple Pixel Noise Model----------------------------
          SGU2(U)=SQRT(16.0+MAX(DATA(IX,IY),0.)
     .    +(0.01*MAX(DATA(IX,IY),0.))**2)
C-------------------------------------------------------------
        EE0=EE0+((((DATA(IX,IY)-SKYBAR)-FU(U)*Z0)/SGU2(U))**2)/RFAC
         U = U + 1
         ENDDO
         ENDDO
C         EEM = EEM + EXP(((FU(U)*Z0-6284)/231.565)**2)!Flux Constraint here
      SSEP = PIXSCALE*SQRT((X2-X3)**2+(Y2-Y3)**2)
      SSEP2 = PIXSCALE*SQRT((X2-X3)**2+(Y2-Y3)**2)
        EE0 = EE0 + EXP(((SSEP-26.50)/(2.50))**2)
     .  + EXP(((SSEP2-44.00)/(10.00))**2)!SSEP Constraint Here
        U = 1 
        DO IX=IXMIN,IXMAX
        DO IY=IYMIN,IYMAX
          FU_MIN(U) = FU(U)
          PPU_MIN(U) = (DATA(IX,IY)-SKYBAR)
            U = U + 1
            ENDDO
            ENDDO
         X1MIN = X1
         Y1MIN = Y1
         X2MIN = X2
         Y2MIN = Y2
         X3MIN = X3
         Y3MIN = Y3
         FMIN = F
         FBMIN = FB
         EMIN = EE0
         ZMIN = Z0
      WRITE (25,*) X1MIN, Y1MIN, X2MIN, Y2MIN, X3MIN, Y3MIN,
     . PIXSCALE*SQRT((X1MIN-X2MIN)**2+(Y1MIN-Y2MIN)**2),
     . PIXSCALE*SQRT((X1MIN-X3MIN)**2+(Y1MIN-Y3MIN)**2), FMIN, FBMIN,
     . ZMIN, (FMIN*ZMIN), ((1-FMIN)*(ZMIN)), (FBMIN*ZMIN), EMIN
C--------------------------------------------------------------
C Begin MCMC LOOP
      DO IT = 1, Steps
      CALL RANDOM_NUMBER(RNARRAY)
      IX1M = IX10 + (200*RNARRAY(1)-100)*(PSTEP*0.01)
      IY1M = IY10 + (200*RNARRAY(2)-100)*(PSTEP*0.01)
      IX2M = IX20 + (200*RNARRAY(3)-100)*(PSTEP*0.01)
      IY2M = IY20 + (200*RNARRAY(4)-100)*(PSTEP*0.01)
      IX3M = IX30 + (200*RNARRAY(5)-100)*(PSTEP*0.01)
03      IY3M = IY30 + (200*RNARRAY(6)-100)*(PSTEP*0.01)
02      IFM = IF0 + (200*RNARRAY(7)-100)*(FSTEP*0.1)
01      IFBM = IFB0 + (200*RNARRAY(8)-100)*(FSTEP*0.1)

C      IF ((IFM .GE. 1000) .OR. (IFM .LE. 0)) THEN
C       GO TO 02
C      ELSE
C      ENDIF
C--------------------------------------------------------------
        X1 = 1.0*IX1M
        Y1 = 1.0*IY1M
        X2 = 1.0*IX2M
        Y2 = 1.0*IY2M
        X3 = 1.0*IX3M
        Y3 = 1.0*IY3M

          F = 0.001*IFM
          FB = 0.001*IFBM
          U = 1
           DO IX=IXMIN,IXMAX
           DO IY=IYMIN,IYMAX
              FU(U)=(F )*USEPSF(IPSTYP, IX-X1, IY-Y1, BRIGHT,
     .        PAR,PSF,NPSF,NPAR,NEXP,NFRAC,DELTAX,DELTAY,DVDXC,DVDYC)
     .        + (1-(F+FB))*USEPSF(IPSTYP, IX-X2, IY-Y2, BRIGHT,
     .        PAR,PSF,NPSF,NPAR,NEXP,NFRAC,DELTAX,DELTAY,DVDXC,DVDYC)
     .        + (FB )*USEPSF(IPSTYP, IX-X3, IY-Y3, BRIGHT,
     .        PAR,PSF,NPSF,NPAR,NEXP,NFRAC,DELTAX,DELTAY,DVDXC,DVDYC)
              PPU(U) = (DATA(IX,IY) - SKYBAR) !raw pixel value - sky
              U = U + 1
              ENDDO
              ENDDO
           PTOT = 0.
           FTOT = 0.
           U = 1
           DO IX=IXMIN,IXMAX
           DO IY=IYMIN,IYMAX 
              PTOT = PTOT + PPU(U)
              FTOT = FTOT + FU(U)
              U = U + 1
              ENDDO
              ENDDO
           U = 1
           ZM = (PTOT/FTOT)
           EEM = 0.0D0
           DO IX=IXMIN,IXMAX
           DO IY=IYMIN,IYMAX
             EEM_CHI2=(((DATA(IX,IY)-SKYBAR)-FU(U)*ZM)/SGU2(U))**2
              EEM=EEM+(((DATA(IX,IY)-SKYBAR)-FU(U)*ZM)/SGU2(U))**2
               U = U + 1
           ENDDO
           ENDDO
C         EEM = EEM + EXP(((FU(U)*Z0-6284)/231.565)**2)!Flux Constraint here
      SSEP = PIXSCALE*SQRT((X2-X3)**2+(Y2-Y3)**2)
      SSEP2 = PIXSCALE*SQRT((X1-X3)**2+(Y1-Y3)**2)
        EEM = EEM + EXP(((SSEP-26.50)/(2.50))**2)
     .  + EXP(((SSEP2-44.00)/(10.00))**2)!SSEP Constraint Here
         EEM = EEM/RFAC
        IF (EEM .LT. EMIN) THEN
              X1MIN = X1
              Y1MIN = Y1
              X2MIN = X2
              Y2MIN = Y2
              X3MIN = X3
              Y3MIN = Y3
              FMIN = IFM * 0.001
              FBMIN = IFBM * 0.001
              EMIN = EEM
              ZMIN = ZM
              LSSEPMIN=PIXSCALE*SQRT((X1MIN-X2MIN)**2+(Y1MIN-Y2MIN)**2)
              BSEPMIN=PIXSCALE*SQRT((X1MIN-X3MIN)**2+(Y1MIN-Y3MIN)**2)
         U = 1                 
         DO IX=IXMIN,IXMAX
         DO IY=IYMIN,IYMAX
        FU_MIN(U) = FU(U)
        PPU_MIN(U) = (DATA(IX,IY)-SKYBAR)
        EEM_CHI2=(((DATA(IX,IY)-SKYBAR)-FU(U)*ZM)/SGU2(U))**2
         U = U + 1
         WRITE(26,*) IX, IY, EEM_CHI2, DATA(IX,IY)
          ENDDO
          ENDDO
         CLOSE(26)
        ELSE
        ENDIF
      IF (EEM .LE. EE0) THEN
        IX10 = IX1M
        IX20 = IX2M
        IX30 = IX3M
        IY10 = IY1M
        IY20 = IY2M
        IY30 = IY3M
        IF0 = IFM
        IFB0 = IFBM
        Z0 = ZM
        EE0 = EEM
        LSSEP = PIXSCALE*SQRT((X1-X2)**2+(Y1-Y2)**2)
        BSEP = PIXSCALE*SQRT((X1-X3)**2+(Y1-Y3)**2)
        ACC_ARRAY = ACC_ARRAY + 1
        WRITE(25,*) X1,Y1,X2,Y2,X3,Y3,LSSEP,BSEP,IFM*0.001,IFBM*0.001,
     .  ZM,((IFM*.001)*(ZM)),((1-IFM*.001)*(ZM)),((IFBM*.001)*(ZM)),EEM
         U = 1
         DO IX=IXMIN,IXMAX
         DO IY=IYMIN,IYMAX
         FU_MIN(U) = FU(U)
         PPU_MIN(U) = (DATA(IX,IY)-SKYBAR)
         U = U + 1
          ENDDO
          ENDDO
      ELSE
        PROB_RAND = RNARRAY(9)
        PROB = EXP(-(EEM-EE0)/2.00)
        IF (PROB .GT. PROB_RAND) THEN
        IX10 = IX1M
        IX20 = IX2M
        IX30 = IX3M
        IY10 = IY1M
        IY20 = IY2M
        IY30 = IY3M
        IF0 = IFM
        IFB0 = IFBM
        Z0 = ZM
        EE0 = EEM
        LSSEP = PIXSCALE*SQRT((X1-X2)**2+(Y1-Y2)**2)
        BSEP = PIXSCALE*SQRT((X1-X3)**2+(Y1-Y3)**2)
        ACC_ARRAY = ACC_ARRAY + 1
        WRITE(25,*) X1,Y1,X2,Y2,X3,Y3,LSSEP,BSEP,IFM*0.001,IFBM*0.001,
     .  ZM,((IFM*.001)*(ZM)),((1-IFM*.001)*(ZM)),((IFBM*.001)*(ZM)),EEM
      ENDIF
      ENDIF
      WRITE(24,*) IX10,IY10,IX20,IY20,IX30,IY30,
     . LSSEP,BSEP,IF0*0.001,IFB0*0.001,Z0,
     . ((IF0*.001)*(Z0)),((1-IF0*.001)*(Z0)),((IFB0*.001)*(Z0)),EE0
      IF (MOD(IT,50000)==0) WRITE(*,"(I6.2)") IT
      ENDDO !End main MCMC iteration
       WRITE(*,*) "Acceptance Rate = ", ACC_ARRAY/Steps
       IF (ACC_ARRAY/Steps .LE. 0.10) WRITE(*,*)
     . "WARNING!! Acceptance rate less than 0.10!"
       IF (ACC_ARRAY/Steps .GE. 0.50) WRITE(*,*)
     . "WARNING!! Acceptance rate larger than 0.50!"
       WRITE(*,*) "      X1               Y1               X2        
     .    Y2               X3               Y3           1-2SEP         
     .   1-3SEP      S1_F_CONTRIB     S3_F_CONTRIB      F_TOTAL
     .     F1              F2           F3              CHI2"
       WRITE(*,*) X1MIN,Y1MIN,X2MIN,Y2MIN,X3MIN,Y3MIN,
     . PIXSCALE*SQRT((X1MIN-X2MIN)**2+(Y1MIN-Y2MIN)**2),
     . PIXSCALE*SQRT((X1MIN-X3MIN)**2+(Y1MIN-Y3MIN)**2),
     .       FMIN,FBMIN,ZMIN,(FMIN*ZMIN),((1-FMIN)*(ZMIN)),
     .       (FBMIN*ZMIN),EMIN
       WRITE(27,*) X1MIN,Y1MIN,X2MIN,Y2MIN,X3MIN,Y3MIN,
     . PIXSCALE*SQRT((X1MIN-X2MIN)**2+(Y1MIN-Y2MIN)**2),
     . PIXSCALE*SQRT((X1MIN-X3MIN)**2+(Y1MIN-Y3MIN)**2),
     .       FMIN,FBMIN,ZMIN,(FMIN*ZMIN),((1-FMIN)*(ZMIN)),
     .       (FBMIN*ZMIN),EMIN
       GO TO 9000
C------------------------------------
C----------END OF MCMC-------------------------------------------------
C------------------------------------
C 2130 CONTINUE
C
C
C
C-----------------------------------------------------------------------
C
C Normal return.
C
 9000 CONTINUE
      CLOSE(24)
      CLOSE(25)
      CLOSE(26)
      CLOSE(27)
      CALL CLFILE (1)
C      CALL CLFILE (11)
      CALL STUPID ('   Done.  ')
      CALL CLFILE (2)
      CALL CLFILE (3)
C      CALL CLFILE (4)
C      CALL CLFILE (5)
      RETURN
C
C-----------------------------------------------------------------------
C
C Irrecoverable errors.
C
 9100 CONTINUE
      CALL STUPID ('Error reading picture.')
      RETURN
C
 9200 CONTINUE
      CALL STUPID ('Not a group file.')
      RETURN
C
 9300 CONTINUE
      CALL STUPID ('Error opening file.')
      RETURN
C
 9900 CONTINUE
      CALL STUPID ('Please type 1, 2, or 3.')
      RETURN
C
      END

C#######################################################################
C
      SUBROUTINE  DAORMV (I, MAXSTR, NSTR, ID, XC, YC, MAG, SKY)
C
C=======================================================================
C
C A simple little subroutine to remove the I-th star from a group.
C The other arguments are obvious.
C
C                        1991 April 18
C
C=======================================================================
C
      IMPLICIT NONE
      INTEGER MAXSTR
      REAL XC(MAXSTR), YC(MAXSTR), MAG(MAXSTR), SKY(MAXSTR)
      INTEGER ID(MAXSTR)
C
      INTEGER I, NSTR
C
C-----------------------------------------------------------------------
C
C If we are trying to delete the last star in the group, all we need to
C do is reduce NSTR by one.  Otherwise, overwrite the I-th star with the
C NSTR-th star, and THEN reduce NSTR by one.
C
      IF (I .EQ. NSTR) GO TO 1000
      ID(I)=ID(NSTR)
      MAG(I)=MAG(NSTR)
      XC(I)=XC(NSTR)
      YC(I)=YC(NSTR)
      SKY(I)=SKY(NSTR)
 1000 NSTR=NSTR-1
      RETURN
      END
