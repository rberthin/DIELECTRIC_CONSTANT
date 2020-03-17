PROGRAM CONV_DIPOLES

    IMPLICIT NONE
    INTEGER :: io, nb_lines, nb_mol, nb_atoms, nb_steps
    INTEGER :: i, m, a, w
    DOUBLE PRECISION :: boxx, boxy, boxz, RESULT, AVERAGE
    DOUBLE PRECISION :: conv, bohr, epsilon0, kbT_298
    DOUBLE PRECISION :: TOT_XX, TOT_YY, TOT_ZZ, VOLUME
    DOUBLE PRECISION :: xdiff, ydiff, zdiff, DIFF_X, DIFF_Y, DIFF_Z
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: SOMME_MUX, SOMME_MUY, SOMME_MUZ, charge
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: step, PX, PY, PZ, TOT_P, P_SQUARED
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: xref, yref, zref, PERM_X, PERM_Y, PERM_Z
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: TOT_X, TOT_Y, TOT_Z
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: MUX, MUY, MUZ, X, Y, Z
    CHARACTER(LEN=25), ALLOCATABLE, DIMENSION (:,:) :: nom
    CHARACTER*200 :: crap
    CHARACTER(LEN=30) :: input_file_dipoles, input_file_data, input_file_pol

    !****** CONSTANTS ******
    conv = 57.21476575  !conv ua en C/m2
    bohr = 0.529177249
    epsilon0 = 8.85418782D-12
    kbT_298 = 4.1143325895999994D-21
    !***********************

    OPEN(unit = 10, file = "epsilon.inpt", status='old', iostat=io)
    
    ! READ INPUT FILE
    !--------------------------
    READ(10,*) input_file_dipoles
    READ(10,*) input_file_data
    READ(10,*) input_file_pol
    READ(10,*) nb_mol
    READ(10,*) nb_atoms
    ALLOCATE(charge(nb_atoms))
    DO a = 1, nb_atoms
       READ(10,*) charge(a)
    ENDDO
    READ(10,*) nb_steps
    !--------------------------
    OPEN(unit = 11, file = input_file_dipoles, status='old', iostat=io)
    OPEN(unit = 12, file = input_file_data, status='old', iostat=io)
    OPEN(unit = 13, file = input_file_pol, status='old', iostat=io)
    OPEN(unit = 20, file = "polarization_total.dat", iostat=io)
    ALLOCATE(MUX(nb_mol, nb_atoms))
    ALLOCATE(MUY(nb_mol, nb_atoms))
    ALLOCATE(MUZ(nb_mol, nb_atoms))
    ALLOCATE(X(nb_mol, nb_atoms))
    ALLOCATE(Y(nb_mol, nb_atoms))
    ALLOCATE(Z(nb_mol, nb_atoms))
    ALLOCATE(SOMME_MUX(nb_mol))
    ALLOCATE(SOMME_MUY(nb_mol))
    ALLOCATE(SOMME_MUZ(nb_mol))
    ALLOCATE(xref(nb_mol))
    ALLOCATE(yref(nb_mol))
    ALLOCATE(zref(nb_mol))
    ALLOCATE(nom(nb_mol, nb_atoms))
    ALLOCATE(PERM_X(nb_mol))
    ALLOCATE(PERM_Y(nb_mol))
    ALLOCATE(PERM_Z(nb_mol))
    ALLOCATE(TOT_X(nb_mol))
    ALLOCATE(TOT_Y(nb_mol))
    ALLOCATE(TOT_Z(nb_mol))
    
    SOMME_MUX = 0.0D0
    SOMME_MUY = 0.0D0
    SOMME_MUZ = 0.0D0
    PERM_X = 0.0D0
    PERM_Y = 0.0D0
    PERM_Z = 0.0D0

    ! READ DIPOLES.OUT
    !----------------------------------------------
    DO i = 1, 5
       READ(11,*) crap
    ENDDO
    DO a = 1, nb_atoms     
      DO m = 1, nb_mol
        READ(11,*) MUX(m,a), MUY(m,a), MUZ(m,a) 
      ENDDO
    ENDDO
    CLOSE(11)
    !----------------------------------------------
    
    ! READ DATA.INPT
    !----------------------------------------------
    DO i = 1, 6
       READ(12,*)
    ENDDO
    READ(12,*) boxx, boxy, boxz
    READ(12,*)
    
    DO a = 1, nb_atoms
      DO m = 1, nb_mol
        READ(12,*) nom(m,a), X(m,a), Y(m,a), Z(m,a)
        IF (a == 2) THEN
          xref(m) = X(m,a)
          yref(m) = Y(m,a)
          zref(m) = Z(m,a)
        ENDIF
      ENDDO
    ENDDO
    CLOSE(12)
    DO m = 1, nb_mol
       DO a = 1, nb_atoms
          xdiff = X(m,a) - xref(m)
          ydiff = Y(m,a) - yref(m)
          zdiff = Z(m,a) - zref(m)
          IF (xdiff > boxx/2) THEN
             X(m,a) = X(m,a) - boxx
          ELSE IF (xdiff < -boxx/2) THEN
             X(m,a) = X(m,a) + boxx
          ENDIF
          IF (ydiff > boxy/2) THEN
             Y(m,a) = Y(m,a) - boxy
          ELSE IF (ydiff < -boxy/2) THEN
             Y(m,a) = Y(m,a) + boxy
          ENDIF
          IF (zdiff > boxz/2) THEN
             Z(m,a) = Z(m,a) - boxz
          ELSE IF (zdiff < -boxz/2) THEN
             Z(m,a) = Z(m,a) + boxz
          ENDIF
          SOMME_MUX(m) = SOMME_MUX(m) + MUX(m,a)
          SOMME_MUY(m) = SOMME_MUY(m) + MUY(m,a)
          SOMME_MUZ(m) = SOMME_MUZ(m) + MUZ(m,a)
          PERM_X(m) = PERM_X(m) + (X(m,a)*charge(a))
          PERM_Y(m) = PERM_Y(m) + (Y(m,a)*charge(a)) 
          PERM_Z(m) = PERM_Z(m) + (Z(m,a)*charge(a))  
       ENDDO
       TOT_X(m) = PERM_X(m) + SOMME_MUX(m)
       TOT_Y(m) = PERM_Y(m) + SOMME_MUY(m)
       TOT_Z(m) = PERM_Z(m) + SOMME_MUZ(m)
    ENDDO
    VOLUME = boxx*boxy*boxz
    TOT_XX = sum(TOT_X)/VOLUME
    TOT_YY = sum(TOT_Y)/VOLUME
    TOT_ZZ = sum(TOT_Z)/VOLUME

    !-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    !-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    ALLOCATE(step(nb_steps))
    ALLOCATE(PX(nb_steps))
    ALLOCATE(PY(nb_steps))
    ALLOCATE(PZ(nb_steps))
    ALLOCATE(TOT_P(nb_steps))

    READ(13,*) 
    READ(13,*)
    READ(13,*)
    READ(13,*) step(1), PX(1), PY(1), PZ(1)
    DIFF_X = (PX(1)-TOT_XX)*conv
    DIFF_Y = (PY(1)-TOT_YY)*conv
    DIFF_Z = (PZ(1)-TOT_ZZ)*conv 
    DO w = 2, nb_steps
       READ(13,*) step(w), PX(w), PY(w), PZ(w)
    ENDDO

    DO w = 1, nb_steps
       PX(w) = (PX(w)*conv) - DIFF_X
       PY(w) = (PY(w)*conv) - DIFF_Y
       PZ(w) = (PZ(w)*conv) - DIFF_Z
       TOT_P(w) = sqrt((PX(w)*PX(w)) + (PY(w)*PY(w)) + (PZ(w)*PZ(w)))
       WRITE(20,*) TOT_P(w)
    ENDDO

END PROGRAM
