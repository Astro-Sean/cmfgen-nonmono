!
! Subroutine designed to facilitae (or check for inconsistencies) in
! the use of CMFGEN.
!
! Routine indicataes to the user whether ions can be deleted from the
! model, or whether it may be necesary to add additional high
! ionization stages.
!
	SUBROUTINE CHECK_IONS_PRESENT(ND)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Aleterd: 28-Jul-2024 - Added data for copper and zinc.
! Altered: 01-May-2024 - Line length of sata statements reduced to 132 for compatability with intel compiler.
! Created: 28-Jan-2009.
!
	INTEGER ND
	REAL(KIND=LDP) MAX_RATIO(NUM_IONS)
!
	REAL(KIND=LDP), PARAMETER ::  LOW_LIMIT=1.0E-10_LDP
	REAL(KIND=LDP), PARAMETER :: HIGH_LIMIT=1.0E-04_LDP
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) HDKT_EV
	REAL(KIND=LDP) T_VAL
	REAL(KIND=LDP) ED_VAL
	REAL(KIND=LDP) ION_FRAC
	REAL(KIND=LDP) INTERNAL_ENERGY
	REAL(KIND=LDP) EXTRA_INT_EN
	REAL(KIND=LDP) ATOM_ION_DENSITY
	INTEGER ISPEC
	INTEGER J,K
	INTEGER ID
	INTEGER DPTH_INDX
	INTEGER LOC_NUM_SPECIES
	LOGICAL FIRST_TIME
!
	INTEGER LUER,ERROR_LU,LUWARN,WARNING_LU
	EXTERNAL ERROR_LU,WARNING_LU
!
! The following table gives IP in eV. Grabbed off the web.
! As its only used to guide the user, the IP's need not be very accurate.
!
	INTEGER, PARAMETER :: NUM_IP_SPECIES=30
	INTEGER, PARAMETER :: NUM_IP_IONS=20
	REAL(KIND=LDP) IP(NUM_IP_IONS,NUM_IP_SPECIES)
	DATA IP(1:1,1)/    13.60_LDP/
	DATA IP(1:2,2)/    24.59_LDP,  54.42_LDP/
	DATA IP(1:3,3)/     5.39_LDP,  75.64_LDP, 122.45_LDP/
	DATA IP(1:4,4)/     9.32_LDP,  18.21_LDP, 153.90_LDP, 217.72_LDP/
	DATA IP(1:5,5)/     8.30_LDP,  25.16_LDP,  37.93_LDP, 259.37_LDP, 340.22_LDP/
	DATA IP(1:6,6)/    11.26_LDP,  24.38_LDP,  47.89_LDP,  64.49_LDP, 392.09_LDP, 489.99_LDP/
	DATA IP(1:7,7)/    14.53_LDP,  29.60_LDP,  47.45_LDP,  77.47_LDP,  97.89_LDP, 552.07_LDP, 667.04_LDP/
	DATA IP(1:8,8)/    13.62_LDP,  35.12_LDP,  54.94_LDP,  77.41_LDP, 113.90_LDP, 138.12_LDP, 739.28_LDP, 871.41_LDP/
	DATA IP(1:9,9)/    17.42_LDP,  34.97_LDP,  62.71_LDP,  87.14_LDP, 114.24_LDP, 157.16_LDP, 185.19_LDP, 953.91_LDP,1103.11_LDP/
	DATA IP(1:10,10)/  21.56_LDP,  40.96_LDP,  63.45_LDP,  97.12_LDP, 126.21_LDP, 157.93_LDP, 207.27_LDP, 
	1                  239.10_LDP,1195.82_LDP,1362.20_LDP/
	DATA IP(1:11,11)/   5.14_LDP,  47.28_LDP,  71.62_LDP,  98.91_LDP, 138.40_LDP, 172.18_LDP, 208.50_LDP, 264.25_LDP, 
	1                   299.86_LDP,1465.11_LDP,1648.71_LDP/
	DATA IP(1:12,12)/   7.65_LDP,  15.04_LDP,  80.14_LDP, 109.27_LDP, 141.26_LDP, 186.76_LDP, 225.02_LDP, 
	1                  265.96_LDP, 328.06_LDP, 367.50_LDP,1761.80_LDP,  1962.66_LDP/
	DATA IP(1:13,13)/   5.99_LDP,  18.83_LDP,  28.45_LDP, 119.99_LDP, 153.83_LDP, 190.48_LDP, 241.76_LDP, 284.65_LDP, 
	1                   330.13_LDP, 398.74_LDP, 442.00_LDP, 2315.02_LDP,2304.14_LDP/
	DATA IP(1:14,14)/   8.15_LDP,  16.35_LDP,  33.49_LDP,  45.14_LDP, 166.77_LDP, 205.26_LDP, 246.46_LDP, 303.54_LDP, 
	1                   351.12_LDP, 401.37_LDP, 476.36_LDP, 523.42_LDP,2437.63_LDP,2673.18_LDP/
	DATA IP(1:15,15)/  10.49_LDP,  19.76_LDP,  30.20_LDP,  51.44_LDP,  65.02_LDP, 220.42_LDP, 263.57_LDP, 309.60_LDP, 
	1                 372.13_LDP, 424.42_LDP, 479.46_LDP, 560.81_LDP, 611.74_LDP,2816.91_LDP,3069.84_LDP/
	DATA IP(1:16,16)/  10.36_LDP,  23.34_LDP,  34.79_LDP,  47.22_LDP,  72.59_LDP,  88.05_LDP, 280.94_LDP, 328.74_LDP, 
	1                 379.55_LDP, 447.50_LDP, 504.84_LDP, 564.44_LDP, 652.22_LDP, 707.01_LDP,3223.78_LDP,3494.19_LDP/
	DATA IP(1:17,17)/  12.97_LDP,  23.82_LDP,  39.61_LDP,  53.47_LDP,  67.80_LDP,  97.03_LDP, 114.19_LDP, 
	1                 348.28_LDP, 400.06_LDP, 455.62_LDP, 529.28_LDP,
	1                 592.00_LDP, 656.71_LDP, 749.76_LDP, 809.40_LDP,3658.52_LDP,3946.30_LDP/
	DATA IP(1:18,18)/  15.76_LDP,  27.63_LDP,  40.74_LDP,  59.81_LDP,  75.02_LDP,  91.01_LDP, 124.32_LDP, 
	1                 143.46_LDP, 422.45_LDP, 478.68_LDP, 538.96_LDP,
	1                 618.26_LDP, 686.10_LDP, 755.74_LDP, 854.77_LDP, 918.03_LDP,4120.88_LDP,4426.22_LDP/
	DATA IP(1:19,19)/   4.34_LDP,  31.63_LDP,  45.81_LDP,  60.91_LDP,  82.66_LDP,  99.39_LDP, 117.56_LDP, 
	1                 154.88_LDP, 175.82_LDP, 503.81_LDP, 564.75_LDP,
	1                 629.42_LDP, 714.62_LDP, 786.65_LDP, 861.06_LDP, 968.02_LDP,1033.42_LDP,4610.85_LDP,4934.04_LDP/
	DATA IP(1:20,20)/   6.11_LDP,  11.87_LDP,  50.91_LDP,  67.27_LDP,  84.50_LDP, 108.78_LDP, 127.17_LDP, 
	1                 147.23_LDP, 188.54_LDP, 211.28_LDP, 591.90_LDP, 657.20_LDP, 726.64_LDP, 817.64_LDP, 894.54_LDP, 
	1                 974.24_LDP,1087.21_LDP,1157.80_LDP,5128.76_LDP,5469.86_LDP/
	DATA IP(1:20,21)/   6.56_LDP,  12.80_LDP,  24.76_LDP,  73.49_LDP,  91.65_LDP, 110.68_LDP, 137.95_LDP, 158.06_LDP, 
	1                  180.03_LDP, 225.17_LDP, 249.80_LDP, 687.36_LDP, 756.69_LDP, 830.80_LDP, 927.50_LDP,1009.48_LDP,
	1                  1094.47_LDP,1212.62_LDP,1287.97_LDP,5674.75_LDP/
	DATA IP(1:20,22)/   6.83_LDP,  13.58_LDP,  27.49_LDP,  43.27_LDP,  99.30_LDP, 119.53_LDP, 140.85_LDP, 170.39_LDP, 
	1                   192.05_LDP, 215.92_LDP, 265.07_LDP, 291.49_LDP, 787.84_LDP, 863.14_LDP, 941.90_LDP,1043.68_LDP,
	1                   1130.74_LDP,1220.91_LDP,1346.32_LDP,1425.40_LDP/
	DATA IP(1:20,23)/   6.75_LDP,  14.66_LDP,  29.33_LDP,  46.71_LDP,  65.28_LDP, 128.13_LDP, 150.59_LDP, 173.39_LDP, 
	1                   205.83_LDP, 230.50_LDP, 255.69_LDP, 308.13_LDP, 336.28_LDP, 895.99_LDP, 976.00_LDP,1060.26_LDP,
	1                   1168.05_LDP,1260.29_LDP,1354.61_LDP,1486.24_LDP/
	DATA IP(1:20,24)/   6.77_LDP,  16.49_LDP,  30.96_LDP,  49.16_LDP,  69.46_LDP,  90.63_LDP, 160.18_LDP, 184.69_LDP, 
	1                   209.25_LDP, 244.39_LDP, 270.82_LDP, 297.97_LDP, 354.77_LDP, 384.16_LDP,1010.62_LDP,1096.54_LDP,
	1                   1184.64_LDP,1298.64_LDP,1396.07_LDP,1495.56_LDP/
	DATA IP(1:20,25)/   7.43_LDP,  15.64_LDP,  33.66_LDP,  51.20_LDP,  72.45_LDP,  95.56_LDP, 119.19_LDP, 194.54_LDP, 
	1                   221.80_LDP, 248.33_LDP, 285.95_LDP, 314.35_LDP, 343.58_LDP, 402.96_LDP, 435.16_LDP,1134.68_LDP,
	1                  1224.02_LDP,1317.30_LDP,1436.49_LDP,1539.09_LDP/
	DATA IP(1:20,26)/   7.90_LDP,  16.19_LDP,  30.65_LDP,  54.83_LDP,  75.04_LDP,  99.08_LDP, 124.99_LDP, 151.11_LDP, 
	1                 233.61_LDP, 262.11_LDP, 290.20_LDP, 330.83_LDP, 360.99_LDP, 392.18_LDP, 457.06_LDP, 489.26_LDP,
	1                1266.51_LDP,1357.72_LDP,1456.18_LDP,1581.59_LDP/
	DATA IP(1:20,27)/   7.88_LDP,  17.08_LDP,  33.50_LDP,  51.30_LDP,  79.49_LDP, 101.98_LDP, 128.93_LDP, 157.85_LDP, 
	1                 186.13_LDP, 275.38_LDP, 304.71_LDP, 335.80_LDP, 379.33_LDP, 411.46_LDP, 443.59_LDP, 511.95_LDP, 
	1                 546.58_LDP,1397.21_LDP,1504.58_LDP,1603.35_LDP/
	DATA IP(1:20,28)/   7.64_LDP,  18.17_LDP,  35.19_LDP,  54.93_LDP,  76.06_LDP, 107.79_LDP, 132.66_LDP, 161.68_LDP, 
	1                 192.78_LDP, 224.59_LDP, 320.98_LDP, 352.38_LDP, 384.51_LDP, 430.12_LDP, 464.32_LDP, 498.52_LDP, 
	1                 571.08_LDP, 607.03_LDP,1541.17_LDP,1647.92_LDP/
	DATA IP(1:20,29)/   7.73_LDP,  20.29_LDP,  36.84_LDP,  57.38_LDP,   79.8_LDP,  103.0_LDP,  139.0_LDP,  166.0_LDP,
	1                 198.0_LDP,  232.2_LDP,  265.33_LDP,  367.0_LDP,  401.0_LDP,  436.0_LDP,   483.1_LDP, 518.7_LDP,
	1                 552.8_LDP,  632.5_LDP,  670.608_LDP,1690.5_LDP/
	DATA IP(1:20,30)/  9.39_LDP,   17.96_LDP,  39.72_LDP,   59.57_LDP,  82.6_LDP,  108.0_LDP,   133.9_LDP,  173.9_LDP,
	1                 203.0_LDP,   238.0_LDP,  274.4_LDP,   310.8_LDP,  417.6_LDP,  453.4_LDP,  490.6_LDP,  540.0_LDP,
	1                 577.8_LDP,   613.3_LDP,  697.5_LDP,  737.37_LDP/

!
	MAX_RATIO=1.0_LDP
	LUER=ERROR_LU()
	LUWARN=WARNING_LU( )
	LOC_NUM_SPECIES=NUM_SPECIES
!
! Check data is available for all species.
!
	DO ISPEC=1,NUM_SPECIES
	  K=NINT(AT_NO(ISPEC))
!	  WRITE(LUER,*)ISPEC,K
	  IF(SPECIES_PRES(ISPEC) .AND. K .GT. NUM_IP_SPECIES)THEN
	    WRITE(LUER,*)'Warning in CHECK_IONS_PRESENT'
	    WRITE(LUER,*)'Data for atomic species with atomic no not availabe'
	    WRITE(LUER,*)'Requested atomic No is ',K
	    WRITE(LUER,*)'Maximum atomic No is ',NUM_IP_SPECIES
	    LOC_NUM_SPECIES=ISPEC-1
	    FLUSH(UNIT=LUER)
	    EXIT 
	  END IF
	END DO
	FLUSH(UNIT=LUER)
!
! Determine ionzation fractions.
!
	DO ISPEC=1,LOC_NUM_SPECIES
	  IF(SPECIES_PRES(ISPEC))THEN
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      T1=0.0_LDP
	      DO K=1,ND
	        T2=SUM(ATM(ID)%XzV_F(:,K))
	        T1=MAX(T1,T2/POP_SPECIES(K,ISPEC))
	      END DO
	      MAX_RATIO(ID)=T1
	    END DO
	    ID=SPECIES_END_ID(ISPEC)-1
	    MAX_RATIO(ID+1)=MAXVAL(ATM(ID)%DXzV_F/POP_SPECIES(:,ISPEC))	
	  END IF
	END DO
!
! Check whether the lowest ionization stages might be omitted.
!
	WRITE(LUWARN,'(/,/,A,A)')' Checking whether some low ionization stages ',
	1                        'may be omitted from the model.'
	DO ISPEC=1,LOC_NUM_SPECIES
	  IF(SPECIES_PRES(ISPEC))THEN
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	      IF(MAX_RATIO(ID) .GT. LOW_LIMIT)THEN
	        DO K=SPECIES_BEG_ID(ISPEC), ID-1
	          WRITE(LUWARN,'(3X,A,T11,A,ES8.2,A,ES8.2)')TRIM(ION_ID(K)),
	1           ' may not need to be include in model. Its maximum fractional abundance of ',
	1           MAX_RATIO(K),' < ',LOW_LIMIT
	        END DO
	        EXIT
	      END IF
	    END DO
	  END IF
	END DO
!
! Check whether the highest ionization stages might be omitted.
!
	WRITE(LUWARN,'(A,A)')' Checking whether some high ionization stages ',
	1                        'may be omitted from the model.'
	DO ISPEC=1,LOC_NUM_SPECIES
	  IF(SPECIES_PRES(ISPEC))THEN
	    DO ID=SPECIES_END_ID(ISPEC),SPECIES_BEG_ID(ISPEC),-1
	      IF(MAX_RATIO(ID) .GT. HIGH_LIMIT)THEN
	        DO K=ID+1,SPECIES_END_ID(ISPEC)
	          IF(K .EQ. SPECIES_END_ID(ISPEC))THEN
	            J=NINT(ATM(K-1)%ZXzV)+1
	          ELSE
	            J=NINT(ATM(K)%ZXzV)
	          END IF
	          WRITE(LUWARN,'(3X,A,T11,A,ES8.2,A,ES8.2)')TRIM(SPECIES_ABR(ISPEC))//TRIM(GEN_ION_ID(J)),
	1             'may not need to be include in model. Its maximum fractional abundance is ',
	1            MAX_RATIO(K),' < ',HIGH_LIMIT
	        END DO
	        EXIT
	      END IF
	    END DO
	  END IF
	END DO
!
! Check to see if lower ionization species need to be included.
!
	WRITE(LUWARN,'(A)')' '
	WRITE(LUWARN,'(A,A)')' Checking whether lower ionization stages ',
	1                             'need to be included in model'
	DO ISPEC=1,LOC_NUM_SPECIES
	  ID=1
	  K=NINT(AT_NO(ISPEC))
	  IF(SPECIES_PRES(ISPEC))ID=SPECIES_BEG_ID(ISPEC)
	  IF(SPECIES_PRES(ISPEC) .AND. ATM(ID)%ZXzV .GT. 1.1_LDP)THEN
	    IF(MAX_RATIO(ID) .GT. 0.1_LDP)THEN
	      J=NINT(ATM(ID)%ZXzV)
	      WRITE(LUWARN,'(3X,A,T11,A,A,A,ES8.1,A,F7.2)')TRIM(SPECIES_ABR(ISPEC))//TRIM(GEN_ION_ID(J)),
	1              'may need to be included in the model (Maximum ',TRIM(ION_ID(ID)),
	1              ' ionization fraction is',MAX_RATIO(ID),');   IP=',IP(J-1,K)
	    ELSE
	      EXIT
	    END IF
	  END IF
	END DO
!
! Now do a quick and dirty check to see whether additional ionization stages should be added.
! We first get the maximum temperature (in SN models this may not be at the inner
! boundary).
!
	T_VAL=0.0_LDP
! Default to the inner boundary: if no T(K) exceeds the initial T_VAL
! (e.g. an all-zero or negative T array), DPTH_INDX and ED_VAL would
! otherwise be left undefined for the ionization checks below.
	DPTH_INDX=ND
	ED_VAL=ED(ND)
	DO K=1,ND
	  IF(T(K) .GT. T_VAL)THEN
	    T_VAL=T(K)
	    ED_VAL=ED(K)
	    DPTH_INDX=K
	  END IF
	END DO
!
! To check whether additional ionization stages need to be included,
! we allow for the possibility that the previous ionization stage
! might not be dominant. Thus we multipy by ION_FRAC when its less
! than 1.
!
	HDKT_EV=4.7994145_LDP/4.13566733_LDP
	FIRST_TIME=.TRUE.
	EXTRA_INT_EN=0.0_LDP
	DO ISPEC=1,LOC_NUM_SPECIES
	  IF(SPECIES_PRES(ISPEC))THEN
	    ID=SPECIES_END_ID(ISPEC)-1
	    J=ATM(ID)%ZXzV+2
	    K=NINT(AT_NO(ISPEC))
	    ION_FRAC=1.0_LDP
	    DO ID=J,MIN(K,NUM_IP_IONS)
!	       WRITE(LUER,*)ISPEC,K,J,ION_FRAC,T_VAL,IP(ID,K); FLUSH(LUER)
	       ION_FRAC=ION_FRAC*(T_VAL**1.5_LDP)*EXP(-HDKT_EV*IP(ID,K)/T_VAL)/2.07E-22_LDP/ED_VAL
	       IF(ION_FRAC .GT. HIGH_LIMIT)THEN
	         IF(FIRST_TIME)THEN
	           WRITE(LUWARN,'(A)')' '
	           WRITE(LUWARN,'(A,A)')' Checking whether additional higher ionization stages ',
	1                             ' need to be included in the model.'
	           WRITE(LUWARN,'(A,A,/,A)')' NB: If XzV needs to be included it will also be necessary to',
	1                             ' include XzIV, as this was',
	1                             '      only included as the ground state.'
	           WRITE(LUWARN,'(A,1X,A,I3,A,ES10.4,A,ES10.4)')' Parameters at check depth:',
	1	          'Depth=',DPTH_INDX,'   T=',T_VAL,'   ED=',ED_VAL
	           FIRST_TIME=.FALSE.
	         END IF
	         WRITE(LUWARN,'(3X,A,T11,A,ES8.1,A)')TRIM(SPECIES_ABR(ISPEC))//TRIM(GEN_ION_ID(ID)),
	1              'may need to be included in the model (IF ~',ION_FRAC,')'
	         IF(ION_FRAC .GT. 1.0_LDP)ION_FRAC=1.0_LDP
	         EXTRA_INT_EN=EXTRA_INT_EN+ION_FRAC*POP_SPECIES(DPTH_INDX,ISPEC)*IP(ID,K)
	       ELSE
	         EXIT
	       END IF
	    END DO
	  END IF
	END DO
!
	ATOM_ION_DENSITY=0.0_LDP
	DO ISPEC=1,NUM_SPECIES
	  IF(SPECIES_PRES(ISPEC))THEN
	    ATOM_ION_DENSITY=ATOM_ION_DENSITY+POP_SPECIES(DPTH_INDX,ISPEC)
	  END IF
	END DO
!
	INTERNAL_ENERGY=0.0_LDP
	DO ISPEC=1,LOC_NUM_SPECIES
	  IF(SPECIES_PRES(ISPEC))THEN
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      T1=0.0_LDP
	      K=DPTH_INDX
	      T2=SUM(ATM(ID)%XzV_F(:,K))
	      J=ATM(ID)%ZXzV-1
	      IF(J .NE. 0)INTERNAL_ENERGY=INTERNAL_ENERGY+T2*SUM(IP(1:J,NINT(AT_NO(ISPEC))))
	    END DO
	    J=J+1
	    ID=SPECIES_END_ID(ISPEC)-1
	    INTERNAL_ENERGY=INTERNAL_ENERGY+ATM(ID)%DXzV_F(DPTH_INDX)*SUM(IP(1:J,NINT(AT_NO(ISPEC))))
	  END IF
	END DO
!
	T1=1.5_LDP*0.86173_LDP*T(DPTH_INDX)*(ED(DPTH_INDX)+ATOM_ION_DENSITY)
	WRITE(LUWARN,'(A)')' '
	WRITE(LUWARN,'(A,ES10.3,A)')'                   Temperature (10^4 K) is ',T(DPTH_INDX)
	WRITE(LUWARN,'(A,ES10.3,A)')'           Atom/ion density (per cm^3)) is ',ATOM_ION_DENSITY
	WRITE(LUWARN,'(A,ES10.3,A)')'                 Electrons per atom/ion is ',ED(DPTH_INDX)/ATOM_ION_DENSITY
	WRITE(LUWARN,'(A,ES10.3,A)')' Total thermal kinetic energy (ev/atom) is ',T1/ATOM_ION_DENSITY
	WRITE(LUWARN,'(A,ES10.3,A)')'  Approximate internal energy (ev/atom) is ',INTERNAL_ENERGY/ATOM_ION_DENSITY
	WRITE(LUWARN,'(A,ES10.3,A)')' App. missing internal energy (ev/atom) is ',EXTRA_INT_EN/ATOM_ION_DENSITY
!
! 4/c . sigma T^4 conerted to eV  (NB: Int J = sigma T^4 / pi).
!
	T1=4.0_LDP*3.14159265_LDP*5.670400E-05_LDP*1.0E+16_LDP/1.60217733E-12_LDP/2.99792458E+10_LDP
	WRITE(LUWARN,'(A,ES10.3,A)')'             Radiation energy (ev/atom) is ',(T1*T(DPTH_INDX)**4)/ATOM_ION_DENSITY
	FLUSH(UNIT=LUER); FLUSH(UNIT=LUWARN)
!
	RETURN
	END
