      program lee   
!    Lee code program for Dense Plasma Focus simulation for GNU ForTran
!       
!    Copyright (C) 2020 Jorge A. García Gallardo
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.      
     implicit none
!
!
! NOTES :	1. Only H, D or T calculations possible, other elements removed
!		2. Tapper capabilities removed.
!		3. Variable names were kept when possible
! 		4. Line labels kept as comments when possible
!		5. Original comments kept when possible
!

	integer, parameter :: dp = SELECTED_REAL_KIND(15)
	character(len=32) :: inputfile
	character(len=32) :: outFile, outParamFile, outRadFile, outtics, outEnergy, outTemp
	character*1 tab
	! constant values
	real(dp):: MU, Pi,bc,mi,MUK
	real(dp):: CON11, CON12, CON2
	real(dp):: FRF, FE
	! tube
	real(dp)::RADA,RADB,Z0,L0,C0,V0,R0
	! gas
	real(dp)::G,MW,P0
	integer::ZN
	real(dp)::dissociatenumber
	! model
	real(dp)::massf,massfr,currf,currfr
	! control
	integer(dp)::Ncon
	real(dp)::D
	character(len=80)::title, description
	! namelist
	namelist/bank/L0,C0,R0
	namelist/tube/RADA,RADB,Z0,V0
	namelist/gas/G,P0,ZN,MW,dissociatenumber
	namelist/model/massf,massfr,currf,currfr
	namelist/control/Ncon, D, title, description
	! other variables
	real(dp)::G1, G2, GCAP
	real(dp)::N0,C,f
	real(dp)::RHO
	real(dp)::I0,T0, TA, ZZCHAR, AL, AA, RESF, BE, BF
	real(dp)::VPINCHCH,TPINCHCH
	real(dp)::ALT
	! 1st phase variables
	real(dp)::T,z,ZZ,ZG
	real(dp)::II, I, V, Io
	real(dp)::AC
	real(dp)::Ipeak,ZZRpeak
	real(dp)::TR,VR,IR
	real(dp)::ZZR,ACR,ZR
	real(dp)::IIR,I0R
	! 2nd phase variables
	real(dp)::CFR, CurrentFactorRatio
	real(dp):: AAg
	integer(dp) :: rowi, FIRSTRADIALROW
	integer(dp) :: zFLAG, gFLAG
	real(dp):: DREAL
	real(dp):: IDELAY, KPDELAY, KSDELAY
	real(dp):: VS, VSR
	real(dp):: SDSPEED, SDDELAYTIME, backhowmanyrows
	integer(dp):: BACKROWNUMBER, delayrow
	real(dp), dimension(1:10000) :: IR_i, KSR_i, KPR_i, SSR_i
	real(dp) :: DR, E1, E2, E3, E4, EINP
	real(dp) :: K1, KP, KPR, KS, KSR
	real(dp) :: peakvp, peakvs, spr, ssr, szr
	real(dp) :: TeV, TM, TRRadial, TRRadialStart
	real(dp) :: VP, VRMax, vsdelay
	real(dp) :: VZ, ZF, ZFR
	! extra variables for v6
	real(dp) :: Dratio, gratio, Pratio, TMmax, Tratio
	real(dp) :: QLN1, QSXR, rbr, rln, rrec, rsxr, Tpinchmin
	real(dp) :: nimax, nimax1, ni1
	! 3rd phase variables
	real(dp) :: CH, IDOT, IDOTKAUS, MUP
	real(dp) :: rp, RRF, RRFMM, RS, TMRS
	real(dp) :: VRF, VRFCMUS, VSV, tc5

	! 4th phase variables
	real(dp):: AB, PRAD
	real(dp):: amin, correctfactor, e5, fbacktrace
	integer(dp):: FLAG, ENDFLAG
	real(dp) :: IDOT1, IDOT2, IDOT3, IDOT4, IDOT5
	real(dp) :: Ipinch
	real(dp) :: NBN, Ni, nipinch, NN, NTN, NTNDOT
	real(dp) :: PBB, PBR, PJ, PLN, PM, PRADS, PREC
	real(dp) :: Q, QBR, QDOT, QJ, QLN, QRAD, QREC, QTOTAL
	real(dp) :: R, ratiobr, ratioln, ratiorec, rpOLD, RPSTART, SDSPEEDSTART
	real(dp) :: sfactor, sflag1, sflag2, sflag3, sigv
	real(dp) :: tevstart, tmb, tpinch, trad, trad1, tstart
	real(dp) :: zfdot, zmax
	! 5th phase variables
	real(dp) :: Ec, Ecap, Edip, EL0, ELp, ELt, EMAGp, Excircuit
	real(dp) :: H, ID, Kmin, L, L0t, L1, Lp, M, MAG
	real(dp) :: SFI0, SFIdip, SFIpeak, SFIpinch, sig, TSlowcompressionphase, ZS

	! read file name from imput argument
	call get_command_argument(1, inputfile)
	open(unit=1, file=trim(inputfile), status='old')
	outParamFile =trim(inputfile(1:index(inputfile,".in")))//"param.out"	
	open(unit=2, file=trim(outParamFile), status='replace')
	outFile =trim(inputfile(1:index(inputfile,".in")))//"out"
	open(unit=3, file=trim(outFile), status='replace')

	outRadFile=trim(inputfile(1:index(inputfile,".in")))//"rad.out"
	open(unit=4, file=trim(outRadFile), status='replace')

	outtics=trim(inputfile(1:index(inputfile,".in")))//"tics.out"
	open(unit=5, file=trim(outtics), status='replace')
	tab = char(9)

	outEnergy=trim(inputfile(1:index(inputfile,".in")))//"energy.out"
	open(unit=10, file=trim(outEnergy), status='replace')


	outTemp=trim(inputfile(1:index(inputfile,".in")))//"temp.out"
	open(unit=11, file=trim(outTemp), status='replace')


	read(1,bank)
	read(1,tube)
	read(1,gas)
	read(1,model)
	read(1,control)
	close(1)
	MU = 1.257D-06
	Pi = 3.142
	bc = 1.38D-23
	mi = 1.67D-27
	MUK = MU / (8.0 * Pi * Pi * bc)
	CON11 = 1.6D-20
	CON12 = 9.999999D-21
	CON2 = 4.6D-31

	FRF = 0.3
	FE = 1 / 3.0
	EINP = 0.0
	FLAG = 0


	ENDFLAG = 0
	NTN = 0.0
	NBN = 0.0
	NN = 0.0
	VRMAX = 0.0
	G1 = 2 / (g + 1)
	G2 = (g - 1) / g
	GCAP = (g + 1) / (g - 1)

	! Calculate ambient number density and some ratios
	N0 = 2.69D+25 * P0 / 760.0

	C = RADB / RADA
	f = z0 / RADA

	! Convert to SI units
	R0 = R0 / 1000.0
	C0 = C0 * 1D-06
	L0 = L0 * 1D-09
	RADB = RADB * 0.01
	RADA = RADA * 0.01
	z0 = z0 * 1D-02
	V0 = V0 * 1000
	RHO = P0 * 2.33D-04 * MW / 4

	TM=0.0

	! Calculate characteristic quantities and scaling parameters
	I0 = V0 / (DSQRT(L0 / C0))
	T0 = DSQRT(L0 * C0)
	TA = DSQRT(4.0 * Pi*Pi *(C*C-1.0) / (MU*DLOG(C)) ) *( (z0*DSQRT(RHO))/(I0/RADA) ) * ((DSQRT(Massf))/CURRF)
	ZZCHAR = z0 / TA
	AL = T0 / TA ! alpha
	AA = DSQRT((g+1.0) * (C * C - 1.0) / DLog(C)) * (f / 2.0) *(DSQRT(massf / massfr)) * (currfr / CURRF)
	RESF = R0 / (DSQRT(L0 / C0))
	BE = 2.0 * (1.0D-07) * DLog(C) * z0 * CURRF / L0
	BF = BE / (DLog(C) * f)
	VPINCHCH = ZZCHAR * AA / f
	TPINCHCH = RADA / VPINCHCH

	! Calculate ratio of characteristic capacitor time to sum of characteristic axial & radial times
	ALT = (AL * AA) / (1 + AA)

	write(2,'("*** The LEE code for plasma focus calculations, CNEA (c) 2019-2020 ***")')
	write(2,'(" ",A)') title
	write(2,'(" ",A/"*** Input data ***")') description
	write(2,'("- Geometry data:")')
	write(2,'("a =",ES14.7," [m];",X,"b =",ES14.7," [m];",X,"zo =",ES14.7," [m]"/)') RADA,RADB,z0
	write(2,'( "- Electrical data:")')
	write(2,'("Co =",ES14.7," [F];",X,"Lo =",ES14.7," [Hy];",X,"Vo =",ES14.7," [V]"/)') C0,L0,V0
	write(2,'("- Gas data:")')
	write(2,'("Po =",ES14.7," [Torr];",X,"gamma =",ES14.7,";",X,"Z =",I3,";"/)') P0,g,ZN
	write(2,'("*** General parameters ***")')
	write(2,'("rho =",ES14.7," [Kg/m3]")') RHO
	write(2,'("- Electrical:")')
	write(2,'("Io =",ES14.7," [A];",X,"to =",ES14.7," [s]; ",X,"ta =",ES14.7," [s]"/)') I0,T0,TA
	write(2,'("va =",ES14.7," []"/)') ZZCHAR
	write(2,'("- Pinch parameters:")')
	write(2,'("Vpinch =",ES14.7," [V];",X,"Tpinch =",ES14.7," [s] "/)') VPINCHCH,TPINCHCH
	write(2,'("- Adimensional parameters:")')
	write(2,'("c = b/a =",ES14.7,";",X,"f = zo/a =",ES14.7,"; ",X,"alpha = to/ta =",ES14.7)') C, f,AL
	write(2,'("alpha_1 =",ES14.7,";",X,"beta =",ES14.7,";",X,"beta_1 =",ES14.7, ";",X,"delta = ",ES14.7)') AA,BE,BF, RESF
	write(2,'("Alt = ",ES14.7,/)') ALT
	if (P0.ge.20) then
		write(2,*) "* WARNING! In real experiments, Pressure above 20 torr will not produce focussing"
	end if
	if (ZN.le.2) then
		if (P0.le.0.05) then
			write(2,*) "* WARNING! In real experiments, Pressure below 0.05 torr may not yield good focus"
		end if
		if (ALT.le.0.68) then
			write(2,*) "* WARNING! Total TRANSIT TIME (axial + radial) MAY BE TOO LONG COMPARED TO effective DISCHARGE Drive TIME"
		end if
		write(*,*) "starting phase 1..."
!		CALL flush(6)
!		CALL flush(0)
		write(2,'(/"*** Phase 1 (axial) parameters ***"/)')
! phase 1: axial phase

		T = 0.0_dp
		z = 0.0_dp
		ZZ = 0.0_dp
		II = 1.0_dp
		I = 0.0_dp
		Io = 0.0_dp
		V=0.0_dp
		AC = AL * DSQRT(1.0_dp / 2.0)
		Ipeak = 0.0_dp
		ZZRpeak = 0.0_dp
		EINP = 0.0_dp
		! Start numerical integration of AXIAL PHASE
		do while(z.lt.1.0)
			! 580
			T = T + D
			ZZ = ZZ + AC * D
			z = z + ZZ * D
			I = I + II * D
			Io = Io + I * D
			! Convert data to real, but convenient units
			! 671
			TR = T * T0 * 1.0D+06
			VR = V * V0 * 1.0D-03
			IR = I * I0 * 1.0D-03
			ZZR = (ZZCHAR / AL) * ZZ * 1.0D-04
			ACR = ((ZZCHAR / AL) / T0) * AC * 1.0D-10
			ZR = z * z0 * 100
			IIR = (II * I0 / T0) * 1.0D-09
			I0R = Io * I0 * T0

			if (IR.gt.Ipeak) then 
				Ipeak = IR
			end if
			if (ZZR.gt.ZZRpeak) then 
				ZZRpeak = ZZR
			end if
			DR = D * T0
			EINP = EINP + (1D-07) * (DLog(C)*ZZR *(1D+04) * IR * IR * (1D+06) * CURRF * CURRF) * DR

			! 680

			ni = massf * N0 * ((g + 1) / (g - 1))
			nimax = N0 * ((g + 1) / (g - 1))

			! Write, in convenient real units, the step-by-step data from numerical integration
			! 771
			write(3,'(9(ES14.7X))')TR,IR,VR,ZR,RADA,RADA,ZZR,0.0,0.0
!			write(4,'(8(ES14.7X))')TR,0.0,0.0,0.0,0.0,0.0,0.0,0.0
			write(10,'(2(ES14.7X))')TR,EINP
			write(11,'(2(ES14.7X))')TR,TM
			V =z * II + ZZ * I
			V = V * BE
			AC = (AL * AL * I * I - ZZ * ZZ) / z
			II = (1.0D+00 - Io - BE * ZZ * I - RESF * I) / (1.0D+00 + BE * z)
			! 800
		end do
		
		write(2,'("I(peak) = ",ES14.7X,"[KA]")') Ipeak
		write(2,'("v_z(peak) = ",ES14.7X,"[cm/us]")') ZZRpeak
		! 815
		ZG = ZZ ! record last value of axial speed


		write(2,'("v_z(last) = ",ES14.7X,"[cm/us]")') ZZR

		write(2,'("(*) PHASE 1 ends at = ",ES14.7X,"[us]")') TR
		write(5,'(I1A1ES14.7)')1,tab,TR
! phase 2 : radial phase 
		write(*,*) "starting phase 2..."
		write(2,'(/"*** Phase 2 (radial) parameters ***"/)')
		CurrentFactorRatio = currfr / CURRF
		CFR = CurrentFactorRatio
		BE = BE * CFR
		BF = BF * CFR
		CURRF = CURRF * CFR

		! Radial phase RI, distances are relative to radius a.
		! KS is shock position, KP is radial piston position, ZF is focus
		! pinch length, all normalized to inner radius a; VS and VP are
		! radial shock and piston speed,VZ is axial pinch length elongating rate
		! Distances, radius and speeds are relative to radius of anode.
		! AS BEFORE QUANTITIES WITH AN R ATTACHED HAVE BEEN RE-COMPUTED AS REAL,
		! i.e. UN-NORMALIZED QUANTITIES EXPRESSED IN USUAL LABORATORY UNITS.
		! FOR THIS PHASE Z=EFFECTIVE CHARGE NUMBER!!!
		rowi = 20

		FIRSTRADIALROW = rowi
		! Reset time increment to finer step-size
		KS = 1
		KP = 1
		ZF = 0.00001
		zFLAG = 0
		gFLAG = 0
		! SET TIME INCREMENT TO HAVE ABOUT 1500 (up to 3000 for high pressure)
		! STEPS IN RADIAL INWARD SHOCK PHASE
		! 950
		DREAL = TPINCHCH / 700
		D = DREAL / T0
		! Set initial 'LookBack' values, for compensation of finite small disturbance speed
		IDELAY = I
		KPDELAY = KP
		KSDELAY = KS
		VSDELAY = -1
		tc5 = 1
		! Set initial value, approximately, for CHARGE NUMBER Z
		! For H2,D2 and He, assume fully ionized with gamma=1.667 during all of radial phase
		z = ZN
		! Start Step-by-step integration of Radial Inward Shock Phase, in non-dimensional units
		! First, compute Inward shock speed
		trradialstart = T * T0
		write(2,'("radial phase starts at = ",ES14.7X,"[us]")') trradialstart * (1D+06)
		! 980
		do while(KS.gt.0.0)
			GCAP = (g + 1.0) / (g - 1.0)

			VS = -AL * AA * IDELAY / (KPDELAY)
			VSR = VS * RADA / T0
			! Real temperature is needed to DETERMINE SMALL DISTURBANCE SPEED
			! FOR COMMUNICATION TIME CORRECTION.
			! Hence the shock speed is re-calculated in SI units, then Plasma 
			! Temp TM is calculated, based on shock theory

			! 1000
			TM = (MW / (8315)) * ((GCAP - 1)/(GCAP * GCAP)) * ((VSR * VSR) / ((1+z) * dissociatenumber))
			TeV = TM / (1.16 * 1D+04)
			! Select Table for G & Z; according to gas
			! (for ZN = 1 or 2)
			! 1005
			! 1080
			G1 = 2 / (g + 1)
			G2 = (g - 1) / g
			if (FLAG.eq.10) then
			end if

			! Next compute Axial elongation speed and Piston speed, using 'lookback' values 
			! to correct for finite small disturbance speed effect

			! 2000
			VZ = -G1 * VS
			K1 = KS / KP

			E1 = G1 * K1 * VSDELAY

			E2 = (1 / g) * (KP / I) * (1 - K1 * K1) * II
			E3 = (G1 / 2) * (KP / ZF) * (1 - K1 * K1) * VZ
			E4 = G2 + (1 / g) * K1 * K1
			VP = (E1 - E2 - E3) / E4
			V = (BE - BF*(DLog(KP/C))* ZF) * II - I*(BF*(ZF/KP)*VP + (BF*(DLog(KP / C))) * VZ)
			II = (1 - Io + BF*I*(ZF/KP)*VP + BF*(DLog(KP/C)) *I*VZ - RESF*I) / (1 + BE - BF*(DLog(KP/C))*ZF)
			! Increment time and Integrate, by linear approx, for I, flowed charge I0, KS, KP and ZF
			T = T + D
			I = I + II * D
			Io = Io + I * D
			KS = KS + VS * D
			KP = KP + VP * D
			ZF = ZF + VZ * D
			! * Re-scales speeds, distances and time to real, convenient units
! 2210
			SSR = VS * (RADA / T0) * 1D-04
			SPR = VP * (RADA / T0) * 1D-04
			SZR = VZ * (RADA / T0) * 1D-04

			! in milimeteres:
			KSR = KS * RADA * 1000
			KPR = KP * RADA * 1000
			ZFR = ZF * RADA * 1000
! 2250
			TR = T * T0 * 1D+06
			VR = V * V0 * 1D-03
			IR = I * I0 * 1D-03
			IIR = (II * I0 / T0) * 1E-09
			DR = D * T0
			if (IR.gt.Ipeak) then
				Ipeak = IR
			end if
			! Obtain Max induced voltage
			if (VR.gt.VRMAX) then
				VRMAX = VR
			end if
			TRRadial = TR * (1D+03) - trradialstart * (1D+09)

			if (SSR.lt.peakvs) then
				peakvs = SSR
			end if
			if (SPR.lt.peakvp) then
				peakvp = SPR
			end if
			! Integrate to find EINP, energy into plasma ie work done by magnetic piston, 
			! by dynamic resistance effect
			EINP = EINP + (1D+03) * (SZR * DLog(1000*RADB / kpr) - (SPR * (zfr/kpr))) *IR*IR * CURRF * CURRF * DR

			ni = massfr * N0 * ((g + 1) / (g - 1))
			nimax = N0 * ((g + 1) / (g - 1))

			write(3,'(9(ES14.7X))')TR,IR,VR,ZFR/10+z0*100,KSR/1000,KPR/1000,SZR,SSR,SPR
			!write(4,'(8(ES14.7X))')TR,0.0,0.0,0.0,0.0,0.0,0.0,0.0
			write(10,'(2(ES14.7X))')TR,EINP
			write(11,'(2(ES14.7X))')TR,TM
			! store 'past' values for compute delayed variables
			IR_i(rowi) = IR
			KSR_i(rowi) = KSR
			KPR_i(rowi) = KPR
			SSR_i(rowi) = SSR
			rowi = rowi + 1

			! To apply finite small disturbance speed correction. Compute propagation 
			! time and the 'lookback' row number
			if (KSR.gt.KPR) then
				IDELAY = I
				KPDELAY = KP
				VSDELAY = VS
				KSDELAY = KS
			else
				SDSPEED = DSQRT(g * dissociatenumber * (1 + z) * bc * TM / (MW * mi))
				SDDELAYTIME = ((KPR - KSR) / 1000) / SDSPEED
				backhowmanyrows = SDDELAYTIME / DR
				BACKROWNUMBER = rowi - NINT(backhowmanyrows)
				if (BACKROWNUMBER.lt.FIRSTRADIALROW) then
					BACKROWNUMBER = FIRSTRADIALROW + 1
				end if
				delayrow = BACKROWNUMBER - 1
				if (delayrow.lt.20) then
					delayrow = 20
				end if
				! Look back to appropriate row to obtain 'lookback' quantities;
				! also non-dimensionalize these quantities
				IDELAY = IR_i(delayrow) / (I0 * 1D-03)
				KPDELAY = KPR_i(delayrow) / (RADA * 1000)
				VSDELAY = SSR_i(delayrow) / ((RADA / T0) * 1D-04)
				KSDELAY = KSR_i(delayrow) / (RADA * 1000)

			end if

			! Check whether inward shock front has reached axia
		! 2314
		end do
		! Put ni1 as the last average ion density on axis before reflected shock starts; nimax1 as last shocked density before RS starts
		ni1 = ni
		nimax1 = nimax
		PLN = 0
		! Inward shock front has reached axis, we have exited from Radial Inward Phase
		! and now go on to the Reflected Shock Phase

		! Put ni1 as the last average ion density on axis before reflected shock starts; nimax1 as last shocked density before RS starts
		write(2,'("EINP="ES14.7X" ")')EINP
		write(2,'("Final plasma Temp. =",ES14.7," [eV]; (",ES14.7," [K] )")')TeV, TeV*11604.45
		write(2,'("(*) PHASE 2 ends at = ",ES14.7X,"[us]")') TR
		write(5,'(I1A1ES14.7)')2,tab,TR
		write(2,'("Phase 2 duration = ",ES14.7X,"[us]")') TRRadial* 1D-03 !(TRRadial * 1D-03 - trradialstart) * 1D+03
! phase 3 : radial reflected shock phase 
		write(*,*) "starting phase 3..."
		write(2,'(/"*** Phase 3 (radial reflected shock) parameters ***"/)')
		VS = VS * RADA / T0
		RS = KS * RADA
		rp = KP * RADA
		ZF = ZF * RADA
		VZ = VZ * RADA / T0
		VP = VP * RADA / T0

		! 2510
		T = T * T0
		D = D * T0
		I = I * I0
		CH = IO * I0 * T0
		IDOT = II * I0 / T0
		! 2800
		! Rem VRF is reflected shock speed taken as a constant value at 0.3 of on-axis forward shock speed
		! Take Reflected shock temperature to be twice on-axis incident shock temperature
		! TMRS = 2 * TM

		! Take strong planar shock approximation (Ref: Robert Gross: The Physics of Strong Shock Waves in Gases 1969,
		! manuscript for Procs of International School of Physics "Enrico Fermi" Course XLVIII, High Energy Density,
		! Varenna, Italy; Academic Press.)
		! However we take RS speed as 0.3 of incident shock speed instead of 0.5
		! [for planar strong shock]as in R Gross; in order to account for diverging radial geometry
		
		gratio = (g + 1) / (g - 1)
		Pratio = 2 + gratio
		Tratio = Pratio * ((gratio + Pratio) / (1 + gratio * Pratio))
		Dratio = Pratio / Tratio
		TMRS = TM * Tratio
		TMmax = TMRS

		TeV = TMRS / (1.16D+04)


		FLAG = 10

		G1 = 2 / (g + 1)
		G2 = (g - 1) / g
		! 3030
		RRF = 0
		FRF = 0.33
		VSV = VS

		VRF = -FRF * VS
		G1 = 2 / (g + 1)
		G2 = (g - 1) / g
		MUP = MU / (2 * Pi)
		VZ = -G1 * VS

		do while(RRF.le.rp)
			! 3080
			T = T + D

			RRF = RRF + VRF * D
			VRFCMUS = VRF * 1E-04
			K1 = 0
			E1 = G1 * K1 * VSV
			E2 = (1 / g) * (rp / I) * (1 - K1 * K1) * IDOT
			E3 = (G1 / 2) * (rp / ZF) * (1 - K1 * K1) * VZ
			E4 = G2 + (1 / g) * K1 * K1

			VP = (E1 - E2 - E3) / E4
			IDOT = (V0 - (CH / C0) - I * R0 - I * CURRF * MUP * ((DLog(RADB / rp)) * VZ - (ZF / rp) * VP)) &
			 	/ (L0 + MUP * CURRF * ((DLog(C)) * z0 * tc5 + (DLog(RADB / rp)) * ZF))

			V = MUP * I * ((DLog(RADB / rp)) * VZ - (ZF / rp) * VP) &
				+ MUP * ((DLog(RADB / rp)) * ZF + (DLog(C)) * z0 * tc5) * IDOT
			V = V * CURRF
			I = I + IDOT * D
			CH = CH + I * D
			rp = rp + VP * D
			! 3240
			ZF = ZF + VZ * D
			! is time increment in secs, DKPR is piston position increment & DZFR length position
			! increment, both in SI units
			if (IR.gt.Ipeak) then
				Ipeak = IR
			end if
			! Obtain Max induced voltage
			if (VR.gt.VRMAX) then
				VRMAX = VR
			end if
			TRRadial = TR *1.0D+03 - trradialstart * (1.0D+09)
			if(SSR.lt.peakvs) then
				peakvs = SSR
			end if
			if(SPR.lt.peakvp) then
				peakvp = SPR
			end if


			! Convert to Real convenient units for print out
			TR = T * 1D+06
			VR = V * 1D-03
			IR = I * 1D-03
			KPR = rp * 1D+03
			ZFR = ZF * 1D+03
			SPR = VP * 1D-04
			SZR = VZ * 1D-04
			IDOTKAUS = IDOT * 1D-09
			! 3486
			RRFMM = RRF * 1D+03

			! Integrate to find EINP, energy dissipated by dynamic resistance effect,
			! which is 0.5 (Ldot) I^2, considering current taking part in the motion
			EINP = EINP + (1D-07) * (SZR * (1D+04) * DLog(1000 * RADB / kpr) &
			- (SPR * (1D+04) * (1000 / kpr) * (zfr / 1000))) * IR * IR * (1D+06) * CURRF * CURRF * D
			! Also integrate for piston work
			!Wpiston = Wpiston + 0.1 * (DZFR * Log(1000 * RADB / kpr) - DKPR * (zfr/kpr)) * IR*IR * CURRF*CURRF

			if (IR.gt.Ipeak) then
				Ipeak = IR
			end if
			! Determine max induced voltage for beam-gas neutron yield computation
			if (VR.gt.VRMAX) then
				VRMAX = VR
			end if
			! use Dratio from RS of strong shocks as described above
			ni = ni1 * Dratio
			nimax = nimax1 * Dratio
			write(3,'(9(ES14.7X))')TR,IR,VR, (ZFR/10+z0*100),RRF,rp,SZR,VRFCMUS,SPR
			!write(4,'(9(ES14.7X))')TR,TM,0.0,0.0,0.0,0.0,0.0,0.0,0.0
			write(10,'(2(ES14.7X))')TR,EINP
			write(11,'(2(ES14.7X))')TR,TM
		end do
		! RS HAS HIT PISTON. RS PHASE ENDS
		write(2,'("Radial phase duration = ",ES14.7X,"[us]")') TR *1.0D+03 - trradialstart * (1.0D+09)
		write(2,'("(*) PHASE 3 ends at = ",ES14.7X,"[us]")') TR
		write(5,'(I1A1ES14.7)')3,tab,TR
		write(2,'("EINP="ES14.7X" ")')EINP
! phase 4 : radiative phase 
		write(*,*) "starting phase 4..."
		write(2,'(/"*** Phase 4 (radiative phase) parameters ***"/)')
		! As RS hits piston, the pressure exerted by the doubly shocked column on the piston shoots
		! up by a factor of approx 6; this will slow the piston down further or even push it back.
		! This effect is included in the following section.
		! However, due to 2-D effect, the over-pressure may not be significant.
		sflag1 = 0
		sflag2 = 0
		sflag3 = 0
		RPSTART = rp
		TeVSTART = TeV
		SDSPEEDSTART = ((g * dissociatenumber * (1 + z) * bc * TM / (MW * mi)))**0.5
		TRAD1 = 0.5 * RPSTART / SDSPEEDSTART
		TRAD1 = 2.0 * TRAD1
		D = DREAL / 1D+08
		! 4002
		TStart = T
		Ipinch = I * CURRF / 1000

		QLN=0.0
		amin = KPR
		zmax = 0
		Tpinch = 0
		nipinch = 0
		write(2,'("- initial phase parameters:")')
		write(2,'("RPSTART = ",ES14.7X,"[]")') RPSTART
		write(2,'("initial plasma Temp. =",ES14.7," [eV]; (",ES14.7," [K] )")')TeVSTART, TeVSTART*11604.45
		write(2,'("SDSPEEDSTART = ",ES14.7X,"[]")') SDSPEEDSTART
		write(2,'("TRAD1 = ",ES14.7X,"[]")') TRAD1
		write(2,'("D = ",ES14.7X,"[]")') D
		write(2,'("TStart = ",ES14.7X,"[]")') TStart
		write(2,'("I_total = ",ES14.7X,"[KA]")') I/1000
		write(2,'("I_pinch = ",ES14.7X,"[KA]")') Ipinch
		write(2,'("amin = ",ES14.7X,"[cm]")') amin

		!4300
		! esto no está en la versión original pero debería estar
		sfactor = 0
		do while((ENDFLAG.ne.1).and.(ENDFLAG.ne.2))


			G1 = 2 / (g + 1)
			G2 = (g - 1) / g

			!Rem Compute Joule heating and radiation terms
			! 4400
			ni = N0 * fe * massfr * (RADA / rp) * (RADA / rp)

			! 4550
			TM = MUK * I * I * CURRF * CURRF / ((1 + z) *dissociatenumber* N0 * RADA * RADA * fe * massfr)

			TeV = TM / (1.16 * 1D+04)
			R = 1290 * z * ZF / (Pi * rp * rp * (TM**(1.5)))

			PJ = R * I * I * CURRF * CURRF
			PBR = -(CON11 * NI) * (TM**(0.5)) * (CON12 * NI) * Pi * (rp * rp) * ZF * (z**3.0)
			PREC = -5.92 * (1D-35) * NI * NI * (z**5.0) * Pi * (rp * rp) * ZF / (TM**(0.5))
			PLN = -(CON2 * NI) * NI * z * (ZN**4.0) * Pi * (rp * rp) * ZF / TM


			! Apply Plasma Self Absorption correction to PBR PREC and PLN:
			! PM is photonic excitation number; AB is absorption corrected factor
			! If AB<1/2.7183, Radiation goes from volume-like PRAD to surface-like PRADS;
			! PRADS has a limit being Blackbody Rad PBB

			PM = 1.66 * (1D-11) * (rp * 100) * (ZN**(0.5)) * (NI * (1D-06)) / (z * (TeV**(1.5)))

			AB = 1.0 + (((1D-14) * (NI * (1D-06)) * z) / (TeV**(3.5)))
			AB = 1.0 / AB
			AB = AB**(1.0 + PM)

			PBR = AB * PBR
			PREC = AB * PREC
			PLN = AB * PLN
			! surface-like emission:
			PRADS = -2.3 * (1D-15) * (ZN**3.5) * (z**0.5) * (TM**4) * 3.142 * rp * (2 * ZF)
			! calibration factor for neon (NX2); got to check for other machines and gases that at cross-over point
			! from volume to surface emission there is a smooth transition in power.
			PRADS = 0.032 * PRADS
			PBB = -5.7 * (1D-08) * (TM**4) * (3.142 * rp * (2 * ZF + rp))
		 	if(ZN.eq.1) then
				if (MW.eq.5) then
				! Rem for D-T (50:50), compute 1. thermonuclear neutron yield component; SIGV computed in m3sec-1
					!4650
					SIGV = 3.68 * (1D-18) * (DExp(-19.94 * (TeV / 1000)**(-(1 / 3)))) * (TeV / 1000)**(-(2.0/3.0))
					!4660
					NTNDOT = 0.5 * NI * NI * 3.142 * (rp * rp) * ZF * SIGV
					NTN = NTN + NTNDOT * D
				else if (MW.eq.2) then
					! GoTo 4700
				else
					! For deuterium, compute 1. thermonuclear neutron yield component:NTN
					! SIGV computed in m3sec-1
					if(TeV.lt.100) then
						!4623
						SIGV = 0
					else if(TeV.eq.100) then
						!4623
						SIGV = 0
					else if((100.lt.TeV).and.(TeV.lt.500)) then
						!4627
						SIGV = (1D-27) * (TeV / 1000.0)**10.0
						!4660
						NTNDOT = 0.5 * NI * NI * 3.142 * (rp * rp) * ZF * SIGV
						NTN = NTN + NTNDOT * D
					else if((500.lt.TeV).and.(TeV.lt.(1D+03))) then
						!4626
						SIGV = 2 * (1D-28) * (TeV / 1000.0)**7.7
						!4660
						NTNDOT = 0.5 * NI * NI * 3.142 * (rp * rp) * ZF * SIGV
						NTN = NTN + NTNDOT * D
					else if(((1D+03).lt.TeV).and.(TeV.lt.(1D+04))) then
						!4625
						SIGV = 2 * (1D-28) * (TeV / 1000.0)**3.63
						!4660
						NTNDOT = 0.5 * NI * NI * 3.142 * (rp * rp) * ZF * SIGV
						NTN = NTN + NTNDOT * D
					else if(TeV.gt.(1D+04)) then
						!4624
						SIGV = 2.4 * (1D-26) * (TeV / 1000.0)**1.55
						!4660
						NTNDOT = 0.5 * NI * NI * 3.142 * (rp * rp) * ZF * SIGV
						NTN = NTN + NTNDOT * D
					end if
				end if				
			end if
			! Calculate rate of net power emission, absorption-corrected
			! 4700
			PRAD = (PBR + PLN + PREC)* AB

			if (sflag1.eq.1) then 
				if(sflag2.eq.1) then			!4720
					PRADS = sfactor * PRADS		!4740
					PRAD = PRADS
					If ((-PRAD).gt.(-PBB)) then	!4745
						PRAD = PBB
					end if
				end if
				if(sflag3.eq.1) then			!4730
					PRADS = sfactor * PRADS		!4740
					PRAD = PRADS
					if((-PRAD).gt.(-PBB)) then	!4745
						PRAD = PBB
					end if
				end if
				if (AB.gt.(1 / 2.7183)) then
					!GoTo 4750
				else
					sfactor = PRAD / PRADS
					sflag3 = 1
					PRADS = sfactor * PRADS		!4740
					PRAD = PRADS
					if ((-PRAD).gt.(-PBB)) then	!4745
						PRAD = PBB
					end if
				end if
			else
				sflag1 = 1
				if (AB.gt.(1 / 2.7183)) then
					!GoTo 4750
				else
					sfactor = 1
					sflag2 = 1
					PRADS = sfactor * PRADS		!4740
					PRAD = PRADS
					if((-PRAD).gt.(-PBB)) then	!4745
						PRAD = PBB
					end if
				end if
			end if
			! 4750
			QDOT = PJ + PRAD	
	 		!Compute slow piston speed
			E2 = (1 / g) * (rp / I) * IDOT
			E3 = (1 / (g + 1)) * (rp / ZF) * VZ
			! E5 term in VP (related to dQ/dt) not corrected.

			correctfactor = dissociatenumber * (1 + z) * NI * bc
			correctfactor = 1
			E5 = (4 * Pi * (g - 1) / (MU * g * ZF)) * ((rp * correctfactor) / (I * I * CURRF * CURRF)) * QDOT
			E4 = (g - 1) / g
			VP = (-E2 - E3 + E5) / E4

			IDOT1 = V0 - CH / C0
			IDOT2 = -MUP * (DLog(RADB / rp)) * VZ * I * CURRF
			IDOT3 = MUP * I * ZF * VP * CURRF / rp
			! 4900
			IDOT4 = -I * (R * CURRF + R0)
			! 4920
			IDOT5 = L0 + MUP * (DLog(C)) * z0 * CURRF + MUP * (DLog(RADB / rp)) * ZF * CURRF

			IDOT = (IDOT1 + IDOT2 + IDOT3 + IDOT4) / IDOT5
			ZFDOT = (((MU * (g + 1)) / (16 * Pi * Pi * RHO)) **(0.5)) * I * CURRF / rp
			! 4970
			V = MUP * I * ((DLog(RADB / rp)) * VZ - (ZF / rp) * VP) &
				+ MUP * (((DLog(RADB / rp)) * ZF) + (DLog(C)) * z0) * IDOT + R * I

			! 4980
			V = V * CURRF
			T = T + D
			I = I + IDOT * D
			CH = CH + I * D
			rpOLD = rp
			rp = rp + VP * D
			SPR = -VP * 1D-04
			!Set Variable time increment to suit both slow and fast piston
			if (SPR.lt.1D+02) then
				D = DREAL / 5
			end if
			if (SPR.eq.1D+02) then
				D = DREAL / 10
			end if
			if (SPR.gt.1D+02) then
				D = DREAL / 100
			end if
			if (SPR.gt.1D+03) then
				D = DREAL / 1D+04
			end if
			if (SPR.gt.1D+04) then
				D = DREAL / 1D+05
			end if
			if (SPR.gt.1D+05) then
				D = DREAL / 1D+06
			end if
			if (SPR.gt.1D+06) then
				D = DREAL / 1D+07
			end if
			if (SPR.gt.1D+07) then
				D = DREAL / 1D+08
			end if
			if (SPR.gt.1D+08) then
				D = DREAL / 1D+09
			end if
			!Rem Set limit for piston position
			! 5047
			if(rp.lt.(0.0000001 * RADA)) then 
				ENDFLAG = 1
			end if
			If( ENDFLAG.eq.1) Then
				!GoTo 7000
			else

				! 5050
				ZF = ZF + ZFDOT * D
				QJ = QJ + PJ * D
				QBR = QBR + PBR * D
				! 5052
				QLN = QLN + PLN * D
				
				! 5053
				if (TM.lt.0.86*1D+06) then
					QREC = QREC + PREC * D
				end if
				! uncorrected for absorption
				QSXR = QLN1 + PLN * D
				
				! Corrected for absorption


				Q = Q + QDOT * D
				! 5056
				QRAD = (QRAD + PRAD * D)

				! estimate proportion of each radiation component using their unabsorbed values:
				! hence estimate absorption corrected QBR, QLN, QREC
				QTOTAL = (QBR + QLN + QREC)
				if( -QTOTAL.lt.0.000001) then
					! GoTo 5070
				else
					rbr = QBR / QTOTAL
					rln = QLN / QTOTAL
					rrec = QREC / QTOTAL
					rsxr = QSXR / QTOTAL
					QBR = rbr * QRAD
					QLN = rln * QRAD
					QREC = rrec * QRAD
				end if
				!5070
				TR = T * 1D+06
				VR = V * 1D-03
				IR = I * 1D-03
				KPR = rp * 1D+03
				ZFR = ZF * 1D+03
				SPR = VP * 1D-04
				SZR = ZFDOT * 1D-04
				IDOTKAUS = IDOT * 1D-09
				TMB = TM
				TRRadial = TR * 1D+03 - trradialstart * (1D+09)

				if (IR.gt.Ipeak) then
					Ipeak = IR
				end if
				! Determine max induced voltage for beam-gas neutron yield computation
				if(VR.gt.VRMAX) then
					VRMAX = VR
				end if
				if (KPR.lt.amin) then
					amin = KPR
				end if
				if(ZFR.gt.zmax) then
					zmax = ZFR
				end if
				if(TM.gt.Tpinch) then
					Tpinch = TM
				end if
				If (NI.gt.nipinch) then
					nipinch = NI
				end if
				if (TM.lt.Tpinchmin) then
					tpinchmin = TM
				end if
				! D is time increment in secs, DKPR is piston position increment & DZFR length position increment, both in SI units
				!  Integrate to find EINP, energy dissipated by dynamic resistance effect, 
				! which is 0.5 (Ldot) I^2, considering current taking part in the motion

				EINP = EINP + (1D-07)*(SZR*(1D+04)*DLog(1000*RADB/KPR) - (SPR*(1D+04)*(1000/KPR)*(ZFR/1000))) &
					* IR * IR * (1D+06) * CURRF * CURRF * D
				write(3,'(8(ES14.7X))')TR,IR,VR,ZFR/10+z0*100,rp,rp,SZR,SPR
				write(4,'(8(ES14.7X))')TR,QJ,QBR,PM,QLN,QREC,PRADS,NTN
				write(10,'(2(ES14.7X))')TR,EINP
				write(11,'(2(ES14.7X))')TR,TM
				! Set limit for duration of radiative phase using transit time of small disturbance across pinch radius
				TRAD = TRAD1
				If (T.gt.(TStart + TRAD)) then
					ENDFLAG = 2
				end if
				!5350	GoTo 4100
			end if
		end do
		!7000 
		write(2,'(/"(!) Slow compression phase stopped either on preset time or on RP limit.")')
		write(2,'(/"- final phase parameters:")')
		write(2,'("QLN="ES14.7X)')QLN
		write(2,'("Final plasma Temp. =",ES14.7," [eV]; (",ES14.7," [K] )")')TeV, TeV*11604.45
		! Slow compression Phase Stopped: Time limit or RP limit.
		write(2,'("ENDFLAG="I3)') ENDFLAG
		TSlowcompressionphase = (T - TStart) * 1E+09
		write(2,'("(*) PHASE 4 ends at ="ES14.7X"[us]")')TR
		write(5,'(I1A1ES14.7)')4,tab,TR
		write(2,'("TSlowcompressionphase ="ES14.7X)')TSlowcompressionphase
		write(2,'("T radial="ES14.7X)') (TRRadial * 1E-03)
		write(2,'("T radial + T star="ES14.7X)') (TRRadial * 1E-03) + trradialstart * (1E+06)
		! Calculate energy in inductances
		Ecap = 0.5 * C0 * V0 * V0
		EL0 = 0.5 * L0 * I * I
		ELt = 0.5 * (MU / (2 * Pi)) * (DLog(C)) * z0 * I * I * CURRF * CURRF
		ELp = 0.5 * (MU / (2 * Pi)) * (DLog(RADB / rp)) * ZF * I * I * CURRF * CURRF
		MAG = MU * I * CURRF / (2 * Pi * rp)
		EMAGp = (MAG * MAG / (2 * MU)) * Pi * rp * rp * ZF
		EL0 = (EL0 / Ecap) * 100
		ELt = (ELt / Ecap) * 100
		ELp = (ELp / Ecap) * 100
		EMAGp = (EMAGp / Ecap) * 100
		SFI0 = (I0 * 1E-03) / (RADA * 100) / (DSqrT(P0))
		SFIpeak = (Ipeak) / (RADA * 100) / (DSqrT(P0))
		SFIdip = (IR) / (RADA * 100) / (DSqrT(P0))

		Kmin = KPR / (RADA * 1000) !7010
		ID = (Ipeak) / (RADA * 100)
		! 7020
		Ec = 0.5 * C0 * (V0 - (CH / C0)) * (V0 - (CH / C0))
		Ec = (Ec / Ecap) * 100
		Excircuit = 100 - (EL0 + ELt + ELp + EMAGp + Ec)
		EINP = (EINP / Ecap) * 100
		! calculate loss of energy from inductances during current dip; ignoring capacitative energy change
		L0t = L0 + (MU / (2 * Pi)) * (DLog(C)) * z0
		Lp = (MU / (2 * Pi)) * (DLog(RADB / rp)) * ZF
		Edip = 0.5 * L0t * CURRF * CURRF * (Ipeak * Ipeak * 1E+06 - I * I) - 0.5 * Lp * CURRF * CURRF * I * I
		Edip = (Edip / Ecap) * 100

		SFIpinch = (Ipinch) / (RADA * 100) / (DSqrT(P0))

		write(2,'("Ipeak="ES14.7X)')Ipeak
		write(2,'("Ipinch="ES14.7X)')Ipinch
		write(2,'("SFIpeak="ES14.7X)')SFIpeak
		write(2,'("SFIpinch="ES14.7X)')SFIpinch
		write(2,'("VRMAX="ES14.7X)')VRMAX
		write(2,'("Kmin="ES14.7X)')Kmin
		write(2,'("Ecap="ES14.7X)')Ecap / 1000
		write(2,'("EINP="ES14.7X)')EINP

		! Calculate  neutron yield in D2; 2 components viz 1. thermonuclear, 2. Beam-gas :


		! 7050 
		If (MW.eq.2) then

		else
			! Computed VRMAX varies typically in range of 30-60kV for small to big devices;
			! too low compared to expt observations;
			! Multiplying by factor 2 will get the range closer to 50-100kV;
			! the range generally observed to be reponsible for beam-target neutrons in PF
			! Multiply VRMAX by factor to get closer to experimental observations; 
			! fine tuned to fit the optimum pressure for Yn for the UNU/ICTP PFF (around 3-3.5 torr at 15 kV)
			
			!VRMAX=42.1
			VRMAX = VRMAX * 3
			if(MW.eq.5) then
				!7100
				sig = (409 + ((((1.076 - 0.01368 * (VRMAX))**2) + 1)**(-1)) * 50200) &
					/ (VRMAX * (DExp(45.95 * VRMAX**(-0.5)) - 1))
				sig = sig * 1D-28

			else
				! For deuterium, compute 2. Beam-gas neutron yield component (ref: NRL Formulary 2006 pg 43)
				 sig = ((((1.177 - (3.08 * (1D-04)) * (VRMAX))**2) + 1)**(-1)) * 482 &
					/ (VRMAX * (DExp(47.88 * VRMAX**(-0.5)) - 1))
				 sig = sig * 1D-28
			end if
			! Calibrate for UNU/ICTP PFF for max neutron yield at optimum pressure as 10^8
			! 7150
			sig = sig * 6.34 * 1D+08
			! Change Calibration to NESSI-like, at expt point of 0.5MA pinch current
			sig = sig / 23.23
			! correct for dissociation number change in June 2016
			! DNchange
			!sig = sig / dissociatenumber
			! Use model Ni I^2 zf^2 LOG(b/rp) VRMAX^-0.5 sig

			NBN = NI * ((Ipinch * 1000)**2) * (ZF**2) * (DLog(RADB / rp)) * (VRMAX**(-0.5)) * sig
			NN = NBN + NTN

		end if
! 7300
		write(2,'("--------------"/"Y_th="ES14.7X" [neutrons] ")')NTN
		write(2,'("Y_bt="ES14.7X" [neutrons] "/"--------------")')NBN


		write(2,'("Ecap="ES14.7X)')Ecap / 1000
		write(2,'("RESF="ES14.7X)')RESF
		write(2,'("C="ES14.7X)')C
		write(2,'("L0="ES14.7X)')L0 * 1E+09
		write(2,'("C0="ES14.7X)')C0 * 1E+06
		write(2,'("R0="ES14.7X)')R0 * 1E+03

		write(2,'("RADB="ES14.7X)')RADB * 100
		write(2,'("RADA="ES14.7X)')RADA * 100

		write(2,'("z0="ES14.7X)')z0 * 100
		write(2,'("V0="ES14.7X)')V0 / 1000
		write(2,'("P0="ES14.7X)')P0
		write(2,'("Spitzer R="ES14.7X "[Ohm]")')R


		write(2,'("ID="ES14.7X " [KA/cm] ")')ID / 1E+06
	



		!ActiveSheet.Cells(17, 22) = VRMAX (before x3 for effective beam energy)

		!write(2,'("Yn="ES14.7X" [neutrons] ")')NN
		write(2,'("Yline="ES14.7X" [J] ")')QLN
		write(2,'("Yline="ES14.7X" % ")')(-100 * QLN) / Ecap
		write(2,'("massf="ES14.7X)')Massf
		write(2,'("currf="ES14.7X)')CURRF
		write(2,'("massfr="ES14.7X)')massfr
		write(2,'("currfr="ES14.7X)')currfr


		write(2,'("trradialstart="ES14.7X" [] ")')trradialstart * (1E+06)

		write(2,'(/"*** Final relevant parameters:")')
		write(2,'("Ipeak="ES14.7X" [KA] ")')Ipeak
		write(2,'("Ipinch start ="ES14.7X" [KA] ")')Ipinch
		write(2,'("peak v_axial="ES14.7X" [cm/us] ")')ZZRpeak
		write(2,'("SFIpeak="ES14.7X)')SFIpeak
		write(2,'("peak radial vs="ES14.7X" [cm/us] ")')peakvs
		write(2,'("peak radial vp="ES14.7X" [cm/us] ")')peakvp
		write(2,'("amin="ES14.7X" [cm] ")')amin / 10
		write(2,'("zmax="ES14.7X" [cm] ")')zmax / 10
		write(2,'("TSlowcompressionphase="ES14.7X" [ns] ")')TSlowcompressionphase
		write(2,'("VRMAX="ES14.7X)')VRMAX/3.0
		write(2,'("Tpinch max ="ES14.7X" [K]")')Tpinch
		write(2,'("Y ="ES14.7X" [neutrons]")')NN
		write(2,'("nipinch="ES14.7X" [1/m3] ")')nipinch
		write(2,'("EINP="ES14.7X" % ")')EINP
		write(2,'("*******************************"/)')

! phase 5 : expandeed axial column phase 
		write(*,*) "starting phase 5..."
		write(2,'(/"*** Phase 5 (expanded axial) parameters ***"/)')
		! Expanded axial phase starts; integrated in normalised quantities
		
		!7380 
		D = 0.005
		BE = BE / CFR

		!7385
		T = T / T0
		I = I / I0
		Io = CH / (I0 * T0)
		ZS = ZF / z0
		ZZ = ZG
		z = 1 + ZS
		L = (DLog(C) + 0.25) / DLog(C)
		H = C * C / (C * C - 1)
		H = DSqrT(H)
		L1 = (DLog(C) + 0.5) / DLog(C)

		!7480
		do while (T.le.2.5)
			T = T + D

			tc5 = 1
			AC = (AL * AL * I * I * L - H * H * ZZ * ZZ) / (1 + H * H * (z - 1))
			II = (1 - Io - BE * I * ZZ * L1 - RESF * I) / (1 + BE * tc5 + BE * L1 * (z - 1))
			ZZ = ZZ + AC * D

			z = z + ZZ * D
			I = I + II * D
			Io = Io + I * D

			M = (1 + (1 / (2 * DLog(C)))) * (z - 1)
			V = BE * ((1 * tc5 + M) * II + I * ZZ * (1 + (1 / (2 * Log(C)))))

			TR = T * T0 * 1E+06
			VR = V * V0 * 1E-03
			IR = I * I0 * 1E-03
			ZZR = (ZZCHAR / AL) * ZZ * 1E-04
			ZR = z * z0 * 100

			write(3,'(9(ES14.7X))')TR,IR,VR,ZR,rp,rp,ZZR,0.0,0.0
			!write(4,'(8(ES14.7X))')TR,0.0,0.0,0.0,0.0,0.0,0.0,0.0
			write(10,'(2(ES14.7X))')TR,0.0
			write(11,'(2(ES14.7X))')TR,0.0
			!Set limit for integration to just over half cycle
		end do
		write(2,'("Final phase time ="ES14.7X" [us] ")')TR
		! "INTEGRATION COMPLETED UP TO DESIRED TIME"
		!9888
		write(*,*) "* All phases completed successfully. Stop. *"
		write(2,'(/"*** All phases completed successfully. ***"/"Stop."/)')
	end if


	close(2)
	close(3)
	close(4)
	close(5)
	close(10)
	close(11)
	end program
