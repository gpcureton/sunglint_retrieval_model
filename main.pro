;------------------------------------------------------------------------------
; The purpose of this program is to plot the modelled glint moments and 
; moment functions, using simulated slope cumulants and cumulant functions as 
; input. This is so that the modelled quantities can be compared to the
; simulated quantities calculated in cumulantFunctionSimulate.pro
;
; Inputs: fileName of simulation results
; Outputs: HDF4 file of modelled results
;
; Author: Geoff Cureton, 29th April, 2009
;------------------------------------------------------------------------------

PRO main,simFileName,modFileName,retFileName,testTau;windowIndex

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Read in input HDF4 file containing simulation results, and initialise   ;;;
	;;; various data structures.                                                ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	simFileName = STRING(simFileName)
	simFileName = STRCOMPRESS(simFileName,/REMOVE_ALL)
	modFileName = STRING(modFileName)
	modFileName = STRCOMPRESS(modFileName,/REMOVE_ALL)
	retFileName = STRING(retFileName)
	retFileName = STRCOMPRESS(retFileName,/REMOVE_ALL)

	simFileID = HDF_SD_START(simFileName, /READ)
	modFileID = HDF_SD_START(modFileName, /READ)

	PRINT,"Input Simulation File:  ",simFileName
	PRINT,"Input Modelling File: ",modFileName 

	PRINT, "Simulation fileID = ",simFileID
	PRINT, "Modelling fileID = ",modFileID

	;;;;;;;;;

	;;; Read the simulated quantities

	PRINT, 'Reading the geometry values...'

	sds_index = hdf_sd_nametoindex(simFileID,"Solar Zenith Angles")
	sds_id = hdf_sd_select(simFileID,sds_index)
	hdf_sd_getdata,sds_id,Solar_Zenith_Angles
	HDF_SD_ENDACCESS, sds_id

	sds_index = hdf_sd_nametoindex(simFileID,"Detector Zenith Angles")
	sds_id = hdf_sd_select(simFileID,sds_index)
	hdf_sd_getdata,sds_id,Detector_Zenith_Angles
	HDF_SD_ENDACCESS, sds_id

	sds_index = hdf_sd_nametoindex(simFileID,"Specular Slopes")
	sds_id = hdf_sd_select(simFileID,sds_index)
	hdf_sd_getdata,sds_id,Specular_Slopes
	HDF_SD_ENDACCESS, sds_id

	sds_index = hdf_sd_nametoindex(simFileID,"Min Slopes")
	sds_id = hdf_sd_select(simFileID,sds_index)
	hdf_sd_getdata,sds_id,Min_Slopes
	HDF_SD_ENDACCESS, sds_id

	sds_index = hdf_sd_nametoindex(simFileID,"Max Slopes")
	sds_id = hdf_sd_select(simFileID,sds_index)
	hdf_sd_getdata,sds_id,Max_Slopes
	HDF_SD_ENDACCESS, sds_id

	;;;;;;;;;

	PRINT, 'Reading some scale values...'

	sds_index = hdf_sd_nametoindex(simFileID,"Cumulant Function 1D length scale")
	sds_id = hdf_sd_select(simFileID,sds_index)
	hdf_sd_getdata,sds_id,lag_1D
	HDF_SD_ENDACCESS, sds_id

	sds_index = hdf_sd_nametoindex(simFileID,"Cumulant Function 2D length scale")
	sds_id = hdf_sd_select(simFileID,sds_index)
	hdf_sd_getdata,sds_id,lag_2D
	HDF_SD_ENDACCESS, sds_id

	;;;;;;;;;

	PRINT, 'Reading the simulated slope and glint moments and cumulants...'

	sds_index = hdf_sd_nametoindex(simFileID,"Slope Moments")
	sds_id = hdf_sd_select(simFileID,sds_index)
	hdf_sd_getdata,sds_id,slopeMoments
	dim_id = hdf_sd_dimgetid(sds_id,0)
	hdf_sd_dimget, dim_id, COUNT=N_moments 
	HDF_SD_ENDACCESS, sds_id

	sds_index = hdf_sd_nametoindex(simFileID,"Glint Moments")
	sds_id = hdf_sd_select(simFileID,sds_index)
	hdf_sd_getdata,sds_id,glintFirstMoments
	HDF_SD_ENDACCESS, sds_id

	sds_index = hdf_sd_nametoindex(simFileID,"Slope Cumulants")
	sds_id = hdf_sd_select(simFileID,sds_index)
	hdf_sd_getdata,sds_id,slopeCumulants
	HDF_SD_ENDACCESS, sds_id

	sds_index = hdf_sd_nametoindex(simFileID,"Glint Cumulants")
	sds_id = hdf_sd_select(simFileID,sds_index)
	hdf_sd_getdata,sds_id,glintCumulants
	HDF_SD_ENDACCESS, sds_id

	;;;;;;;;;

	PRINT, 'Reading the simulated slope and glint second moment functions...'

	sds_index = hdf_sd_nametoindex(simFileID,"Average Slope Second Moment Function")
	sds_id = hdf_sd_select(simFileID,sds_index)
	hdf_sd_getdata,sds_id,slopeSecondMomentFunction
	dim_id = hdf_sd_dimgetid(sds_id,0)
	hdf_sd_dimget, dim_id, COUNT=N 
	HDF_SD_ENDACCESS, sds_id
																		  
	sds_index = hdf_sd_nametoindex(simFileID,"Average Glint Second Moment Function")
	sds_id = hdf_sd_select(simFileID,sds_index)
	hdf_sd_getdata,sds_id,glintSecondMomentFunction
	HDF_SD_ENDACCESS, sds_id

	;;;;;;;;;

	PRINT, 'Reading the simulated slope and glint second cumulant functions...'

	sds_index = hdf_sd_nametoindex(simFileID,"Average Slope Second Cumulant Function")
	sds_id = hdf_sd_select(simFileID,sds_index)
	hdf_sd_getdata,sds_id,slopeSecondCumulantFunction
	HDF_SD_ENDACCESS, sds_id
																		  
	sds_index = hdf_sd_nametoindex(simFileID,"Average Glint Second Cumulant Function")
	sds_id = hdf_sd_select(simFileID,sds_index)
	hdf_sd_getdata,sds_id,glintSecondCumulantFunction
	HDF_SD_ENDACCESS, sds_id

	;;;;;;;;;

	PRINT, 'Reading the simulated slope and glint third moment functions...'

	sds_index = hdf_sd_nametoindex(simFileID,"Average Slope Third Moment Function")
	sds_id = hdf_sd_select(simFileID,sds_index)
	hdf_sd_getdata,sds_id,slopeThirdMomentFunction
	dim_id = hdf_sd_dimgetid(sds_id,0)
	hdf_sd_dimget, dim_id, COUNT=NN 
	HDF_SD_ENDACCESS, sds_id
																		  
	sds_index = hdf_sd_nametoindex(simFileID,"Average Glint Third Moment Function")
	sds_id = hdf_sd_select(simFileID,sds_index)
	hdf_sd_getdata,sds_id,glintThirdMomentFunction
	dim_id = hdf_sd_dimgetid(sds_id,2)
	hdf_sd_dimget, dim_id, COUNT=N_angles
	HDF_SD_ENDACCESS, sds_id

	;;;;;;;;;

	PRINT, 'Reading the slope and glint third cumulant functions...'
																		  
	sds_index = hdf_sd_nametoindex(simFileID,"Average Slope Third Cumulant Function")
	sds_id = hdf_sd_select(simFileID,sds_index)
	hdf_sd_getdata,sds_id,slopeThirdCumulantFunction
	HDF_SD_ENDACCESS, sds_id
																		  
	sds_index = hdf_sd_nametoindex(simFileID,"Average Glint Third Cumulant Function")
	sds_id = hdf_sd_select(simFileID,sds_index)
	hdf_sd_getdata,sds_id,glintThirdCumulantFunction
	HDF_SD_ENDACCESS, sds_id

	;;;;;;;;;

	PRINT, "Closing HDF access to file ",simFileName
	HDF_SD_END, simFileID

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;;; Read the modelled quantities

	sds_index = hdf_sd_nametoindex(modFileID,"Modelled Glint Third Moment Function")
	sds_id = hdf_sd_select(modFileID,sds_index)
	hdf_sd_getdata,sds_id,M3_L_modelled
	HELP,M3_L_modelled
	HDF_SD_ENDACCESS, sds_id

	PRINT, "Closing HDF access to file ",modFileName
	HDF_SD_END, modFileID

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	d2r = !DPI/180.D
	r2d = 180.D/!DPI

	PRINT,""
	PRINT,"N = ",N
	PRINT,"NN = ",NN
	PRINT,"N_angles = ",N_angles
	PRINT,"N_moments = ",N_moments
	PRINT,'delta_x (1D) = ',lag_1D[1],' meters', Format='(A,F20.10,A)'
	PRINT,'delta_x (2D) = ',lag_2D[1],' meters', Format='(A,F20.10,A)'
	delta_x = lag_1D[1]
	PRINT,'Solar Zenith Angles (radians) = ',Solar_Zenith_Angles
	PRINT,'Solar Zenith Angles (degrees) = ',Solar_Zenith_Angles*r2d
	PRINT,'Detector Zenith Angles (radians) = ',Detector_Zenith_Angles
	PRINT,'Detector Zenith Angles (degrees) = ',Detector_Zenith_Angles*r2d
	PRINT,'Minimum Slopes = ',Min_Slopes
	PRINT,'Specular Slopes = ',Specular_Slopes
	PRINT,'Maximum Slopes = ',Max_Slopes

	PRINT,""
	PRINT,"Simulated Slope Moments(N_moments) =      ",slopeMoments
	PRINT,"Simulated Glint First Moments(N_angles) = "
	PRINT,glintFirstMoments
	PRINT,""
	PRINT,"Simulated Slope Cumulants(N_angles) =     ",slopeCumulants
	PRINT,"Simulated Glint Cumulants(N_angles,N_moments) = "
	PRINT,glintCumulants
	PRINT,"Simulated Glint Cumulants(N_moments,N_angles) = "
	PRINT,TRANSPOSE(glintCumulants)
	PRINT,""

	x_max = FLOAT(N-1L)*delta_x ; meters
	k_max = 2.D*!DPI/delta_x
	delta_k = k_max/DOUBLE(N-1L)
	k_N = DOUBLE(N/2L)*delta_k

	PRINT,'x_max =   ',x_max,' meters', Format='(A,F20.10,A)'
	PRINT,'k_max =   ',k_max,' meters^{-1}', Format='(A,F20.10,A)'
	PRINT,'delta_k = ',delta_k,' meters^{-1}', Format='(A,F20.10,A)'
	PRINT,'Nyquist Wavenumber = ',k_N,' meters^{-1}',Format='(A,F20.12,A)'

	;;; Turn on error reporting
	!EXCEPT=2

	xwinsize=800
	ywinsize=450

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Initialise the COMMON blocks, for passing data to other functions       ;;;
	;;; and subroutines which cannot be easily passed as parameters             ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;;; Common block for the integration limits
	COMMON GEOMCOMMON, xiMin, xiMax

	;;; Common Block to pass slope cumulants to procedures and functions
	;;; which cannot passed as parameters to those procedures and functions.
	COMMON SLOPECUMULANTS,kappa_xi_1,kappa_xi_2,kappa_xi_3,kappa_xi_11,kappa_xi_21,kappa_xi_12, $
		kappa_xi_110,kappa_xi_101,kappa_xi_011, $
		kappa_xi_210,kappa_xi_120, $
		kappa_xi_201,kappa_xi_102, $
		kappa_xi_012,kappa_xi_021, $
		kappa_xi_111

	;;; Common Block containing the glint moments
	COMMON GLINTMOMENTS,mu_L_1,mu_L_2,mu_L_3,mu_L_11,mu_L_21,mu_L_12,mu_L_111

	;;; Common Block to pass lambda variables to procedures and functions
	;;; which cannot passed as parameters to those procedures and functions.
	COMMON LAMBDA,lambda2_xi_Matrix,lambda3_xi_Matrix,lambda_xi_Det,lambda2_xi_invMatrix,lambda3_xi_invMatrix

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Initialise the various structures, which are used to pass data          ;;;
	;;; in an encapsulated fashion.                                             ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;;; Initialise the GEOM structure, which will be passed to the glint moment routine(s)

	GEOM = {N_angles:0L,source_angle:DINDGEN(N_angles),xi_min:DINDGEN(N_angles),$
		xi_0:DINDGEN(N_angles), xi_max:DINDGEN(N_angles)}

	beta = 0.68D*d2r

	GEOM.N_angles = N_angles
	GEOM.source_angle = Solar_Zenith_Angles
	gamma = (GEOM.source_angle - Detector_Zenith_Angles)/2.D
	GEOM.xi_0 = Specular_Slopes
	GEOM.xi_min = Min_Slopes
	GEOM.xi_max = Max_Slopes
	dxi = GEOM.xi_max - GEOM.xi_min

	xiMin = GEOM.xi_min
	xiMax = GEOM.xi_max

	;;; Initialise the SCALE structure, which will contain the various values of 
	;;; spatial and wavenumber increments

	SCALE = {N:123L,delta_x:0.D,x_max:0.D,k_max:0.D,delta_k:0.D,NN:64L}

	x = DBLARR(N/2L+1L)
	;x = DINDGEN(N/2L+1L)*delta_x
	x = lag_1D[0L:N/2L]

	SCALE.N = N
	SCALE.NN = NN
	SCALE.delta_x = delta_x
	SCALE.x_max = x_max
	SCALE.k_max = k_max
	SCALE.delta_k = delta_k

	;;; Initialise the SLOPE structure, which contains the slope moments and 
	;;; second and third slope moment functions
	SLOPE = {moments:DBLARR(N_moments),secondMomentFunction:DBLARR(N/2L+1L),thirdMomentFunction:DBLARR(NN/2L+1L,NN/2L+1L)}

	SLOPE.moments = slopeMoments
	SLOPE.secondMomentFunction = slopeSecondMomentFunction[0:N/2L]
	SLOPE.thirdMomentFunction = slopeThirdMomentFunction[0:NN/2L,0:NN/2L]

	;;; Initialise the GLINT structure, which contains the first glint moment, 
	;;; second glint moment function and third glint moment function for each
	;;; geometry defined in GEOM.

	thirdMomentFunc = {thirdMomentFunction:DBLARR(NN/2L+1L,NN/2L+1L)}
	thirdMomentFunctionGeom = REPLICATE(thirdMomentFunc, N_angles)

	GLINT = {geomIndex:0L,firstMoment:DBLARR(N_angles),secondMomentFunction:DBLARR(N_angles,N/2L+1L),$
		thirdMomentFunctionGeom:thirdMomentFunctionGeom}
	GLINT.firstMoment = glintFirstMoments
	GLINT.secondMomentFunction = TRANSPOSE(glintSecondMomentFunction[0:N/2L,*])
	;HELP,glintSecondMomentFunction
	;HELP,GLINT.secondMomentFunction
	;HELP,glintSecondCumulantFunction
	;HELP,glintThirdMomentFunction
	;HELP,glintThirdCumulantFunction

	FOR geometry=0L,N_angles-1L DO BEGIN
		GLINT.thirdMomentFunctionGeom[geometry].thirdMomentFunction = $
			glintThirdMomentFunction[0:NN/2L,0:NN/2L,geometry]
		;HELP,GLINT.thirdMomentFunctionGeom[geometry].thirdMomentFunction
	ENDFOR

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Initialise the simulated slope moments and moment functions, and the    ;;;
	;;; simulated slope cumulants and cumulant functions from quantities read   ;;;
	;;; from the input file.                                                    ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	PRINT,"Initialising simulated slope variables..."

	;;; The simulated slope moments and cumulants
	mu_xi_sim = DBLARR(N_moments)
	mu_xi_sim = slopeMoments
	kappa_xi_sim = DBLARR(N_moments)
	kappa_xi_sim = slopeCumulants

	;;; The second order slope moment and cumulant functions
	M2_xi_sim = DBLARR(N/2L+1L)
	C2_xi_sim = DBLARR(N/2L+1L)
	M2_xi_sim = SLOPE.secondMomentFunction
	C2_xi_sim = slopeSecondCumulantFunction[0L:N/2L]

	;;; The main slope third moment and cumulant functions
	M3_xi_sim = DBLARR(NN/2L+1L,NN/2L+1L)
	C3_xi_sim = DBLARR(NN/2L+1L,NN/2L+1L)
	M3_xi_sim = SLOPE.thirdMomentFunction
	C3_xi_sim = slopeThirdCumulantFunction[0:NN/2L,0:NN/2L]

	;;; The special slope third moment and cumulant functions, along
	;;; the tau_1 and the (tau_1 = tau_2) axes
	M21_xi_sim = DBLARR(NN/2L+1L)
	C21_xi_sim = DBLARR(NN/2L+1L)
	M12_xi_sim = DBLARR(NN/2L+1L)
	C12_xi_sim = DBLARR(NN/2L+1L)

	FOR j=0L,NN/2L DO BEGIN
		;;; The third slope moment and cumulant functions along the tau_1 axis
		M21_xi_sim[j] = M3_xi_sim[0L,j]
		C21_xi_sim[j] = C3_xi_sim[0L,j]

		;;; The third slope moment and cumulant functions along the (tau_1 = tau_2) axis
		M12_xi_sim[j] = M3_xi_sim[j,j]
		C12_xi_sim[j] = C3_xi_sim[j,j]
	ENDFOR


	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Initialise the simulated glint moments and moment functions, and the    ;;;
	;;; simulated glint cumulants and cumulant functions from quantities read   ;;;
	;;; from the input file.                                                    ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	PRINT,"Initialising simulated glint variables..."

	;;; The simulated first glint moments and cumulants
	mu_L1_sim = DBLARR(N_angles)
	mu_L1_sim = GLINT.firstMoment
	kappa_L_sim = DBLARR(N_angles,N_moments)
	kappa_L_sim = TRANSPOSE(glintCumulants)

	;;; The second order glint moment and cumulant functions
	M2_L_sim = DBLARR(N_angles,N/2L+1L)
	C2_L_sim = DBLARR(N_angles,N/2L+1L)
	M2_L_sim = GLINT.secondMomentFunction
	C2_L_sim = TRANSPOSE(glintSecondCumulantFunction[0L:N/2L,*])

	;;; The main glint third moment and cumulant functions
	M3_L_sim = DBLARR(N_angles,NN/2L+1L,NN/2L+1L)
	C3_L_sim = DBLARR(N_angles,NN/2L+1L,NN/2L+1L)
	FOR geometry=0L,N_angles-1L DO BEGIN
		M3_L_sim[geometry,*,*] = GLINT.thirdMomentFunctionGeom[geometry].thirdMomentFunction
	ENDFOR
	C3_L_sim = TRANSPOSE(glintThirdCumulantFunction[0L:NN/2L,0L:NN/2L,*])
	
	;;; The special glint third moment and cumulant functions, along
	;;; the tau_1 and the (tau_1 = tau_2) axes
	M21_L_sim = DBLARR(N_angles,NN/2L+1L)
	C21_L_sim = DBLARR(N_angles,NN/2L+1L)
	M12_L_sim = DBLARR(N_angles,NN/2L+1L)
	C12_L_sim = DBLARR(N_angles,NN/2L+1L)

	FOR geometry=0L,N_angles-1L DO BEGIN
		FOR j=0L,NN/2L DO BEGIN
			;;; The third glint moment function along the tau_1 axis
			M21_L_sim[geometry,j] = M3_L_sim[geometry,0L,j]
			C21_L_sim[geometry,j] = C3_L_sim[geometry,0L,j]

			;;; The third glint moment function along the (tau_1 = tau_2) axis
			M12_L_sim[geometry,j] = M3_L_sim[geometry,j,j]
			C12_L_sim[geometry,j] = C3_L_sim[geometry,j,j]
		ENDFOR
	ENDFOR
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Initialise the modelled glint moments and moment functions, and the     ;;;
	;;; modelled glint cumulants and cumulant functions from quantities read    ;;;
	;;; from the input file.                                                    ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	PRINT,"Initialising modelled glint variables..."

	;;; The special glint third moment and cumulant functions, along
	;;; the tau_1 and the (tau_1 = tau_2) axes
	M21_L_modelled = DBLARR(N_angles,NN/2L+1L)
	M12_L_modelled = DBLARR(N_angles,NN/2L+1L)

	FOR geometry=0L,N_angles-1L DO BEGIN
		FOR j=0L,NN/2L DO BEGIN
			;;; The third glint moment function along the tau_1 axis
			M21_L_modelled[geometry,j] = M3_L_modelled[0L,j,geometry]

			;;; The third glint moment function along the (tau_1 = tau_2) axis
			M12_L_modelled[geometry,j] = M3_L_modelled[j,j,geometry]
		ENDFOR
	ENDFOR

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;		Retrieve the slope cumulants from the        ;;;
	;;;     the simulated glint moments.                  ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	PRINT,'#####################################################',Format='(/A/)'
	PRINT,'Retrieving the slope cumulants...',Format='(A/)'

	fitted = [1,1,1]
	;kappa_xi_guess = [0.D , 0.2D^2.D , 0.D]
	kappa_xi_guess = [0.D , 0.005D , 0.D]

	;;; Retrieve the slope cumulants from the simulated glint moments
	mu_L1_fitted = DBLARR(N_angles)
	slopeCumulant_retrieve,GEOM,mu_L1_sim,kappa_xi_guess,$
		fitted,sigma,chisq,mu_L1_fitted

	;;; Save the retrieved slope cumulants
	kappa_xi_retrieved = kappa_xi_guess

	;;; Calculate the fitted glint cumulants from the fitted glint
	;;; moments for each geometry
	kappa_L_fitted = DBLARR(N_angles,N_moments)
	glintCumulantsFromMoments,mu_L1_fitted,kappa_L_fitted

	PRINT,"Simulated Slope Cumulants(N_moments) =     ",kappa_xi_sim,FORMAT='(A,5E)'
	PRINT,"Retrieved Slope Cumulants(N_moments) =     ",kappa_xi_retrieved,FORMAT='(A,5E)'
	PRINT,"Simulated Glint Moments(N_angles) = ",mu_L1_sim,FORMAT='(A,5E)'
	PRINT,"Fitted Glint Moments(N_angles)    = ",mu_L1_fitted,FORMAT='(A,5E)'

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Model the glint second moment functions for the various geometries,
	;;; using the simulated slope cumulants, and the modelled glint moments associated
	;;; with those cumulants, for comparison with the simulated glint second moment 
	;;; functions for the various geometries
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	PRINT,'#####################################################',Format='(/A/)'
	PRINT,'Retrieving the slope second cumulant function...',Format='(A/)'

	;;; Define the lag variable using the basic spatial increment
	tau = DBLARR(N/2L+1L)
	tau = DINDGEN(N/2L+1L)*delta_x

	;;; Define the lag values at which the glint second moment function will be modelled
	;;; (not including the case of zero lag)
	N_tau_samples = 2L*30L   ;;; Must be an even number!!! We have this many samples 
	                         ;;; in each scale region.
	tauLo = tau[1]
	tauMid = 10.D
	tauHi = tau[N/2]

	tau_index = LONARR(2L*N_tau_samples+1L)
	tau_reduced = DBLARR(2L*N_tau_samples+1L)
	tau_index[0L] = 0L
	tau_reduced[0L] = 0.D

	createTauAbcissa,N,N_tau_samples,delta_x,tauLo,tauMid,tauHi,tau,tau_index,tau_reduced

	M2_L_temp1 = DBLARR(N_angles)
	M2_L_temp2 = DBLARR(N_angles)
	M2_L_fitted = DBLARR(N_angles,N_tau_samples+1L)
	C2_L_fitted = DBLARR(N_angles,N_tau_samples+1L)
	C2_xi_retrieved = DBLARR(N_tau_samples+1L)

	;;; Copy the retrieved slope cumulants to the common block
	kappa_xi_1 = kappa_xi_retrieved[0]
	kappa_xi_2 = kappa_xi_retrieved[1]
	kappa_xi_3 = kappa_xi_retrieved[2]

	;;; Fix the value of kappa_xi_21 and kappa_xi_12 
	;;; and initialise the guesses
	fitted = [1,0,0]
	kappa_xi_guess = [0.95D*kappa_xi_2,0.D,0.D]
	
	;;; Initialise the GLINTMOMENTS COMMON block with the simulated glint moments,
	;;; so they can be used in the retrieval of the second slope moment and cumulant functions
	mu_L_1 = mu_L1_sim

	;;; For zero lag
	C2_xi_retrieved[0] = 1.D
	M2_L_fitted[*,0] = 1.D

	FOR tau_sample=1L,N_tau_samples DO BEGIN

		M2_L_temp1 = M2_L_sim[*,tau_index[tau_sample]]
		
		slopeCumulantFunction_retrieve2,GEOM,M2_L_temp1,kappa_xi_guess,$
			fitted,sigma,chisq,M2_L_temp2

		M2_L_fitted[*,tau_sample] = M2_L_temp2

		C2_xi_retrieved[tau_sample] = kappa_xi_guess[0]/kappa_xi_retrieved[1]

		PRINT,'tau_index:',tau_index[tau_sample],$
			'   tau_reduced:',tau_reduced[tau_sample],$
			'   kappa_xi_2_ret:  ',kappa_xi_retrieved[1],$
			'   kappa_11_xi_ret:  ',kappa_xi_guess[0],$
			'   C2_xi_retrieved:  ',C2_xi_retrieved[tau_sample]
	ENDFOR

	;;; Interpolate the retrieved second slope cumulant function to the delta_x spatial increment
	N_interp = N/2+1
	PRINT,"N_interp = ",N_interp
	C2_xi_retrieved_interp = DBLARR(N_interp)
	HELP,C2_xi_retrieved_interp
	C2_xi_retrieved_interp = INTERPOL(C2_xi_retrieved,tau_reduced,tau[0L:N_interp-1L],/SPLINE)

	C2_xi_retrieved = DBLARR(N)
	C2_xi_retrieved[0L:N/2L] = C2_xi_retrieved_interp[0L:N/2L]
	C2_xi_retrieved[N/2L+1L:N-1L] = REVERSE(C2_xi_retrieved_interp[1L:N/2L-1L])
	C2_xi_retrieved[0L] = 1L

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Use a cubic spline to smooth things out a bit

	; Calculate interpolating cubic spline: 
	;C2_xi_temp = SPL_INIT(tau, C2_xi_retrieved) 
	;result = SPL_INTERP(tau, C2_xi_retrieved[0L:N/2L], C2_xi_temp, tau)
	;C2_xi_retrieved[0L:N/2L] = C2_xi_temp

	;C2_xi_temp1 = DBLARR(100)
	;C2_xi_temp2 = DBLARR(100)
	;tau_temp = DBLARR(100)
	;HELP,tau_temp
	;HELP,C2_xi_temp1
	;HELP,C2_xi_temp2

	;tau_temp[0:99] = tau[0L:99]
	;C2_xi_temp1[0:99] = C2_xi_retrieved[0L:99]

	;SPLINE_P, tau_temp, C2_xi_temp1, tau_temp, C2_xi_temp2,TAN0=[1,-30],TAN1=[0.4,-0.05]
	;C2_xi_temp1 = C2_xi_retrieved
	;C2_xi_temp2 = SMOOTH(C2_xi_retrieved,3)
	;C2_xi_retrieved = C2_xi_temp2

	;windowIndex=-1
	;windowIndex++
	;WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Spline interpolation',RETAIN=2

	;PLOT,C2_xi_temp1,xrange=[0,50],linestyle=0,/ystyle,charsize=chsize
	;OPLOT,C2_xi_temp2,linestyle=2
	;OPLOT,C2_xi_sim,linestyle=3

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	PRINT,"C2_xi_retrieved[0]       = ",C2_xi_retrieved[0]
	PRINT,"C2_xi_retrieved[1]       = ",C2_xi_retrieved[1]
	PRINT,"C2_xi_retrieved[N/2L-1L] = ",C2_xi_retrieved[N/2L-1L]
	PRINT,"C2_xi_retrieved[N/2L]    = ",C2_xi_retrieved[N/2L]
	PRINT,"C2_xi_retrieved[N/2L+1L] = ",C2_xi_retrieved[N/2L+1L]
	PRINT,"C2_xi_retrieved[N-1L]    = ",C2_xi_retrieved[N-1L]

	PRINT,"C2_xi_sim[0]       = ",C2_xi_sim[0]
	PRINT,"C2_xi_sim[1]       = ",C2_xi_sim[1]
	PRINT,"C2_xi_sim[N/2L-1L] = ",C2_xi_sim[N/2L-1L]
	PRINT,"C2_xi_sim[N/2L]    = ",C2_xi_sim[N/2L]

	;;; Interpolate the fitted second glint moment function to the delta_x spatial increment
	N_interp = N/2+1
	PRINT,"N_interp = ",N_interp
	M2_L_fitted_interp = DBLARR(N_angles,N_interp)
	HELP,M2_L_fitted_interp

	FOR geometry=0L,N_angles-1L DO BEGIN
		M2_L_fitted_interp[geometry,*] = INTERPOL(M2_L_fitted[geometry,*],tau_reduced,tau[0L:N_interp-1L],/SPLINE)
	ENDFOR
	
	M2_L_fitted = DBLARR(N_angles,N)
	M2_L_fitted[*,0L:N/2L] = M2_L_fitted_interp[*,0L:N/2L]
	M2_L_fitted[*,N/2L+1L:N-1L] = REVERSE(M2_L_fitted_interp[*,1L:N/2L-1L],2)
	M2_L_fitted[*,0L] = 1L

	PRINT,"M2_L_fitted[*,0]       = ",M2_L_fitted[*,0]
	PRINT,"M2_L_fitted[*,1]       = ",M2_L_fitted[*,1]
	PRINT,"M2_L_fitted[*,N/2L-1L] = ",M2_L_fitted[*,N/2L-1L]
	PRINT,"M2_L_fitted[*,N/2L]    = ",M2_L_fitted[*,N/2L]
	PRINT,"M2_L_fitted[*,N/2L+1L] = ",M2_L_fitted[*,N/2L+1L]
	PRINT,"M2_L_fitted[*,N-1L]    = ",M2_L_fitted[*,N-1L]
	PRINT,""

	;;; Compute the modelled linear second glint cumulant function using the simulated slope cumulants
	;;; and the modelled glint cumulants

	C2_L_fitted = DBLARR(N_angles,N)
	FOR geometry=0L,N_angles-1L DO BEGIN
		C2_L_fitted[geometry,*] = $
			mu_L1_fitted[geometry]*(M2_L_fitted[geometry,*]-mu_L1_fitted[geometry])/$
				kappa_L_fitted[geometry,1]
	ENDFOR
	C2_L_fitted[*,0L] = 1L

	PRINT,"C2_L_fitted[*,0]       = ",C2_L_fitted[*,0]
	PRINT,"C2_L_fitted[*,1]       = ",C2_L_fitted[*,1]
	PRINT,"C2_L_fitted[*,N/2L-1L] = ",C2_L_fitted[*,N/2L-1L]
	PRINT,"C2_L_fitted[*,N/2L]    = ",C2_L_fitted[*,N/2L]
	PRINT,"C2_L_fitted[*,N/2L+1L] = ",C2_L_fitted[*,N/2L+1L]
	PRINT,"C2_L_fitted[*,N-1L]    = ",C2_L_fitted[*,N-1L]

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;		Plot the simulated and retrieved second slope cumulant functions	;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	SET_PLOT,"X"
	!P.FONT = -1
	!P.MULTI = 0
	chsize = 1.5
	windowIndex=-1

 	;;;;;;;;
 
 	abcissaTex="\tau (m)" 	;;; tau_{1} (m)
 	ytitleTex="C^{(11)}_{L}(\tau)"		;;; (tau_{1},0)
 
	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex++
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Retrieved second slope cumulant function',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)

 	PLOT,C2_xi_retrieved,xtitle=abcissaStr,xrange=[0.,20],$
 		ytitle=ytitleStr,yrange=[-0.0,1.0],linestyle=0,/ystyle,charsize=chsize
	OPLOT,C2_xi_sim,linestyle=2
	
 	;;;;;;;;
 
 	abcissaTex="\tau (m)" 	;;; tau_{1} (m)
 	ytitleTex="C^{(11)}_{L}(\tau)"		;;; (tau_{1},0)
 
	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex++
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='(SIM-RET) second slope cumulant function',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)

 	PLOT,C2_xi_sim-C2_xi_retrieved,xtitle=abcissaStr,xrange=[0.,20],$
 		ytitle=ytitleStr,yrange=[-0.1,0.1],linestyle=0,/ystyle,charsize=chsize
	OPLOT,C2_xi_sim-C2_xi_retrieved,psym=2
	
 	;;;;;;;;
 
 	abcissaTex="\tau (m)" 	;;; tau_{1} (m)
 	ytitleTex="M^{(11)}_{L}(\tau)"		;;; (tau_{1},0)
 
	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex++
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Fitted second glint moment function',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)

 	PLOT,M2_L_fitted[0,*],xtitle=abcissaStr,xrange=[0.,N-1],$
 		ytitle=ytitleStr,yrange=[-0.01,0.03],linestyle=0,/ystyle,charsize=chsize
	FOR geometry=0L,N_angles-1L DO BEGIN
		OPLOT,M2_L_fitted[geometry,*],linestyle=geometry
	ENDFOR

 	;;;;;;;;
 
 	abcissaTex="\tau (m)" 	;;; tau_{1} (m)
 	ytitleTex="C^{(11)}_{L}(\tau)"		;;; (tau_{1},0)
 
	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex++
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Fitted second glint cumulant function',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)

 	PLOT,C2_L_fitted[0,*],xtitle=abcissaStr,xrange=[0.,N-1],$
 		ytitle=ytitleStr,yrange=[-0.01,0.03],linestyle=0,/ystyle,charsize=chsize
	FOR geometry=0L,N_angles-1L DO BEGIN
		OPLOT,C2_L_fitted[geometry,*],linestyle=geometry
	ENDFOR

	;END

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Model the glint third moment functions along a couple of slices, for the 
	;;; various geometries, using the simulated slope cumulants, and the modelled 
	;;; glint moments associated with those cumulants, for comparison with the 
	;;; simulated glint second moment 
	;;; functions for the various geometries
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	PRINT,'#####################################################',Format='(/A/)'
	PRINT,'Retrieving the third slope cumulant function along the tau_1 axis (21)...',Format='(/A/)'

	M21_L_temp1 = DBLARR(N_angles)
	M21_L_temp2 = DBLARR(N_angles)
	M21_L_fitted = DBLARR(N_angles,NN)
	C21_L_fitted = DBLARR(N_angles,NN)

	C21_xi_retrieved = DBLARR(NN)
	C12_xi_retrieved = DBLARR(NN)

	;;; For zero lag
	C21_xi_retrieved[0] = 1.D
	M21_L_fitted[*,0] = 1.D
	C21_L_fitted[*,0] = 1.D

	;;; Initialise the SLOPECUMULANTS COMMON block with the simulated slope cumulants,
	;;; so they can be used in the modelling of the third glint moment and cumulant functions

	kappa_xi_1 = kappa_xi_retrieved[0]
	kappa_xi_2 = kappa_xi_retrieved[1]
	;kappa_xi_3 = kappa_xi_retrieved[2]
	kappa_xi_3 = kappa_xi_sim[2]

	;;; Set the fitted flags to fix the second slope cumulant function,
	;;; and allow kappa_xi_21 and kappa_xi_12 to vary.
	fitted = [0,1,1]

	;;; Initialise the slope cumulant function guesses
	kappa_xi_guess = [kappa_xi_11,$
		              0.95*kappa_xi_3,$
					 -0.95*kappa_xi_3]

	;testTau = 5L
	;testTau = 7L
	;testTau = 10L
	;testTau = 16L
	;testTau = 32L
	;testTau = 64L
	;testTau = 128L
	;testTau = 256L

	PRINT,""
	PRINT,"kappa_xi_sim       = ",kappa_xi_sim,FORMAT='(A,3E)'
	PRINT,"kappa_xi_retrieved = ",kappa_xi_retrieved,FORMAT='(A,3E)'
	PRINT,""

	;FOR tau1=1L,NN/2L DO BEGIN
	FOR tau1=1L,testTau DO BEGIN
		PRINT,"tau1 = ",tau1

		;;; Pass the second glint moment function at this lag
		;;; to the common block, which can be accessed 
		;;; by  thirdGlintMomentFunction21()
		mu_L_11 = M2_L_sim[*,tau1]

		;;; Use the modelled glint function as the data to be fitted to
		M21_L_temp1 = M21_L_modelled[*,tau1]
		;;; Use the simulated glint function as the data to be fitted to
		;M21_L_temp1 = M21_L_sim[*,tau1]

		;;; Initialise the second slope cumulant function using the 
		;;; retrieved value at this lag, which is fixed in this
		;;; retrieval
		kappa_xi_guess[0] = kappa_xi_2 * C2_xi_retrieved[tau1]

		PRINT,'Simulated slope cumulants (11,21,12):   ',$
			kappa_xi_sim[1]*C2_xi_sim[tau1],$
			kappa_xi_sim[2]*C21_xi_sim[tau1],$
			kappa_xi_sim[2]*C12_xi_sim[tau1],FORMAT='(A,3E)'
		PRINT,'Sim norm slope cumulants (11,21,12) :   ',$
			C2_xi_sim[tau1],$
			C21_xi_sim[tau1],$
			C12_xi_sim[tau1],FORMAT='(A,3E)'
		PRINT,'Initial slope cumulants (11,21,12):     ',kappa_xi_guess,FORMAT='(A,3E)'

		;;; Retrieve kappa_xi_21 and kappa_xi_12 slope cumulant functions
		slopeCumulantFunction_retrieve21,GEOM,M21_L_temp1,kappa_xi_guess,$
			fitted,sigma,chisq,M21_L_temp2

		M21_L_fitted[*,tau1] = M21_L_temp2

		;;; We know that kappa_xi_12 = -kappa_xi_21, and we set 
		;;; appropriately
		kappa_xi_guess[2] = -1.D*kappa_xi_guess[1]

		PRINT,'Retrieved slope cumulants (11,21,12):   ',kappa_xi_guess,FORMAT='(A,3E)'
		
		;;; Compute the normalised slope cumulant functions
		C21_xi_retrieved[tau1] = kappa_xi_guess[1]/kappa_xi_3
		C12_xi_retrieved[tau1] = kappa_xi_guess[2]/kappa_xi_3

		PRINT,'Retr norm slope cumulants (11,21,12):   ',C2_xi_retrieved[tau1],$
			C21_xi_retrieved[tau1],C12_xi_retrieved[tau1],FORMAT='(A,3E)'

		PRINT,"M21_L_mod[geom] = ",M21_L_temp1,FORMAT='(A,5E)'
		PRINT,"M21_L_fit[geom] = ",M21_L_temp2,FORMAT='(A,5E/)'
	ENDFOR

 	;;;;;;;;
 
 	postScriptOut=0L
 	abcissaTex="\tau_{1} (m)" 	;;; tau_{1} (m)
 	ytitleTex="C^{(111)}_{\xi}(\tau_{1},0)"		;;; (tau_{1},0)
 
	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex++
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Retrieved Third slope cumulant function (21) - tau_1',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)
 
 	PLOT,C21_xi_sim,xtitle=abcissaStr, $
		xrange=[1L,NN/2],$
 		ytitle=ytitleStr, $
		yrange=[1.5*MIN(C21_xi_sim),1.5*MAX(C21_xi_sim)], $
 		linestyle=0,/ystyle,/xstyle,charsize=chsize
 	OPLOT,C12_xi_sim,linestyle=3
 	OPLOT,C21_xi_retrieved,psym=4
 	OPLOT,C12_xi_retrieved,psym=5
 
 	;;;;;;;;

 	postScriptOut=0L
 	abcissaTex="\tau_{1} (m)" 	;;; tau_{1} (m)
 	ytitleTex="M^{(111)}_{L}(\tau_{1},0)"		;;; (tau_{1},0)

	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex++
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Modelled Third glint moment function (21) - tau_1',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)
 
 	PLOT,M21_L_modelled[0,*],xtitle=abcissaStr, $
		xrange=[1L,NN/2],$
 		ytitle=ytitleStr, $
		yrange=[0.,0.02], $
		;yrange=[-500.,500.], $
 		linestyle=0,/ystyle,/xstyle,charsize=chsize
 
	FOR geometry=1L,N_angles-1L DO BEGIN
		OPLOT,M21_L_modelled[geometry,*],linestyle=geometry
	ENDFOR

 	;;;;;;;;

	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex++
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Fitted Third glint moment function (21) - tau_1',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)
 
 	PLOT,M21_L_fitted[0,*],xtitle=abcissaStr, $
		xrange=[1L,NN/2],$
 		ytitle=ytitleStr, $
		yrange=[0.,0.02], $
		;yrange=[-500.,500.], $
 		linestyle=0,/ystyle,/xstyle,charsize=chsize
 
	FOR geometry=1L,N_angles-1L DO BEGIN
		OPLOT,M21_L_fitted[geometry,*],linestyle=geometry
	ENDFOR

 	;;;;;;;;
 
	PRINT,'#####################################################',Format='(/A/)'
	PRINT,'Retrieving the third slope cumulant function along the tau_1=tau_2 axis (12)...',Format='(/A/)'

	M12_L_temp1 = DBLARR(N_angles)
	M12_L_temp2 = DBLARR(N_angles)
	M12_L_fitted = DBLARR(N_angles,NN)
	C12_L_fitted = DBLARR(N_angles,NN)

	;;; For zero lag
	C12_xi_retrieved[0] = 1.D
	M12_L_fitted[*,0] = 1.D
	C12_L_fitted[*,0] = 1.D

	;;; Initialise the SLOPECUMULANTS COMMON block with the simulated slope cumulants,
	;;; so they can be used in the modelling of the third glint moment and cumulant functions

	kappa_xi_1 = kappa_xi_retrieved[0]
	kappa_xi_2 = kappa_xi_retrieved[1]
	;kappa_xi_3 = kappa_xi_retrieved[2]
	kappa_xi_3 = kappa_xi_sim[2]

	;;; Set the fitted flags to fix the second slope cumulant function,
	;;; and allow kappa_xi_21 and kappa_xi_12 to vary.
	fitted = [0,0,1]

	;;; Initialise the slope cumulant function guesses
	kappa_xi_guess = [kappa_xi_11,$
		              kappa_xi_3 * C21_xi_retrieved[1],$
					  kappa_xi_3 * C12_xi_retrieved[1]]

	PRINT,""
	PRINT,"kappa_xi_sim       = ",kappa_xi_sim,FORMAT='(A,3E)'
	PRINT,"kappa_xi_retrieved = ",kappa_xi_retrieved,FORMAT='(A,3E)'
	PRINT,""

	;FOR tau=1L,NN/2L DO BEGIN
	FOR tau=1L,testTau DO BEGIN
		PRINT,"tau = ",tau

		;;; Pass the second glint moment function at this lag
		;;; to the common block, which can be accessed 
		;;; by  thirdGlintMomentFunction21()
		mu_L_11 = M2_L_sim[*,tau]

		;;; Use the modelled glint function as the data to be fitted to
		M12_L_temp1 = M12_L_modelled[*,tau]
		;;; Use the simulated glint function as the data to be fitted to
		;M12_L_temp1 = M12_L_sim[*,tau]

		;;; Initialise the second slope cumulant function using the 
		;;; retrieved value at this lag, which is fixed in this
		;;; retrieval.
		kappa_xi_guess[0] = kappa_xi_2 * C2_xi_retrieved[tau]

		;;; Initialise the kappa_xi_21 slope cumulant function using the 
		;;; retrieved value at this lag, which is fixed in this
		;;; retrieval.
		kappa_xi_guess[1] = kappa_xi_3 * C21_xi_retrieved[tau]

		PRINT,'Simulated slope cumulants (11,21,12):   ',$
			kappa_xi_sim[1]*C2_xi_sim[tau],$
			kappa_xi_sim[2]*C21_xi_sim[tau],$
			kappa_xi_sim[2]*C12_xi_sim[tau],FORMAT='(A,3E)'
		PRINT,'Sim norm slope cumulants (11,21,12) :   ',$
			C2_xi_sim[tau],$
			C21_xi_sim[tau],$
			C12_xi_sim[tau],FORMAT='(A,3E)'
		PRINT,'Initial slope cumulants (11,21,12):     ',kappa_xi_guess,FORMAT='(A,3E)'

		;;; Retrieve kappa_xi_12 slope cumulant function
		slopeCumulantFunction_retrieve12,GEOM,M12_L_temp1,kappa_xi_guess,$
			fitted,sigma,chisq,M12_L_temp2

		M12_L_fitted[*,tau] = M12_L_temp2

		;;; We know that kappa_xi_21 = -kappa_xi_12, and we set 
		;;; appropriately
		kappa_xi_guess[1] = -1.D*kappa_xi_guess[2]

		PRINT,'Retrieved slope cumulants (11,21,12):   ',kappa_xi_guess,FORMAT='(A,3E)'
		
		;;; Compute the normalised slope cumulant functions
		C21_xi_retrieved[tau] = kappa_xi_guess[1]/kappa_xi_3
		C12_xi_retrieved[tau] = kappa_xi_guess[2]/kappa_xi_3

		PRINT,'Retr norm slope cumulants (11,21,12):   ',C2_xi_retrieved[tau],$
			C21_xi_retrieved[tau],C12_xi_retrieved[tau],FORMAT='(A,3E)'

		PRINT,"M12_L_mod[geom] = ",M12_L_temp1,FORMAT='(A,5E)'
		PRINT,"M12_L_fit[geom] = ",M12_L_temp2,FORMAT='(A,5E/)'
	ENDFOR

 	;;;;;;;;
 
 	postScriptOut=0L
 	abcissaTex="\tau (m)" 	;;; tau (m)
 	ytitleTex="C^{(111)}_{\xi}(\tau,\tau)"		;;; (tau,tau)
 
	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex++
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Retrieved Third slope cumulant function (12) - tau',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)
 
 	PLOT,C21_xi_sim,xtitle=abcissaStr, $
		xrange=[1L,NN/2],$
 		ytitle=ytitleStr, $
		yrange=[1.5*MIN(C21_xi_sim),1.5*MAX(C21_xi_sim)], $
 		linestyle=0,/ystyle,/xstyle,charsize=chsize
 	OPLOT,C12_xi_sim,linestyle=3
 	OPLOT,C21_xi_retrieved,psym=4
 	OPLOT,C12_xi_retrieved,psym=5
 
 	;;;;;;;;

 	postScriptOut=0L
 	abcissaTex="\tau_{1} (m)" 	;;; tau_{1} (m)
 	ytitleTex="M^{(111)}_{L}(\tau_{1},0)"		;;; (tau_{1},0)

	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex++
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Modelled Third glint moment function (12) - tau_1',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)
 
 	PLOT,M12_L_modelled[0,*],xtitle=abcissaStr, $
		xrange=[1L,NN/2],$
 		ytitle=ytitleStr, $
		yrange=[0.,0.02], $
		;yrange=[-500.,500.], $
 		linestyle=0,/ystyle,/xstyle,charsize=chsize
 
	FOR geometry=1L,N_angles-1L DO BEGIN
		OPLOT,M12_L_modelled[geometry,*],linestyle=geometry
	ENDFOR

 	;;;;;;;;

	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex++
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Fitted Third glint moment function (12) - tau_1',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)
 
 	PLOT,M12_L_fitted[0,*],xtitle=abcissaStr, $
		xrange=[1L,NN/2],$
 		ytitle=ytitleStr, $
		yrange=[0.,0.02], $
		;yrange=[-500.,500.], $
 		linestyle=0,/ystyle,/xstyle,charsize=chsize
 
	FOR geometry=1L,N_angles-1L DO BEGIN
		OPLOT,M12_L_fitted[geometry,*],linestyle=geometry
	ENDFOR

	;END
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Model the glint third moment functions for the various geometries,
	;;; using the simulated slope cumulant function, and the modelled glint 
	;;; moments associated with those cumulants, for comparison with the 
	;;; simulated glint second moment functions for the various geometries.
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	PRINT,'#####################################################',Format='(/A/)'
	PRINT,'Retrieving the slope third cumulant function...',Format='(/A/)'

	M3_L_temp1 = DBLARR(N_angles)
	M3_L_temp2 = DBLARR(N_angles)
	M3_L_fitted = DBLARR(NN,NN,N_angles)
	C3_L_fitted = DBLARR(NN,NN,N_angles)

	C3_xi_retrieved = DBLARR(NN,NN)

	;;; For zero lag
	M3_L_fitted[0L,0L,*] = 1.D
	C3_L_fitted[0L,0L,*] = 1.D
	C3_xi_retrieved[0L,0L] = 1.D

	;;; Initialise the SLOPECUMULANTS COMMON block with the simulated slope cumulants,
	;;; so they can be used in the modelling of the second glint moment and cumulant functions

	kappa_xi_1 = kappa_xi_retrieved[0]
	kappa_xi_2 = kappa_xi_retrieved[1]
	kappa_xi_3 = kappa_xi_retrieved[2]
	;kappa_xi_1 = kappa_xi_sim[0]
	;kappa_xi_2 = kappa_xi_sim[1]
	;kappa_xi_3 = kappa_xi_sim[2]

	;;; Set the fitted flags to fix kappa_xi_21 and kappa_xi_12,
	;;; and allow kappa_xi_111 to vary.
	;fitted = [0,0,1]

	;;; Initialise the slope cumulant function guesses
	kappa_xi_guess = [0.9*kappa_xi_3]
	;kappa_xi_guess = [0.9]

	;;; Initialise the GLINTMOMENTS COMMON block with the modelled glint moments
	mu_L_1 = mu_L1_sim

	;;; Compute the modelled glint third moment function,
	;;; avoiding the Nyquist line at 
	FOR tau1=1L,testTau DO BEGIN
		;FOR tau2=1L,tau1-1L DO BEGIN
		;IF (C12_xi_retrieved[tau1] < 0.D) THEN kappa_xi_guess=[0.95*C12_xi_retrieved[tau1]]
		;IF (C12_xi_retrieved[tau1] > 0.D) THEN kappa_xi_guess=[0.95*C12_xi_retrieved[tau1]]
		;kappa_xi_guess=[0.95*C12_xi_retrieved[tau1]]
		guessSlope = kappa_xi_3*(C12_xi_retrieved[tau1]-C21_xi_retrieved[tau1])/DOUBLE(tau1)
		FOR tau2=tau1-1L,1L,-1 DO BEGIN

			PRINT,'tau1:',tau1,'   tau2:',tau2
			
			;;; Use the modelled glint function as the data to be fitted to
			M3_L_temp1 = M3_L_modelled[tau1,tau2,*]
			;;; Use the simulated glint function as the data to be fitted to
			;M3_L_temp1 = M3_L_sim[tau1,tau2,*]

			;;; The second order slope cumulants
			kappa_xi_110 = kappa_xi_2*C2_xi_retrieved[tau1]
			kappa_xi_101 = kappa_xi_2*C2_xi_retrieved[tau2]
			kappa_xi_011 = kappa_xi_2*C2_xi_retrieved[ABS(tau2-tau1)]
			;kappa_xi_110 = kappa_xi_2*C2_xi_sim[tau1]
			;kappa_xi_101 = kappa_xi_2*C2_xi_sim[tau2]
			;kappa_xi_011 = kappa_xi_2*C2_xi_sim[ABS(tau2-tau1)]

			;;; The third order slope cumulants. Since we are calculating
			;;; in the region (tau1>tau2), the lag values (tau2-tau1) will be
			;;; -ve. Since kappa_xi_021(-tau) = -kappa_xi_021(tau), we use
			;;; the slope cumulant values at +ve lag, multiplied by -1. 
			kappa_xi_210 =      kappa_xi_3*C21_xi_retrieved[tau1]
			kappa_xi_120 =      kappa_xi_3*C12_xi_retrieved[tau1]
			kappa_xi_201 =      kappa_xi_3*C21_xi_retrieved[tau2]
			kappa_xi_102 =      kappa_xi_3*C12_xi_retrieved[tau2]
			kappa_xi_021 = -1.D*kappa_xi_3*C21_xi_retrieved[ABS(tau2-tau1)]
			kappa_xi_012 = -1.D*kappa_xi_3*C12_xi_retrieved[ABS(tau2-tau1)]
			;kappa_xi_210 =      kappa_xi_3*C21_xi_sim[tau1]
			;kappa_xi_120 =      kappa_xi_3*C12_xi_sim[tau1]
			;kappa_xi_201 =      kappa_xi_3*C21_xi_sim[tau2]
			;kappa_xi_102 =      kappa_xi_3*C12_xi_sim[tau2]
			;kappa_xi_021 = -1.D*kappa_xi_3*C21_xi_sim[ABS(tau2-tau1)]
			;kappa_xi_012 = -1.D*kappa_xi_3*C12_xi_sim[ABS(tau2-tau1)]

			kappa_xi_guess = [kappa_xi_3*C21_xi_retrieved[tau1]] ;; Good?
			;kappa_xi_guess = [guessSlope*DOUBLE(tau2) + kappa_xi_3*C21_xi_retrieved[tau1]] ;; 
			;kappa_xi_guess = [kappa_xi_3*C3_xi_sim[tau1,tau2]]
			;kappa_xi_guess = [C3_xi_sim[tau1,tau2]]
			;kappa_xi_guess = [prevRow[tau1]]

			PRINT,'Simulated slope cumulants      :   ',kappa_xi_1,kappa_xi_2,kappa_xi_3
			PRINT,'Simulated slope cumulants (111):   ',kappa_xi_sim[2]*C3_xi_sim[tau1,tau2],FORMAT='(A,E)'
			PRINT,'Sim norm slope cumulants  (111):   ',C3_xi_sim[tau1,tau2],FORMAT='(A,E)'
			PRINT,'Initial slope cumulants   (111):   ',kappa_xi_guess[0],FORMAT='(A,E)'

			;;; Retrieve kappa_xi_21 and kappa_xi_12 slope cumulant functions
			slopeThirdCumulantFunction_retrieve,GEOM,M3_L_temp1,kappa_xi_guess,$
				sigma,chisq,M3_L_temp2

			M3_L_fitted[tau1,tau2,*] = M3_L_temp2
			M3_L_fitted[tau2,tau1,*] = M3_L_temp2

			PRINT,'Retrieved slope cumulants (111):   ',kappa_xi_guess[0],FORMAT='(A,E)'

			;;; Compute the normalised slope cumulant functions
			C3_xi_retrieved[tau1,tau2] = kappa_xi_guess[0]/kappa_xi_3
			C3_xi_retrieved[tau2,tau1] = kappa_xi_guess[0]/kappa_xi_3

			PRINT,'Ret norm slope cumulants (111):    ',C3_xi_retrieved[tau1,tau2],FORMAT='(A,E)'

			PRINT,"M3_L_mod[geom] = ",TRANSPOSE(M3_L_temp1),FORMAT='(A,5E)'
			PRINT,"M3_L_fit[geom] = ",M3_L_temp2,FORMAT='(A,5E/)'
		ENDFOR
	ENDFOR

	;;; Plug in the values of C3_xi_retrieved along the tau1 and tau1=tau2 axes...
	FOR tau=1L,testTau DO BEGIN
		C3_xi_retrieved[tau,0L]  = C21_xi_retrieved[tau]
		C3_xi_retrieved[0L,tau]  = C21_xi_retrieved[tau]
		C3_xi_retrieved[tau,tau] = C12_xi_retrieved[tau]
	ENDFOR

	;;; Fill the rest of the retrieved slope third cumulant function array
	biCovarianceSymmetry,C3_xi_retrieved,NN

	;;; Plug in the values of M3_L_fitted along the tau1 and tau1=tau2 axes...
	FOR tau1=1L,testTau DO BEGIN
		M3_L_fitted[tau1,0L,*] = M21_L_fitted[*,tau1]
		M3_L_fitted[0L,tau1,*] = M21_L_fitted[*,tau1]
		M3_L_fitted[tau1,tau1,*] = M12_L_fitted[*,tau1]
	ENDFOR

	;; Fill the rest of the fitted third moment function array
	temp_3 = DBLARR(NN,NN)
	FOR geometry=0L,N_angles-1L DO BEGIN
		temp_3 = M3_L_fitted[*,*,geometry]
		biCovarianceSymmetry,temp_3,NN
		M3_L_fitted[*,*,geometry] = temp_3
	ENDFOR

	;;; Compute the fitted third glint cumulant function
	FOR geometry=0L,N_angles-1L DO BEGIN
		FOR tau1=0L,testTau DO BEGIN
			FOR tau2=0L,tau1 DO BEGIN
				C3_L_fitted[tau1,tau2,geometry] = (mu_L1_fitted[geometry]*M3_L_fitted[tau1,tau2,geometry] $
					- (mu_L1_fitted[geometry]^2.D)*(M2_L_fitted[geometry,tau1] $
					+ M2_L_fitted[geometry,tau2] $
					+ M2_L_fitted[geometry,ABS(tau2-tau1)]) $
					+ 2.D*(mu_L1_fitted[geometry]^3.D))/kappa_L_fitted[geometry,2L]

				C3_L_fitted[tau2,tau1,geometry] = C3_L_fitted[tau1,tau2,geometry]
			ENDFOR
		ENDFOR
		C3_L_fitted[0L,0L,geometry] = 1.D
	ENDFOR

	;; Fill the rest of the fitted third cumulant function array
	FOR geometry=0L,N_angles-1L DO BEGIN
		temp_3 = C3_L_fitted[*,*,geometry]
		biCovarianceSymmetry,temp_3,NN
		C3_L_fitted[*,*,geometry] = temp_3
	ENDFOR

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;	Glint Third Moment Function                 ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;;--- 3D plot setup
	xwinsize=800
	ywinsize=450
	
	rotate_z=-20.
	rotate_x=10.
	z_axis=2
	ysize=1.*ywinsize
	xsize=ROUND(1.2*ysize)

	;plotLim = 16L
	;plotLim = 32L
	;plotLim = 64L

	plotLim = testTau

	bispecMin = MIN(DOUBLE(glintThirdMomentFunction[*,*,0]))
	bispecMax = MAX(DOUBLE(glintThirdMomentFunction[*,*,0]))
	PRINT,"glint third moment function minimum: Angle "+STRING(0)+" = ",bispecMin
	PRINT,"glint third moment function maximum: Angle "+STRING(0)+" = ",bispecMax

	FOR geometry=0,N_angles-1L DO BEGIN
		!P.MULTI = 0
		windowIndex++
		window,windowIndex,xsize = 2*xsize,ysize = 2*ysize,title='Modelled Glint Third Moment Function: geometry:'+STRING(geometry),RETAIN=2

		shiftFunc = SHIFT(M3_L_modelled[*,*,geometry],NN/2-1L,NN/2-1L)
		
		;SURFACE,DOUBLE(shiftFunc[NN/2-1L-plotLim:NN/2-1L+plotLim,NN/2-1L-plotLim:NN/2-1L+plotLim]),AZ=rotate_z,$
		SURFACE,DOUBLE(shiftFunc),AZ=rotate_z,$
			;zrange=[bispecMin,bispecMax], $
			zrange=[-0.01,0.20], $
			/xstyle,/ystyle,/zstyle,$
			charsize=3.5,xtitle='i',ytitle='j',ZAXIS=z_axis
	ENDFOR
 
	;;;;;;;;;;;;;;;

	FOR geometry=0,N_angles-1L DO BEGIN
		!P.MULTI = 0
		windowIndex++
		window,windowIndex,xsize = 2*xsize,ysize = 2*ysize,title='Fitted Glint Third Moment Function: geometry:'+STRING(geometry),RETAIN=2

		shiftFunc = SHIFT(M3_L_fitted[*,*,geometry],NN/2-1L,NN/2-1L)
		
		;SURFACE,DOUBLE(shiftFunc[NN/2-1L-plotLim:NN/2-1L+plotLim,NN/2-1L-plotLim:NN/2-1L+plotLim]),AZ=rotate_z,$
		SURFACE,DOUBLE(shiftFunc),AZ=rotate_z,$
			;zrange=[bispecMin,bispecMax], $
			zrange=[-0.01,0.20], $
			/xstyle,/ystyle,/zstyle,$
			charsize=3.5,xtitle='i',ytitle='j',ZAXIS=z_axis
	ENDFOR

	;;;;;;;;;;;;;;;


	!P.MULTI = 0
	windowIndex++
	window,windowIndex,xsize = 2*xsize,ysize = 2*ysize,title='Simulated Slope Third Cumulant Function Surface',RETAIN=2

	bispecMin = MIN(DOUBLE(slopeThirdCumulantFunction))
	bispecMax = MAX(DOUBLE(slopeThirdCumulantFunction))
	PRINT,"Simulated slope third moment function minimum:",bispecMin
	PRINT,"Simulated slope third moment function maximum:",bispecMax

	shiftFunc = SHIFT(slopeThirdCumulantFunction,NN/2-1L,NN/2-1L)
	
	;SURFACE,DOUBLE(shiftFunc[NN/2-1L-plotLim:NN/2-1L+plotLim,NN/2-1L-plotLim:NN/2-1L+plotLim]),AZ=rotate_z,$
	SURFACE,DOUBLE(shiftFunc),AZ=rotate_z,$
		zrange=[bispecMin,bispecMax], $
		/xstyle,/ystyle,/zstyle,$
		charsize=3.5,xtitle='i',ytitle='j',ZAXIS=z_axis

	;;;;;;;;;;;;;;;

	!P.MULTI = 0
	windowIndex++
	window,windowIndex,xsize = xsize,ysize = ysize,title='Simulated Slope Third Cumulant Function Contour',RETAIN=2

	contourMin = MIN(DOUBLE(slopeThirdCumulantFunction))
	contourMax = MAX(DOUBLE(slopeThirdCumulantFunction))
	contourRange=contourMax-contourMin
	PRINT,"Simulated slope third moment function minimum:",contourMin
	PRINT,"Simulated slope third moment function maximum:",contourMax
	
	;CONTOUR,DOUBLE(shiftFunc[NN/2-1L-plotLim:NN/2-1L+plotLim,NN/2-1L-plotLim:NN/2-1L+plotLim]), $
	CONTOUR,DOUBLE(shiftFunc), $
		/xstyle,/ystyle,/zstyle, $
		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(50L)/DOUBLE(49L)

	;;;;;;;;;;;;;;;

	!P.MULTI = 0
	windowIndex++
	window,windowIndex,xsize = 2*xsize,ysize = 2*ysize,title='Retrieved Slope Third Cumulant Function Surface',RETAIN=2

	shiftFunc = SHIFT(C3_xi_retrieved,NN/2-1L,NN/2-1L)
	
	;SURFACE,DOUBLE(shiftFunc[NN/2-1L-plotLim:NN/2-1L+plotLim,NN/2-1L-plotLim:NN/2-1L+plotLim]),AZ=rotate_z,$
	SURFACE,DOUBLE(shiftFunc),AZ=rotate_z,$
		zrange=[bispecMin,bispecMax], $
		;zrange=[-0.01,0.04], $
		/xstyle,/ystyle,/zstyle,$
		charsize=3.5,xtitle='i',ytitle='j',ZAXIS=z_axis

	;;;;;;;;;;;;;;;

	!P.MULTI = 0
	windowIndex++
	window,windowIndex,xsize = xsize,ysize = ysize,title='Retrieved Slope Third Cumulant Function Contour',RETAIN=2

	;CONTOUR,DOUBLE(shiftFunc[NN/2-1L-plotLim:NN/2-1L+plotLim,NN/2-1L-plotLim:NN/2-1L+plotLim]), $
	CONTOUR,DOUBLE(shiftFunc), $
		/xstyle,/ystyle,/zstyle, $
		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(50L)/DOUBLE(49L)

	;END

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;   Open the output HDF file   ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	PRINT,"Open output filename: ",retFileName

	retFile = fileinfo(retFileName)
	help, fileinfo(retFileName), /structure
	PRINT, retFile

	IF (NOT retFile.EXIST) THEN BEGIN
		;;; Create and open file using SD interface
		fileID = HDF_SD_START(retFileName, /CREATE)
		;fileID = HDF_OPEN(fileName, /CREATE,/WRITE)
		PRINT, 'Created new HDF file: ',retFileName
	ENDIF ELSE BEGIN
		;;; Create and open file using SD interface
		PRINT, 'HDF file ',retFileName,' exists, opening...'
		fileID = HDF_SD_START(retFileName, /RdWr)
		;fileID = HDF_OPEN(fileName, /WRITE)
		PRINT, 'Opened HDF file ',retFileName,' for reading and writing'
	ENDELSE

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;   Add some attributes to the file   ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	IF (NOT retFile.EXIST) THEN BEGIN
		PRINT,"Writing global attributes to ",retFileName
		HDF_SD_ATTRSET, fileID, 'DATE', SYSTIME()
		HDF_SD_ATTRSET, fileID, 'EXPERIMENT', 'sunglintRetrievalModel'
		HDF_SD_ATTRSET, fileID, 'NAME', 'Geoff Cureton'
		HDF_SD_ATTRSET, fileID, 'EMAIL ADDRESS', 'geoff.cureton@physics.org'
	ENDIF ELSE BEGIN
		HDF_SD_ATTRSET, fileID, 'DATE', SYSTIME()+" Zulu"
	ENDELSE

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;   Add some datasets to the file   ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	   
	IF (NOT retFile.EXIST) THEN BEGIN

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the geometry information to global attributes, and variables   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

		PRINT, 'Writing geometry angles and slopes...'

		sourceAngleID  = HDF_SD_CREATE(fileID, "Solar Zenith Angles", [N_angles], /FLOAT)
		HDF_SD_ADDDATA, sourceAngleID, Solar_Zenith_Angles
		numAnglesDimID = HDF_SD_DIMGETID(sourceAngleID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, sourceAngleID 

		detectorAngleID  = HDF_SD_CREATE(fileID, "Detector Zenith Angles", [N_angles], /FLOAT)
		HDF_SD_ADDDATA, detectorAngleID,Detector_Zenith_Angles
		numAnglesDimID = HDF_SD_DIMGETID(detectorAngleID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, detectorAngleID 

		specularSlopeID  = HDF_SD_CREATE(fileID, "Specular Slopes", [N_angles], /FLOAT)
		HDF_SD_ADDDATA, specularSlopeID,Specular_Slopes
		numAnglesDimID = HDF_SD_DIMGETID(specularSlopeID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, specularSlopeID 

		minSlopeID  = HDF_SD_CREATE(fileID, "Min Slopes", [N_angles], /FLOAT)
		HDF_SD_ADDDATA, minSlopeID,Min_Slopes
		numAnglesDimID = HDF_SD_DIMGETID(minSlopeID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, minSlopeID 

		maxSlopeID  = HDF_SD_CREATE(fileID, "Max Slopes", [N_angles], /FLOAT)
		HDF_SD_ADDDATA, maxSlopeID,Max_Slopes
		numAnglesDimID = HDF_SD_DIMGETID(maxSlopeID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, maxSlopeID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the retrieved slope cumulants            ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		PRINT, 'Writing the retrieved slope cumulants...'

		slopeCumulantID = HDF_SD_CREATE(fileID, 'Retrieved Slope Cumulants',[N_moments],    /FLOAT) 
		HDF_SD_ADDDATA, slopeCumulantID, kappa_xi_retrieved
		numCumulantsDimID = HDF_SD_DIMGETID(slopeCumulantID, 0)
		HDF_SD_DIMSET, numCumulantsDimID, LABEL='Number of Cumulants', NAME='N_moments'
		HDF_SD_ENDACCESS, slopeCumulantID
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the fitted glint moments and cumulants   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		PRINT, 'Writing the Fitted glint moments and cumulants...'
		
		glintMomentID = HDF_SD_CREATE(fileID, 'Fitted Glint First Moments',     [N_angles], /FLOAT)
		HDF_SD_ADDDATA, glintMomentID, mu_L1_fitted
		numAnglesDimID = HDF_SD_DIMGETID(glintMomentID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, glintMomentID
		
		glintCumulantID = HDF_SD_CREATE(fileID, 'Fitted Glint Cumulants',     [N_angles,N_moments], /FLOAT)
		HDF_SD_ADDDATA, glintCumulantID, kappa_L_fitted
		numAnglesDimID = HDF_SD_DIMGETID(glintCumulantID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		numCumulantsDimID = HDF_SD_DIMGETID(glintCumulantID, 1)
		HDF_SD_DIMSET, numCumulantsDimID, LABEL='Number of Cumulants', NAME='N_moments'
		HDF_SD_ENDACCESS, glintCumulantID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Set 1D wavenumber scale   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		wavenumber = FINDGEN(N)*delta_k
		wavenumberID  = HDF_SD_CREATE(fileID, "Power Spectrum wavenumber scale", [N], /FLOAT)
		HDF_SD_ADDDATA, wavenumberID,  wavenumber
		HDF_SD_ATTRSET, wavenumberID,  'units', 'meters^{-1}'
		HDF_SD_ATTRSET, wavenumberID,  'increment', delta_k
		wavenumberDimID = HDF_SD_DIMGETID(wavenumberID, 0)
		HDF_SD_DIMSET, wavenumberDimID, LABEL='Data length', NAME='N', SCALE=wavenumber, UNIT='meters^{-1}'
		HDF_SD_ENDACCESS, wavenumberID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Set 1D lag scale   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		length = FINDGEN(N)*delta_x
		lengthID  = HDF_SD_CREATE(fileID, "Cumulant Function 1D length scale", [N], /FLOAT)
		HDF_SD_ADDDATA, lengthID,  length
		HDF_SD_ATTRSET, lengthID,  'units', 'meters'
		HDF_SD_ATTRSET, lengthID,  'increment', delta_x
		lengthDimID = HDF_SD_DIMGETID(lengthID, 0)
		HDF_SD_DIMSET, lengthDimID, LABEL='Data length', NAME='N', SCALE=length, UNIT='meters'
		HDF_SD_ENDACCESS, lengthID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the Retrieved Slope Second Cumulant Function   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		PRINT, 'Writing retrieved slope second cumulant function ...'
		
		slopeSecondCumulantFunctionID = HDF_SD_CREATE(fileID, "Retrieved Slope Second Cumulant Function"    , [N], /FLOAT)
		HDF_SD_ADDDATA, slopeSecondCumulantFunctionID, C2_xi_retrieved
		HDF_SD_ATTRSET, slopeSecondCumulantFunctionID, 'long_name', $
			'Retrieved Slope Second Cumulant Function, secondCumulantFunction[0 ... N-1]'
		secondCumulantFunctionDimID = HDF_SD_DIMGETID(slopeSecondCumulantFunctionID, 0)
		HDF_SD_DIMSET, secondCumulantFunctionDimID, LABEL='Second Cumulant Function Data length', NAME='N', SCALE=length, UNIT='meters'
		HDF_SD_ENDACCESS, slopeSecondCumulantFunctionID
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the Fitted Glint Second Moment Function   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		PRINT, 'Writing fitted glint second moment function ...'
		
		glintSecondMomentFunctionID = HDF_SD_CREATE(fileID, "Fitted Glint Second Moment Function"    , [N_angles, N], /FLOAT)
		HDF_SD_ADDDATA, glintSecondMomentFunctionID, M2_L_fitted
		HDF_SD_ATTRSET, glintSecondMomentFunctionID, 'long_name', $
			'Fitted glint Second Moment Function, secondMomentFunction[0 ... N_angles-1][0 ... N-1]'
		numAnglesDimID = HDF_SD_DIMGETID(glintSecondMomentFunctionID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		secondMomentFunctionDimID = HDF_SD_DIMGETID(glintSecondMomentFunctionID, 1)
		HDF_SD_DIMSET, secondMomentFunctionDimID, LABEL='Second Moment Function Data length', NAME='N', SCALE=length, UNIT='meters'
		HDF_SD_ENDACCESS, glintSecondMomentFunctionID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the Fitted Glint Second Cumulant Function   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		PRINT, 'Writing fitted glint second cumulant function ...'
		
		glintSecondCumulantFunctionID = HDF_SD_CREATE(fileID, "Fitted Glint Second Cumulant Function"    , [N_angles, N], /FLOAT)
		HDF_SD_ADDDATA, glintSecondCumulantFunctionID, C2_L_fitted
		HDF_SD_ATTRSET, glintSecondCumulantFunctionID, 'long_name', $
			'Fitted glint Second Cumulant Function, secondCumulantFunction[0 ... N_angles-1][0 ... N-1]'
		numAnglesDimID = HDF_SD_DIMGETID(glintSecondCumulantFunctionID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		secondCumulantFunctionDimID = HDF_SD_DIMGETID(glintSecondCumulantFunctionID, 1)
		HDF_SD_DIMSET, secondCumulantFunctionDimID, LABEL='Second Cumulant Function Data length', NAME='N', SCALE=length, UNIT='meters'
		HDF_SD_ENDACCESS, glintSecondCumulantFunctionID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Set 2D wavenumber scale   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		wavenumber2 = FINDGEN(NN)*delta_k
		wavenumber2ID  = HDF_SD_CREATE(fileID, "Bispectrum wavenumber scale", [NN], /FLOAT)
		HDF_SD_ADDDATA, wavenumber2ID,  wavenumber2
		HDF_SD_ATTRSET, wavenumber2ID,  'units', 'meters^{-1}'
		HDF_SD_ATTRSET, wavenumber2ID,  'increment', delta_k
		wavenumber2DimID = HDF_SD_DIMGETID(wavenumber2ID, 0)
		HDF_SD_DIMSET, wavenumber2DimID, LABEL='Data length', NAME='NN', SCALE=wavenumber2, UNIT='meters^{-1}'
		HDF_SD_ENDACCESS, wavenumber2ID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Set 2D lag scale   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		length2 = FINDGEN(NN)*delta_x
		length2ID  = HDF_SD_CREATE(fileID, "Cumulant Function 2D length scale", [NN], /FLOAT)
		HDF_SD_ADDDATA, length2ID,  length2
		HDF_SD_ATTRSET, length2ID,  'units', 'meters'
		HDF_SD_ATTRSET, length2ID,  'increment', delta_x
		length2DimID = HDF_SD_DIMGETID(length2ID, 0)
		HDF_SD_DIMSET, length2DimID, LABEL='Data length', NAME='NN', SCALE=length2, UNIT='meters'
		HDF_SD_ENDACCESS, length2ID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the Retrieved Slope Third Cumulant Function    ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		PRINT, 'Writing retrieved slope and third cumulant function ...'
		
		slopeThirdMomentFunctionID = HDF_SD_CREATE(fileID, "Retrieved Slope Third Cumulant Function"    , [NN,NN], /FLOAT)
		HDF_SD_ADDDATA, slopeThirdMomentFunctionID, C3_xi_retrieved
		HDF_SD_ATTRSET, slopeThirdMomentFunctionID, 'long_name', $
			'Retrieved Slope Third Cumulant Function, ThirdMomentFunction[0 ... NN-1][0 ... NN-1]'
		thirdMomentFunctionDimID = HDF_SD_DIMGETID(slopeThirdMomentFunctionID, 0)
		HDF_SD_DIMSET, thirdMomentFunctionDimID, LABEL='Third Moment Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		thirdMomentFunctionDimID = HDF_SD_DIMGETID(slopeThirdMomentFunctionID, 1)
		HDF_SD_DIMSET, thirdMomentFunctionDimID, LABEL='Third Moment Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		HDF_SD_ENDACCESS, slopeThirdMomentFunctionID
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the Fitted Glint Third Moment Function   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		PRINT, 'Writing fitted glint third moment function ...'
		
		glintThirdMomentFunctionID = HDF_SD_CREATE(fileID, "Fitted Glint Third Moment Function"    , [NN,NN, N_angles], /FLOAT)
		HDF_SD_ADDDATA, glintThirdMomentFunctionID, M3_L_fitted
		HDF_SD_ATTRSET, glintThirdMomentFunctionID, 'long_name', $
			'Fitted Glint Third Moment Function, ThirdMomentFunction[0 ... NN-1][0 ... NN-1][0 ... N_angles-1]'
		thirdMomentFunctionDimID = HDF_SD_DIMGETID(glintThirdMomentFunctionID, 0)
		HDF_SD_DIMSET, thirdMomentFunctionDimID, LABEL='Third Moment Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		thirdMomentFunctionDimID = HDF_SD_DIMGETID(glintThirdMomentFunctionID, 1)
		HDF_SD_DIMSET, thirdMomentFunctionDimID, LABEL='Third Moment Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		numAnglesDimID = HDF_SD_DIMGETID(glintThirdMomentFunctionID, 2)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, glintThirdMomentFunctionID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the Fitted Glint Third Cumulant Function   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		PRINT, 'Writing fitted glint third cumulant function ...'
		
		glintThirdCumulantFunctionID = HDF_SD_CREATE(fileID, "Fitted Glint Third Cumulant Function"    , [NN,NN, N_angles], /FLOAT)
		HDF_SD_ADDDATA, glintThirdCumulantFunctionID, C3_L_fitted
		HDF_SD_ATTRSET, glintThirdCumulantFunctionID, 'long_name', $
			'Fitted Glint Third Cumulant Function, ThirdCumulantFunction[0 ... NN-1][0 ... NN-1][0 ... N_angles-1]'
		thirdCumulantFunctionDimID = HDF_SD_DIMGETID(glintThirdCumulantFunctionID, 0)
		HDF_SD_DIMSET, thirdCumulantFunctionDimID, LABEL='Third Cumulant Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		thirdCumulantFunctionDimID = HDF_SD_DIMGETID(glintThirdCumulantFunctionID, 1)
		HDF_SD_DIMSET, thirdCumulantFunctionDimID, LABEL='Third Cumulant Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		numAnglesDimID = HDF_SD_DIMGETID(glintThirdCumulantFunctionID, 2)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, glintThirdCumulantFunctionID

	ENDIF ELSE BEGIN

	ENDELSE

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;   Close the output HDF file   ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	PRINT, 'Write Operation Completed'
	HDF_SD_END, fileID
	PRINT, '*** File Closed ***'
	PRINT, ''

END
