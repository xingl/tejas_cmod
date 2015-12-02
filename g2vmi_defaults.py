def set_defaults(self):
	print 'Using settings for Valanju'
	dic = {	# --- VMEC input Parameters ---
		'use_OMP':False,
		'OMP_NUM_THREADS':8,
		'LFORBAL':'F',
		'NSTEP':300,
		'NS_ARRAY':[97, 257, 513],
		'FTOL_ARRAY':[1.0e-06, 1.0e-08, 1.0e-10],
		'NITER_ARRAY':[5000, 30000, 80000],
		'use_PRECON': False,
		'PRECON_TYPE':'GMRES',
		'PREC2D_THRESHOLD':1e-13,
		# --- Grid Parameters ---
		'LASYM':'T',
		'NFP':1,
		'NTOR':4,
		'NTHETA':80,
		'NZETA':48,
		'APHI':[3.191756, -4.669766, 3.471685, -0.993675],
		'use_APHI':False,
		# --- Boundary Parameters ---
		'LFREEB':'F',
		'NVACSKIP': 3,
		'MPOL':100,
		'useDESCUR':True,
		#'mgrid': None  # filled in below
		# --- Current Parameters ---		 		
		'useIcoil':False,
		'useCcoil':False,
		'useEF':False,
		# --- Profile Parameters ---
		'NCURR':1,
		'akima_ac_flag':'density'}
		
	if dic['useEF']:
		dic['mgrid'] = 'None'
	else:
		dic['mgrid'] = 'None'
	return dic



### code snippets from g2vmi.py -> fit_profiles ###
# current density profile: I = integral(jtor * dV/dpsi * dpsi/ds * ds) = integral(js * ds)  -->  js
# dpsipol/dpsitor = 1/q  ->  dpsi/ds = 1/q * psitor_norm/psipol_norm = 1/q * psitor(Sep.)/2pi(psipol(Sep.) - psipol(Axis))
# not used any more, instead js comes as derivative of Itor
#dpsids = abs(self.psitordic['psitor1D'][-1])/(2*pi*(self.ep.siBry - self.ep.siAxis)) / self.ep.PROFdict['q_prof'] 
#js = jtor * Vprime * dpsids
#f_js = scinter.UnivariateSpline(self.ep.PSIdict['psiN1D'], js, s = 0)

#js = f_js(psi)
#if(self.dic['akima_ac_flag'] == 'density'):
#	ac_f = js								# current density: js = jtor * dV/dpsi * dpsi/ds
#else:
#	ac_f = append(0, cumtrapz(js, s))		# integrated current: Itor = integral(js * ds)

