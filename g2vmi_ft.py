#! /usr/bin/python
"""
g2vmi.py
Author: A. Wingen
Date: Mar. 20. 2014
----------------------------------------------------------------------------------------
Creates a VMEC input file
Usage: for help and how to run from command line, run     
                      python g2vmi.py -h
Initialize:
  import VMEC.Python.g2vmi as g2vmi
  GTI = g2vmi.g2vmi(gfile, fluxLimit, Terms, dic, path, Time, tag)
Input:
  gfile (string):     EFIT geqdsk file name or (full or rel.) pathname
  fluxLimit (float):  percent of poloidal flux surface to use as VMEC outer boundary
  Terms (string):     number of polynomial terms for AC,AI,AM; use '000000' for Akima splines
  dic (dictionary):   optional; user provided input parameters; default: self.set_defaults
  path (string):      optional; path to write VMEC input file to; default is '.' (cwd)
  Time (int):         optional; time for I- & C-coil currents; default is g-file time; 
                      ignored, if diiidsup.in file already exists in path
  tag (string):       optional; arbitrary tag, added to input file name; default is None
Output:
  VMEC input file: 'path/input.' + fitFLAG + '.' + gfile[1::] + '_' + tag
      with fitFLAG = 'akima' for Akima splines, or fitFLAG = 'vmec' for polynomial fit
#
Note:
  calls either DESCUR to optimize the boundary Fourier series or
  Fortran xgtovmi to obtain optimized boundary Fourier series only, otherwise self-sustained.
Usefull hint:
  create a python script called g2vmi_defaults.py somewhere in your PYTHONPATH (outside of VMEC/Python),
  which is the EXACT COPY of the member function set_defaults(self). There you can set and change 
  your own preferences, without changing the source code of this class. The class uses it, 
  if it exists. If not then the built-in set_defaults function is used.
"""

import os
import sys
import copy
import getpass
import scipy.interpolate as scinter
from scipy.integrate import cumtrapz, simps
from  numpy import *
from time import strftime

import EFIT.equilParams_class as EPC
from Misc.deriv import deriv

#-----------------------------------------------------------------------------------------
#--- g2vmi class -------------------------------------------------------------------------
class g2vmi:
	def __init__(self, gfile, fluxLimit, Terms, dic = None, path = '.', Time = None, tag = None):
		# --- input variables ---
		idx = gfile[::-1].find('/')	# returns location of last '/' in gfile or -1 if not found
		if(idx == -1):
			gpath = '.'
		else:
			idx *= -1
			gpath = gfile[0:idx - 1]	# path without a final '/'
			gfile = gfile[idx::]
				
		idx = gfile.find('.')
		fmtstr = '0' + str(idx-1) + 'd'
		shot, time = int(gfile[1:idx]), int(gfile[idx+1::])
		gfile = 'g' + format(shot,fmtstr) + '.' + format(time,'05d')

		fluxLimit = float(fluxLimit)
		Terms = format(int(Terms),'06d')
		
		if(Terms == format(int(000000),'06d')): fitFLAG = 'akima'
		else: fitFLAG = 'vmec'
		
		# --- set member variables ---
		self.gpath = gpath			# path to g-file and a-file
		self.shot = shot			# shot number
		self.TIME = time			# time of g-file
		self.gfile = gfile			# g-file name
		self.fluxLimit = fluxLimit	# percentage of normalized poloidal flux to mark VMEC boundary
		self.Terms = Terms			# number of terms in polynomial fit for AC,AI,AM; 000000 switches to akima splines
		self.fitFLAG = fitFLAG		# flag for poly fit or akima splines
		self.path = path			# path to write VMEC input file		  
		self.filename = 'input.' + fitFLAG + '.' + gfile[1::]	# raw name (no tags) for VMEC input file
				
		self.dic = {}				# dictionary of input file data
		if not (dic == None): self.dic = dic
		else: 
			try:
				import g2vmi_defaults; reload(g2vmi_defaults)
				self.dic = g2vmi_defaults.set_defaults(self) 
			except:
				self.dic = self.set_defaults()
				
		if not (Time == None): self.time = Time	  
		else: self.time = time		# time to get I- & C-coil currents from, if different from gfile time	
				
		if not (tag == None): self.tag = tag
		else: self.tag = None		# tag to append on input file name

		# --- load g-file and get toroidal flux ---
		print 'Loading g-file: ', gpath + '/' + gfile
		self.ep = EPC.equilParams(gpath + '/' + gfile)
		self.psitordic = self.ep.getTorPsi()
	
		# --- read a-file & locate currents ---
		afile = 'a' + format(shot,fmtstr) + '.' + format(time,'05d')
		self.afile = {'name':afile}
		
		if os.path.isfile(gpath + '/' + afile ):
			print 'Loading a-file: ', gpath + '/' + afile 
		
			with open(gpath + '/' + afile, 'r') as f:
				a_input = []
				for line in f:		# read to end of file
					a_input.append(line)
			
			self.afile['data'] = a_input
			
			# locate F-coil & E-coil currents
			for i,line in enumerate(a_input):
				line = line.strip().split()
				if(len(line) == 4):
					line = [float(line[k]) for k in xrange(4)]
					if(line[0] == int(line[0])) & (line[1] == int(line[1])) & (line[2] == int(line[2])) & (line[3] == int(line[3])):
						F_coil_idx = i + int(ceil((line[0] + line[1]) / 4.0)) + 1
						E_coil_idx = i + int(ceil((line[0] + line[1] + line[2]) / 4.0)) + 1
						break
			self.afile['F_coil_idx'] = F_coil_idx
			self.afile['E_coil_idx'] = E_coil_idx
		elif (not self.dic['mgrid']) | (self.dic['mgrid'] == 'None') | (self.dic['mgrid'] == 'none'):
			print 'no a-file found, but no mgrid file either, so set to fixed boundary mode'
			self.dic['mgrid'] = None
			self.dic['LFREEB'] = 'F' 
		else:
			raise RuntimeError('a-file not found')
		

	#-------------------------------------------------------------------------------------
	#--- Member Functions ----------------------------------------------------------------

	# --- set defaults in dic ---
	# sets the self.dic member variable
	def set_defaults(self):
		print 'Using default settings'
		dic = {	# --- Runtime Parameters ---
			'use_OMP':True,
			'OMP_NUM_THREADS':4,
			'LFORBAL':'F',
			'NSTEP':300,
			'NS_ARRAY':[13, 25, 49, 97, 193, 385],
			'FTOL_ARRAY':[1.0e-8, 1.0e-9, 1.0e-9, 1.0e-10, 1.0e-10, 1.0e-16],
			'NITER_ARRAY':[1500,3000,8000,10000],
			'use_PRECON':True,
			'PRECON_TYPE':'GMRES',
			'PREC2D_THRESHOLD':1e-12,
			# --- Grid Parameters ---
			'LASYM':'T',
			'LRFP':'F',	# not used
			'NFP':1,
			'NTOR':4,
			'NTHETA':101,
			'NZETA':48,
			'APHI':[3.191756, -4.669766, 3.471685, -0.993675],
			'use_APHI':True,
			# --- Boundary Parameters ---
			'LFREEB':'F',
			'NVACSKIP': 3,
			'MPOL':25,
			'useDESCUR':True,
			#'mgrid': None  # filled in below
			# --- Current Parameters ---		 		
			'coilDICT':{'smTime':25,
						'iCoilBranch':'blessed',
						'path':self.path,
						'quiet':True},
			'useIcoil':True,
			'useCcoil':True,
			'useEF':False,
			# --- Profile Parameters ---
			'NCURR':1,
			'akima_ac_flag':'current'}
			
		if dic['useEF']:
			dic['mgrid'] = '/home/shared/vmec/makegrid/raw/single_fp_ef/mgrid_d3d_efbic_kp48.nc'
		else:
			dic['mgrid'] = '/home/shared/vmec/makegrid/raw/single_fp_noef/mgrid_d3d_efbic_kp48.nc'
		return dic

	# --- get I-coil & C-coil currents ---
	# Ccoils = ['C79', 'C139', 'C199', 'C259', 'C319', 'C19']
	# Icoils = ['iu30', 'iu90', 'iu150', 'iu210', 'iu270', 'iu330',
	#		    'il30', 'il90', 'il150', 'il210', 'il270', 'il330']
	def get_CIcoils(self):
		print 'Reading I-coil & C-coil currents ...',
		if (not self.dic['useIcoil']) & (not self.dic['useCcoil']):
			print 'Coils not used'
			return zeros(3), zeros(12)
		elif os.path.isfile(self.dic['coilDICT']['path'] + '/diiidsup.in'):	# read diiidsup.in if it already exists
			print 'from diiidsup.in file'
		else:																# else retrieve from MDS+
			print '... from MDS+ database'
			from Misc.Getd3dat import Getd3dat
			Getd3dat(self.shot, self.time, smoothingTime = self.dic['coilDICT']['smTime'], 
					 iCoilDataSet = self.dic['coilDICT']['iCoilBranch'], quiet = self.dic['coilDICT']['quiet'], 
					 path = self.dic['coilDICT']['path'])
		
		_, Ccoils, Icoils = read_diiidsup(path = self.dic['coilDICT']['path'] + '/')
		if not self.dic['useIcoil']: Icoils *= 0
		if not self.dic['useCcoil']: Ccoils *= 0
		
		return Ccoils, Icoils
	
	
	# --- get F-coil currents ---
	# Fcoils = ['F1A', 'F2A', 'F3A', 'F4A', 'F5A', 'F6A', 'F7A', 'F8A', 'F9A',
	#           'F1B', 'F2B', 'F3B', 'F4B', 'F5B', 'F6B', 'F7B', 'F8B', 'F9B']
	def get_Fcoils(self):
		if not self.dic['mgrid']: return zeros(18)
		print 'Reading F-coil currents from a-file'
		Fturns = array([58,58,58,58,58,55,55,58,55,		# number of turns in Fcoils
						58,58,58,58,58,55,55,58,55])
		idx = self.afile['F_coil_idx']
		Fcoils = (self.afile['data'][idx].strip() 
				+ self.afile['data'][idx + 1].strip() 
				+ self.afile['data'][idx + 2].strip()
				+ self.afile['data'][idx + 3].strip() 
				+ self.afile['data'][idx + 4])
		Fcoils = array([float(cur) for cur in split_data(Fcoils)])
		if not Fcoils.shape[0] == 18:
			print 'F-coil currents not found'
			return zeros(18)
		
		# check if Fcoils are already divided by Fturns or not
		if (mean(abs(Fcoils)) > 1e+4): Fcoils /= Fturns
			
		return Fcoils
		
		
	# --- get E-coil currents ---
	def get_Ecoils(self):
		if not self.dic['mgrid']: return 0, 0
		print 'Reading E-coil currents from a-file'
		idx = self.afile['E_coil_idx']
		Ecoils = self.afile['data'][idx].strip() + self.afile['data'][idx + 1]
		Ecoils = array([float(cur) for cur in split_data(Ecoils)])
		eca = mean(Ecoils[[0,2,4]])
		ecb = mean(Ecoils[[1,3,5]])

		return eca, ecb
		
		
	# --- Get B-coil current ---
	# either use B-coil current from review+
	# ->  D3D has 24 B-coils with 6 turns each; Bcoil is current in one turn, so use 6 times the Bcoil current
	# or as is done here:
	# use toroidal magnetic field Bt0 from g-file
	def get_Bcoils(self):
		print 'Getting B-coil current'
		R0 = self.ep.data.get('rcentr')		# Major Radius of torus in [m]
		Bt0 = self.ep.data.get('bcentr')	# toroidal magnetic field at R0 in [T]
		Bcoil = Bt0*R0/2e-7					# total coil current, necessary to generate Bt0
		Bcoil /= 24.0						# D3D (and VMEC) has 24 B-coils
		
		return Bcoil


	# --- get profiles ---
	def fit_profiles(self):
		print 'Fitting profiles ...',
		if (self.fitFLAG == 'akima'):
			print 'akima splines'		
		else:
			print 'polynomial'
		
		# get Flux Surfaces
		self.FluxSurfList = self.ep.get_allFluxSur()
		
		# current profile: I = integral(jtor dV) = integral(jtor * dV/dpsi * dpsi)
		jtor = self.ep.jtor_profile(self.FluxSurfList)*1e+6		# toroidal current density
		Vprime = self.ep.volume2D(self.FluxSurfList)['Vprime']	# Plasma 2D-Volume Derivative: dV/dpsi
		
		# total current, using cumulative Trapezoid Rule
		Itor = append(0, cumtrapz(jtor * Vprime, self.ep.PSIdict['psiN1D']))		
		f_Itor = scinter.UnivariateSpline(self.ep.PSIdict['psiN1D'], Itor, s = 0)
		
		# VMEC toroidal flux knots, 90 total
		#s = linspace(0, 1, 90)
		s = linspace(0, 0.8, 45)					# 1/2 of knots between 0 and 0.8
		s = append(s[0:-1], linspace(0.8, 1, 46))	# 1/2 of knots between 0.8 and 1.0
		
		# get normalized poloidal flux psi that matches VMEC toroidal flux s, with s = 1 <=> psi = fluxLimit
		f_psitor = scinter.UnivariateSpline(self.ep.PSIdict['psiN1D'], self.psitordic['psitorN1D'], s = 0)	# EFIT based conversion function
		psi = s * self.fluxLimit							# limit normalized poloidal flux to fluxLimit
		x = f_psitor(psi) / f_psitor(self.fluxLimit)		# x is renormalized (x = 0 -> 1) toroidal flux of psi = 0 -> fluxLimit
		f_psiN = scinter.UnivariateSpline(x, psi, s = 0)	# new conversion function, based on x
		psi = f_psiN(s)										# normalized poloidal flux, matching VMEC s
		
		# current profile
		Itor = f_Itor(psi)
		if(self.dic['akima_ac_flag'] == 'density'):
			ac_f = deriv(Itor, s)	# current density: js = dItor/ds
		else:
			ac_f = Itor				# integrated current: Itor = integral(js * ds)

		# pressure and iota profiles
		am_f = self.ep.PROFdict['pfunc'](psi)
		ai_f = 1.0/self.ep.PROFdict['qfunc'](psi)		
		profiles = {'s':s, 'ac_f':ac_f, 'am_f':am_f, 'ai_f':ai_f}
		
		# Polynomial
		if not (self.fitFLAG == 'akima'):
			profiles['AC'] = polyfit(s, js, int(Terms[0:2]) - 1)[::-1]
			profiles['AI'] = polyfit(s, ai_f, int(Terms[2:4]) - 1)[::-1]
			profiles['AM'] = polyfit(s, am_f, int(Terms[4:6]) - 1)[::-1]

		return profiles


	# --- get Boundary Fourier series ---
	def get_bndy(self):
		print 'Getting boundary Fourier series',
		if self.dic['useDESCUR']:
			print '... from DESCUR'
			# choosing the proper angle distribution is critical			
			# use equal arc poloidal angle
#			th = linspace(0, 2*pi, 1000)
			th = linspace(0, 2*pi, 100000)
			R, Z = self.ep.flux_surface(self.fluxLimit, len(th), th)
			Rth, Zth = deriv(R,th), deriv(Z,th)
			integ = antiderivative(sqrt(Rth**2 + Zth**2), th)
			integ /= abs(integ).max()
			thp = 2*pi*integ
			f_th = scinter.UnivariateSpline(thp, th, s = 0)
		
			# get surface
#			N_th = 257
			N_th = 10000

			theta = f_th(linspace(0, 2*pi, N_th + 1))[0:-1]
			R_bndy, Z_bndy = self.ep.flux_surface(self.fluxLimit, N_th, theta)

			# run DESCUR: write Bndy to file, which is input to xcurve -> call xcurve
			out = self.run_descur(R_bndy, Z_bndy)

#			# locate Fourier coeffs in xcurve output
#			for i in xrange(len(out)):
#				if (len(out[i]) > 0):
#					if (out[i][0] == 'RBC(0,0)'):
#						idx = i
#						break

			# get boundary Fourier coefficients 
###			MPOL = 20			# descur returns 20 modes at max
			MPOL = self.dic['MPOL']		# descur returns 20 modes at max
			RBC = zeros(MPOL)
			RBS = zeros(MPOL)
			ZBC = zeros(MPOL)
			ZBS = zeros(MPOL)

			for j in xrange(MPOL):
#				i = idx + j
				coeff = genfromtxt(self.path+'/coeff.dat')
				RBC[j] = coeff[j][0]
				RBS[j] = coeff[j][1]
				ZBC[j] = coeff[j][2]
				ZBS[j] = coeff[j][3]
	
		else:
			out = self.run_xgtovmi()
			
			# locate RBC(0,0)
			for i in xrange(len(out)): 
				if (out[i][0] == 'RBC(0,0)'): 
					idx = i
					break
			
			# get boundary Fourier coefficients 
			MPOL = 20	# min([self.dic['MPOL'], 20])	# xgtovmi returns 20 modes at max
			RBC = zeros(MPOL)
			RBS = zeros(MPOL)
			ZBC = zeros(MPOL)
			ZBS = zeros(MPOL)
	
			for j in xrange(MPOL):
				i = idx + 2*j
				RBC[j] = out[i][2]
				RBS[j] = out[i][5]
				ZBC[j] = out[i+1][2]
				ZBS[j] = out[i+1][5]

		return {'RBC':RBC, 'RBS':RBS, 'ZBC':ZBC, 'ZBS':ZBS}
		
	
	# --- run & read from Fortran DESCUR: xcurve ---
	def run_descur(self, R_bndy, Z_bndy):
		from subprocess import call, STDOUT
		
		# write boundary to a temporary file
		N = len(R_bndy)
		with open(self.path + '/bndy.dat', 'w') as f:
			f.write(str(N) + ' 1 1\n')
			for i in xrange(N):
				f.write(format(R_bndy[i], ' 13.7E') + ' ' + format(Z_bndy[i], ' 13.7E') + '\n')
				
		os.system('python test.py')
		# run it: xcurve is interactive, so use pipe to pass in answers to input questions
		#FNULL = open(os.devnull, 'w')	# be quiet
		#call('echo -e "4\nV\n0\nbndy.dat\nn\n" | xcurve', shell = True, cwd = self.path, stdout = FNULL, stderr = STDOUT)
		
		# read it
#		with open(self.path + '/outcurve', 'r') as f:
#			out = []
#			for line in f:		# read to end of file
				#Modified drh
#				temp = line.replace("(  ","(").replace(",  ",",").replace("="," = ").replace(", R"," ,R").replace(", Z"," ,Z").replace(", ",",")
#				out.append(temp.strip().split())
				
		# garbage collect
###		call(['rm', 'bndy.dat', 'outcurve', 'plotout'], cwd = self.path)
	
#		return out
	
	
	# --- run & read from Fortran xgtovmi ---
	def run_xgtovmi(self):
		from subprocess import call, STDOUT
		print '... from xgtovmi'
		
		# copy a- an g-file if not already there
		remove_gfile = False
		remove_afile = False
		if not os.path.isfile(self.path + '/' + self.gfile ):
			call(['cp', self.gpath + '/' + self.gfile , '.'], cwd = self.path)
			remove_gfile = True
		if not os.path.isfile(self.path + '/' + self.afile['name'] ):
			call(['cp', self.gpath + '/' + self.afile['name'] , '.'], cwd = self.path)
			remove_afile = True
		
		# run it
		FNULL = open(os.devnull, 'w')	# be quiet
		call(['xgtovmi', self.gfile, str(self.fluxLimit), '/cps', self.Terms], cwd = self.path, stdout = FNULL, stderr = STDOUT)
		
		# read it
		with open(self.path + '/' + self.filename, 'r') as f:
			out = []
			for line in f:		# read to end of file
				out.append(line.strip().split())
				
		# garbage collect
		if remove_gfile:
			call(['rm', self.gfile], cwd = self.path)
		if remove_afile:
			call(['rm', self.afile['name']], cwd = self.path)
		call(['rm', 'plotout', 'pgplot.ps', 'map_' + self.gfile + '.nc', 'fort.' + self.gfile[1::]], cwd = self.path)
		call(['rm', 'input.' + self.fitFLAG + '.' + self.gfile[1::]], cwd = self.path)
		
		return out


	# --- evaluate Boundary Fourier series ---
	def ev_bndy(self, u, bndy):
		N = len(u)
		MPOL = len(bndy['RBC'])
		xm = arange(MPOL)
		R = zeros(N)
		Z = zeros(N)
		for i in xrange(N):
			R[i] = sum(bndy['RBC'] * cos(xm*u[i]) + bndy['RBS'] * sin(xm*u[i]))
			Z[i] = sum(bndy['ZBC'] * cos(xm*u[i]) + bndy['ZBS'] * sin(xm*u[i]))

		return R, Z


	# --- get enclosed flux at the edge ---
	def get_phiedge(self, bndy, N = 256):
		print 'Computing enclosed flux'

		# get surface
		u = linspace(0, 2*pi, N + 1)[0:-1]
		R_bndy, Z_bndy = self.ev_bndy(u, bndy)
		r_bndy = self.ep.__get_r__(R_bndy, Z_bndy)
		theta = self.ep.__get_theta__(R_bndy, Z_bndy)

		# evaluate integral(Btor * r * dr)
		intR = zeros(N + 1)
		for i in xrange(N):
			r = append(arange(0, r_bndy[i], 0.0001), r_bndy[i])
			R = r * cos(theta[i]) + self.ep.rmaxis
			Z = r * sin(theta[i]) + self.ep.zmaxis
			Btor = self.ep.BtFunc.ev(R,Z)
			intR[i] = simps(Btor*r, r)
		
		theta = append(theta, theta[0])			# make theta a full 2pi circle
		intR[-1] = intR[0]						# periodicity!!!			
		idx = where(abs(diff(theta)) > 5)[0]	# locate 2pi jump in theta
		theta[idx + 1 ::] += 2*pi				# make theta increase monotonically
		phiedge = simps(intR, theta)			# integral(integral(Btor * r * dr) dtheta)

		return phiedge


	# --- write VMEC input file ---
	def write(self, phiedge, coils, profiles, bndy):
		dic = self.dic
		if not (self.tag == None):
			self.filename += '_' + str(self.tag)
		
		#Fcoils = coils['Fcoils']
		#Ccoils = coils['Ccoils']
		#Icoils = coils['Icoils']
		#eca = coils['eca']
		#ecb = coils['ecb']
		#Bcoil = coils['Bcoil']
	
		RBC = bndy['RBC']
		RBS = bndy['RBS']
		ZBC = bndy['ZBC']
		ZBS = bndy['ZBS']
	
		print 'Writing input file ...',
		with open(self.path + '/' + self.filename, 'w') as f:
			f.write(" &INDATA\n")
			f.write("!----- Runtime Parameters -----------------------------------------------" + "\n")
			#if dic['use_OMP']: 
			#	f.write("  OMP_NUM_THREADS = " + str(dic['OMP_NUM_THREADS']) + "\n")
			f.write("  DELT = 0.25\n")
			#f.write("  TCON0 = 2.0\n")
			f.write("  LFORBAL = " + dic['LFORBAL'] + "\n")
			f.write("  NSTEP = " + str(dic['NSTEP']) + "\n")		
			f.write("  NS_ARRAY = " + ' '.join([str(item) for item in dic['NS_ARRAY']]) + "\n")
			f.write("  FTOL_ARRAY = " + ' '.join([str(item) for item in dic['FTOL_ARRAY']]) + "\n")
			f.write("  NITER_ARRAY = " + ' '.join([str(item) for item in dic['NITER_ARRAY']]) + "\n")
			if dic['use_PRECON']:
				f.write("  PRECON_TYPE = '" + dic['PRECON_TYPE'] + "'\n")
				f.write("  PREC2D_THRESHOLD = " + str(dic['PREC2D_THRESHOLD']) + "\n")
			else:
				f.write("!  PRECON_TYPE = '" + dic['PRECON_TYPE'] + "'\n")
				f.write("!  PREC2D_THRESHOLD = " + str(dic['PREC2D_THRESHOLD']) + "\n")

			f.write("!----- Grid Parameters --------------------------------------------------" + "\n")
			f.write("  LASYM = " + dic['LASYM'] + "\n")
			#f.write("  LRFP = " + dic['LRFP'] + "\n")
			f.write("  NFP = " + str(dic['NFP']) + "\n")		
			#f.write("  NTOR = " + str(dic['NTOR']) + "\n")
			#f.write("  NTHETA = " + str(dic['NTHETA']) + "\n")
			#f.write("  NZETA = " + str(dic['NZETA']) + "\n")
			#if dic['use_APHI']:
			#	f.write("  APHI = " + ' '.join([str(item) for item in dic['APHI']]) + "\n")
			#else:
			#	f.write("!  APHI = " + ' '.join([str(item) for item in dic['APHI']]) + "\n")
			f.write("  PHIEDGE = " + str(phiedge) + "\n")

			f.write("!----- Boundary Parameters ----------------------------------------------" + "\n")
			f.write("  LFREEB = " + dic['LFREEB'] + "\n")
			#f.write("  MGRID_FILE = '" + str(dic['mgrid']) + "'\n")
			f.write("  NVACSKIP = " + str(dic['NVACSKIP']) + "\n")

			f.write("!----- All Currents -----------------------------------------------------" + "\n")
			f.write("  CURTOR = " + str(self.ep.data.get('current')) + "\n")
			#f.write("  !-- F-coils --\n")
			#f.write("  EXTCUR(01) = " + format(Fcoils[0],' 13.7E') + "  EXTCUR(02) = " + format(Fcoils[9],' 13.7E') + "\n")  #  1: f1a  2: f1b
			#f.write("  EXTCUR(03) = " + format(Fcoils[1],' 13.7E') + "  EXTCUR(04) = " + format(Fcoils[10],' 13.7E') + "\n") #  3: f2a  4: f2b 
			#f.write("  EXTCUR(05) = " + format(Fcoils[2],' 13.7E') + "  EXTCUR(06) = " + format(Fcoils[11],' 13.7E') + "\n") #  5: f3a  6: f3b 
			#f.write("  EXTCUR(07) = " + format(Fcoils[3],' 13.7E') + "  EXTCUR(08) = " + format(Fcoils[12],' 13.7E') + "\n") #  7: f4a  8: f4b 
			#f.write("  EXTCUR(09) = " + format(Fcoils[4],' 13.7E') + "  EXTCUR(10) = " + format(Fcoils[13],' 13.7E') + "\n") #  9: f5a 10: f5b 
			#f.write("  EXTCUR(11) = " + format(Fcoils[5],' 13.7E') + "  EXTCUR(12) = " + format(Fcoils[14],' 13.7E') + "\n") # 11: f6a 12: f6b 
			#f.write("  EXTCUR(13) = " + format(Fcoils[6],' 13.7E') + "  EXTCUR(14) = " + format(Fcoils[15],' 13.7E') + "\n") # 13: f7a 14: f7b 
			#f.write("  EXTCUR(15) = " + format(Fcoils[7],' 13.7E') + "  EXTCUR(16) = " + format(Fcoils[16],' 13.7E') + "\n") # 15: f8a 16: f8b 
			#f.write("  EXTCUR(17) = " + format(Fcoils[8],' 13.7E') + "  EXTCUR(18) = " + format(Fcoils[17],' 13.7E') + "\n") # 17: f9a 18: f9b 
			#f.write("  !-- B-coils --\n")
			#if 'scaled' in str(dic['mgrid']): 
			#	f.write("  EXTCUR(19) = " + format(abs(Bcoil),' 13.7E') + "\n")
			#else:
			#	f.write("  EXTCUR(19) = " + format(Bcoil,' 13.7E') + "\n")
			#f.write("  !-- E-coils --\n")
			#f.write("  EXTCUR(20) = " + format(eca,' 13.7E') + "  EXTCUR(21) = " + format(ecb,' 13.7E') + "\n")				  # 20: eca  21: ecb
			#f.write("  !-- I-coils --\n")
			#f.write("  EXTCUR(22) = " + format(Icoils[5],' 13.7E') + "  EXTCUR(23) = " + format(Icoils[4],' 13.7E') + "\n")   # 22: iu330    23: iu270
			#f.write("  EXTCUR(24) = " + format(Icoils[3],' 13.7E') + "  EXTCUR(25) = " + format(Icoils[2],' 13.7E') + "\n")   # 24: iu210    25: iu150
			#f.write("  EXTCUR(26) = " + format(Icoils[1],' 13.7E') + "  EXTCUR(27) = " + format(Icoils[0],' 13.7E') + "\n")   # 26: iu90     27: iu30
			#f.write("  EXTCUR(28) = " + format(Icoils[11],' 13.7E') + "  EXTCUR(29) = " + format(Icoils[10],' 13.7E') + "\n") # 28: il330    29: il270
			#f.write("  EXTCUR(30) = " + format(Icoils[9],' 13.7E') + "  EXTCUR(31) = " + format(Icoils[8],' 13.7E') + "\n")   # 30: il210    31: il150
			#f.write("  EXTCUR(32) = " + format(Icoils[7],' 13.7E') + "  EXTCUR(33) = " + format(Icoils[6],' 13.7E') + "\n")   # 32: il90     33: il30
			#f.write("  !-- C-coils --\n")
			#f.write("  EXTCUR(34) = " + format(Ccoils[0],' 13.7E') + "  EXTCUR(35) = " + format(Ccoils[1],' 13.7E') + "\n")   # 34: ccoil79  35: ccoil139
			#f.write("  EXTCUR(36) = " + format(Ccoils[2],' 13.7E') + "\n") 													  # 36: ccoil199

			f.write("!----- Profiles ---------------------------------------------------------" + "\n")
			f.write("  NCURR = " + str(dic['NCURR']) + "\n")
		
			if(self.fitFLAG == 'akima'):
				f.write("  pmass_type = 'Akima_spline'\n")
				output = write_array('am_aux_s', profiles['s'])
				for line in output: f.write(line)
				output = write_array('am_aux_f', profiles['am_f'])
				for line in output: f.write(line)
		
				if(dic['akima_ac_flag'] == 'density'):
					f.write("  pcurr_type = 'Akima_spline_Ip'\n")
				else:
					f.write("  pcurr_type = 'Akima_spline_I'\n")
				output = write_array('ac_aux_s', profiles['s'])
				for line in output: f.write(line)
				output = write_array('ac_aux_f', profiles['ac_f'])
				for line in output: f.write(line)
		
				f.write("  piota_type = 'Akima_spline'\n")
				output = write_array('ai_aux_s', profiles['s'])
				for line in output: f.write(line)
				output = write_array('ai_aux_f', profiles['ai_f'])
				for line in output: f.write(line)
			else:
				output = write_array('AI', profiles['AI'])
				for line in output: f.write(line)

				output = write_array('AC', profiles['AC'])
				for line in output: f.write(line)

				output = write_array('AM', profiles['AM'])
				for line in output: f.write(line)

			f.write("!----- Boundary ---------------------------------------------------------" + "\n")
			f.write("  MPOL = " + str(dic['MPOL']) + "\n")
			for m in xrange(len(RBC)):
				f.write("  RBC(0," + str(m) + ") = " + format(RBC[m],' 13.7E') + "  RBS(0," + str(m) + ") = " + format(RBS[m],' 13.7E') + "\n")
				f.write("    ZBC(0," + str(m) + ") = " + format(ZBC[m],' 13.7E') + "  ZBS(0," + str(m) + ") = " + format(ZBS[m],' 13.7E') + "\n")
			f.write("  RAXIS = " + str(self.ep.rmaxis) + "\n")
			f.write("  ZAXIS = " + str(self.ep.zmaxis) + "\n")

			# finish the input file
			f.write("/\n")
			f.write("/\n")
			f.write("!----- Final comments ---------------------------------------------------" + "\n")
			f.write("! mapcode boundary fraction is " + str(self.fluxLimit) + "\n")
			f.write("! file " + self.filename + " generated by\n")
			f.write("! g2vmi for " + getpass.getuser() + " on " + strftime("%c") + "\n")
			#f.write("! Coil ordering for VMEC differs from that for EFIT\n")
			#f.write("!      1:      f1a      2:     f1b\n")
			#f.write("!      3:      f2a      4:     f2b\n")
			#f.write("!      5:      f3a      6:     f3b\n")
			#f.write("!      7:      f4a      8:     f4b\n")
			#f.write("!      9:      f5a      10:    f5b\n")
			#f.write("!      11:     f6a      12:    f6b\n")
			#f.write("!      13:     f7a      14:    f7b\n")
			#f.write("!      15:     f8a      16:    f8b\n")
			#f.write("!      17:     f9a      18:    f9b\n")
			#f.write("!      19:     1oR\n")
			#f.write("!      20:     eca      21:    ecb\n")
			#f.write("!      22:     iu330    23:    iu270\n")
			#f.write("!      24:     iu210    25:    iu150\n")
			#f.write("!      26:     iu90     27:    iu30\n")
			#f.write("!      28:     il330    29:    il270\n")
			#f.write("!      30:     il210    31:    il150\n")
			#f.write("!      32:     il90     33:    il30\n")
			#f.write("!      34:     ccoil79\n")
			#f.write("!      35:     ccoil139\n")
			#f.write("!      36:     ccoil199\n")
			f.write("! initial call used to create this file\n")
			f.write("! g2vmi.py " + self.gpath + '/' + self.gfile + " " + str(self.fluxLimit) + " -t " + self.Terms + "\n")
	
		print 'wrote VMEC input file to dir: ' + self.path + '/'
	
	
	# --- plot boundary and input profiles ---
	def plot_input(self):
		print 'Plotting input ...'
		import VMEC.Python.plotScripts as ps
		import pylab as plt
		plt.ion()
		ps.plot_input_all(self.path + '/' + self.filename, self.gpath + '/' + self.gfile, FluxSurfList = self.FluxSurfList)
		plt.ioff()
		plt.show()


# ----------------------------------------------------------------------------------------
# --- End of class -----------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# --- F = antiderivative(y, x) -----------------------------------------------------------
# calculates the antiderivative F of y with respect to x at each point x[i], i = 0,...,N-1
# so F'(x) = y(x) with F(x[0]) = 0 (Note, this fixes the intergartional constant)
# x does not need to be equally spaced; algorithm uses trapezoid rule 
# F[-1] = F[N-1] = \int\limits_{x[0]}^(x[N-1]) y(x) dx
def antiderivative(y, x):
    F = zeros(x.shape)
    for i in xrange(1, x.size): 
        F[i] = F[i-1] + 0.5*(x[i] - x[i-1])*(y[i] + y[i-1])
            
    return F


# ----------------------------------------------------------------------------------------
# --- write_array ------------------------------------------------------------------------
# puts array into a formated string list, suitable for writing to file
def write_array(name, s, cols = 5):
	N = len(s)
	output = []
	output.append("  " + name + "  =  " + format(s[0],' 13.7E') + "\n")
	idx = 1
	while (idx < N):
		line =  '    ' + format(s[idx],' 13.7E')
		idx += 1
		for k in xrange(cols -1):
			if(idx < N): 
				line += '  ' + format(s[idx],' 13.7E')
			idx += 1
		output.append(line + '\n')
	return output


# ----------------------------------------------------------------------------------------
# --- split_data -------------------------------------------------------------------------
# In the g-file there is typically no space between data points,
# if the value is negative (space is occupied by the '-' sign).
# Each row in the file is read in as a single string.
# This routine searches for the 'e' in each floating point value
# and splits the string 'line' accordingly
def split_data(line):
	length = len(line)
	index = []

	# find all the 'E' in the string, each number is in floating representation and therefore has one
	for j in range(0,length):
		if ((line[j] == 'E') | (line[j] == 'e')) & (j < length-4):
			index.append(j+4)
		
	# Split up line into fragments using the positions of 'E'
	line3 = []
	line3.append(line[0:index[0]].strip())  # first data# .strip() removes any blanks in front 
	for k in range(0,len(index)-1):
		line3.append(line[index[k]:index[k+1]].strip()) # omitt line[index[-1]:length] = '\n', so stop before 
	return line3


# ----------------------------------------------------------------------------------------
# --- read_diiidsup ----------------------------------------------------------------------
# reads the diiidsup.in file in path and returns the coil currents
def read_diiidsup(path = './'):
	with open(path + 'diiidsup.in', 'r') as f:
		# skip 4 lines
		for i in xrange(4): head = f.readline()
		
		# Fcoils
		Fcoils = []
		for i in xrange(2):
			line = f.readline().split()
			Fcoils.append([float(x) for x in line])
		Fcoils = [num for elem in Fcoils for num in elem]
		Fcoils = array(Fcoils)
		
		# skip a line
		head = f.readline()
		
		# Ccoils
		Ccoils = []
		line = f.readline().split()
		Ccoils.append([float(x) for x in line])
		Ccoils = [num for elem in Ccoils for num in elem]
		Ccoils = array(Ccoils)

		# skip a line
		head = f.readline()
		
		# Icoils
		Icoils = []
		for i in xrange(2):
			line = f.readline().split()
			Icoils.append([float(x) for x in line])
		Icoils = [num for elem in Icoils for num in elem]
		Icoils = array(Icoils)
	
	return Fcoils, Ccoils, Icoils


# ----------------------------------------------------------------------------------------
# --- replace_akima_current --------------------------------------------------------------
# replaces ac_aux_f in input file with the current density profile from wout file; keeps ac_aux_s the same
# wout (string): wout file name
def replace_akima_current(filename, wout, type = 'density'):
	import VMEC.Python.wout_class as WC
	# read input file
	print 'Modify vmec input file ...'
	with open(filename, 'r') as f:
		input = []
		for line in f:		# read to end of file
			input.append(line)

	# locate pcurr_type, ac_aux_s and ac_aux_f lines in input; akima splines only
	for i in xrange(len(input)):
		if (input[i].strip()[0:5] == 'NCURR'):
			input[i] = "  NCURR = 1\n"					# switch to force current profile
		if(input[i].strip()[0:10] == 'pcurr_type'):
			idx = i
		if (input[i].strip()[0:8] == 'ac_aux_s'):
			s_idx = i
		if (input[i].strip()[0:8] == 'ac_aux_f'):
			f_idx = i
			break

	# read knot locations in s to use same knots as before
	s = [float64(input[s_idx].strip().split()[2])]
	for line in input[s_idx+1:f_idx]:
		line = line.strip().split()
		for x in line:
			s.append(float64(x))
	s = array(s)
	N = len(s)
	
	# set current flag
	if(type == 'density'): 
		input[idx] = "  pcurr_type = 'Akima_spline_Ip'\n"
	else: 
		input[idx] = "  pcurr_type = 'Akima_spline_I'\n"

	# get current profile from wout and interpolate
	w = WC.Wout(wout).getAll()
	jcurv = w['spljcurv'](s)

	# current profile
	js = jcurv * 2*pi
	if(type == 'density'):
		f = js								# current density: js
	else:
		f = append(0, cumtrapz(js, s))		# integrated current: I = integral(js * ds)

	# reset f values
	output = write_array('ac_aux_f', f)
	for i, line in enumerate(output): 
		input[f_idx + i] = line
		
	# write file
	with open(filename, 'w') as f:
		for i in xrange(len(input)): f.write(input[i])


# ----------------------------------------------------------------------------------------
# --- replace_akima_iota -----------------------------------------------------------------
# replaces ai_aux_f in input file with the iota profile from wout file; keeps ai_aux_s the same
# if LFREEB == True: switch input file to free bndy mode and rename it, replace keyword 
# 'fixed' or 'fixd' with 'free' in filename
# wout (string): wout file name
def replace_akima_iota(filename, wout, LFREEB = False):
	import VMEC.Python.wout_class as WC
	# read input file
	print 'Modify vmec input file ...'
	with open(filename, 'r') as f:
		input = []
		for line in f:		# read to end of file
			input.append(line)

	# locate ai_aux_s line in original input; akima splines only
	for i in xrange(len(input)):
		if (input[i].strip()[0:6] == 'LFREEB') & LFREEB:
			input[i] = "  LFREEB = T\n"				# switch to free boundary mode
		if (input[i].strip()[0:5] == 'NCURR'):
			input[i] = "  NCURR = 0\n"				# switch to force iota profile
		if (input[i].strip()[0:8] == 'ai_aux_s'):
			s_idx = i
		if (input[i].strip()[0:8] == 'ai_aux_f'):
			f_idx = i
			break

	# read knot locations in s to use same knots as before
	s = [float64(input[s_idx].strip().split()[2])]
	for line in input[s_idx+1:f_idx]:
		line = line.strip().split()
		for x in line:
			s.append(float64(x))
	s = array(s)
	N = len(s)
	
	# get iota profile from wout and interpolate
	w = WC.Wout(wout).getAll()
	f = w['spliotaf'](s)
	
	# reset f values
	output = write_array('ai_aux_f', f)
	for i, line in enumerate(output): 
		input[f_idx + i] = line
	
	# write file
	if LFREEB:
		if('fixed' in filename): filename = filename.replace('fixed', 'free')
		elif('fixd' in filename): filename = filename.replace('fixd', 'free')
	with open(filename, 'w') as f:
		for i in xrange(len(input)): f.write(input[i])


# ----------------------------------------------------------------------------------------
# --- main -------------------------------------------------------------------------------
def main(gfile, fluxLimit, Terms = '000000', path = '.', Time = None, plotit = False, tag = None):
	GTI = g2vmi(gfile, fluxLimit, Terms, path = path, Time = Time, tag = tag)

	# get coil currents
	Ccoils, Icoils = GTI.get_CIcoils()
	Fcoils = GTI.get_Fcoils()
	eca, ecb = GTI.get_Ecoils()
	Bcoil = GTI.get_Bcoils()
	coils = {'Ccoils':Ccoils, 'Icoils':Icoils, 'Fcoils':Fcoils, 'eca':eca, 'ecb':ecb, 'Bcoil':Bcoil}
	
	# fit profiles
	profiles = GTI.fit_profiles()
	
	# get boundary Fourier series
	bndy = GTI.get_bndy()
	
	# get enclosed toroidal flux
	phiedge = GTI.get_phiedge(bndy)

	# write VMEC input file
	GTI.write(phiedge, coils, profiles, bndy)
	
	# check input profiles & boundary
	if plotit:
		GTI.plot_input()


# ----------------------------------------------------------------------------------------
# --- Launch main() ----------------------------------------------------------------------
if __name__ == '__main__':
	import argparse
	import textwrap
	parser = argparse.ArgumentParser(description = 'Creates a VMEV input file', 
				formatter_class = argparse.RawDescriptionHelpFormatter,
				epilog = textwrap.dedent('''\
                Examples: g2vmi.py g148712.04101 0.9995
                          g2vmi.py /path/to/gfile/g148712.04101 0.9995
                          g2vmi.py ./g148712.04101 0.9995 -tag hallo'''))

	parser.add_argument('gfile', help = 'EFIT geqdsk file name or (full or rel.) pathname', type = str)
	parser.add_argument('fluxLimit', help = 'percent of poloidal flux surface to use as VMEC outer boundary', type = float)
	parser.add_argument('-t', '--terms', help = 'number of polynomial terms for AC,AI,AM; use 00 for Akima splines; default is 000000', type = str, default = '000000')
	parser.add_argument('-tag', help = 'arbitrary tag, added to input file name; default is None', type = str, default = None)
	parser.add_argument('-p', '--plotit', help = 'plot input', action = 'store_true', default = False)
	parser.add_argument('-time', help = 'time for I- & C-coil currents; default is g-file time; ignored, if diiidsup.in file already exists', type = int, default = None)
	parser.add_argument('-path', help = 'path to write VMEC input file to; default is ./', type = str, default = '.')
	args = parser.parse_args()
	
	main(args.gfile, args.fluxLimit, Terms = args.terms, path = args.path, Time = args.time, plotit = args.plotit, tag = args.tag)

	
