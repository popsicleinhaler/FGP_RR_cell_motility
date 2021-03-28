from auto import *

#B = [2,4.26,10,20,25,30,40,50,100]

f = load('chem_Dawn_Kathy')

pars = ['kon1','koff1','kon2','koff2','kfs','L','n','B','I_rho','L_rho','delta_rho','L_R','I_R','delta_R','alpha_R','delta_P','I_K','k_X','k_G','k_C','PIX','Paxtot','m','alpha_PAK','Rho_Square','Rac_Square','Pax_Square','PAKtot']
#        4lp    none     2lp    none  none  4LP na  na    2LPs   2 rev    2rev       2lps  2rev   2lps      3lps     4lps      2rev  2rev 2rev  2rev   2rev   2rev   none   2rev        2rev           2lps       2lps       nothing -> ran longer: 4lps

#       'kon1'    'koff1'   'kon2'  'koff2'  'kfs'      'L'     'n'    'B'   'I_rho' 'L_rho' 'delta_rho' 'L_R'    'I_R' 'delta_R','alpha_R','delta_P',  'I_K',    'k_X',   'k_G'   'k_C'   'PIX'   'Paxtot'  'm'  'alp_PAK''Rh_Sq'  'Ra_Sq'  'Pax_Sq'  'PAKtot'
bmx = [[1759,200],[9,700],[393,700],[9,1310],[9,4629],[300,1000],[1,1],[9,220],[40,60],[9,110],  [9,47],  [60,105],[9,60],[30,1]    ,[9,800], [170,145],  [9,2206],[9,1504],[9,2343],[9,490],[9,153], [9,115], [1,1], [9,65],[9,124],[60,1783],[140,134], [8400,1]]

def fse(plot=False, par2='kfs', lpnum=4, to2=False, par2v=5):
	i = pars.index(par2)

	if to2:
		# becareful of when par2v is smaller than the default
		r0 = run(f,ICP=[par2],UZSTOP={par2:par2v}) + run(f, DS='-', ICP=[par2], NMX=220)
		#r0 = run(f,UZR={'delta_P':0.04},ICP=['delta_P'],NMX=200) + run(f,UZR={'delta_P':0.04}, DS='-', ICP=['delta_P'],NMX=200)
		r0 = run(r0('UZ1'),ICP=['delta_P'],NMX=500) + run(r0('UZ1'), DS='-', ICP=['delta_P'],NMX=200)

	#r0 = run(f, DS='-', ICP=['B'],UZSTOP={'B':b})
	#\,PAR={8:b})
	#pl(r0, stability=True,bifurcation_y='Rac', top_title='1par in B')
	#wait()
	else:
		r0 = run(f,ICP=['delta_P'],NMX=170) + run(f, DS='-', ICP=['delta_P'],NMX=170)
		#\,PAR={8:b})

	if plot:
		pl(r0, stability=True,bifurcation_y='Rac', top_title='1par in delta_P')
		wait()
	
	lps = r0('LP')
	r2f = []
	r2r = []
	counter = 0
	for lp in lps:
		if counter > (lpnum-1):
			break
		r2f += run(lp,ISW=2, ICP=['delta_P',par2], NMX=200 , SP=['LP0','HB0','BP0'])
		#\,PAR={8:b})
	#plot(r2f,bifurcation_y='kfs',top_title='forward LPs in r2')
	#wait()
		r2r += run(lp,ISW=2, DS='-',ICP=['delta_P',par2], NMX=700, SP=['LP0','HB0','BP0'])
		#\,PAR={8:b})
		counter += 1
	#plot(r2r,bifurcation_y='kfs',top_title='bakcward LPs in r2')
	#print(r2f+r2r)
	#wait()
	#plot(r2r+r2f,bifurcation_y='kfs',top_title='cont LPs in kfs and delta_P')
	return r0,r2f,r2r


def sef(b=2):
	if b >= 2:
		if b==2:
			b += 0.01
		
		#pl(r0, stability=True,bifurcation_y='Rac', top_title='1par with B = '+str(b))
		#wait()
	if b < 2:
		r0 = run(f,ICP=['B'],NMX=220) + run(f, DS='-', ICP=['B'], UZSTOP={'B':b})
		r0 = run(r0('UZ1'),ICP=['delta_P'],NMX=500) + run(r0('UZ1'), DS='-', ICP=['delta_P'],NMX=200)
		#pl(r0, stability=True,bifurcation_y='Rac', top_title='1par with B = '+str(b))
		#wait()
	lps = r0('LP')
	r2f = []
	r2r = []
	counter = 0
	for lp in lps:
		if counter > 3:
			break
		r2f += run(lp,ISW=2, ICP=['delta_P','kfs','B'], UZSTOP={'kfs':100}, SP=['LP0','HB0','BP0'])
		r2r += run(lp,ISW=2, DS='-',ICP=['delta_P','kfs','B'], SP=['LP0','HB0','BP0'],NMX=700)
		counter += 1
	##print(r2f+r2r)
	#wait()
	#plot(r2r+r2f+r0,bifurcation_y='kfs',top_title='cont LPs in kfs and delta_P for B:'+str(b))
	return r2f+r2r


	#\,PAR={8:b})
	#pl(r0, stability=True,bifurcation_y='Rac', top_title='1par in B')
	#wait()
	#
	# 25 -> f:500, r:9881






	#
		#\,PAR={8:b})

	#plot(r2f,bifurcation_y='kfs',top_title='forward LPs in r2')
	#wait()


	#plot(r2r,bifurcation_y='kfs',top_title='bakcward LPs in r2')







def par3lp():
	r0,r2f,r2r = fse()
	lps = r0('LP')


	r2f = []
	r2r = []
	counter = 0
	for lp in lps:
		if counter > 3:
			break
		r2f += run(lp,ISW=2, ICP=['delta_P','kfs','B'], UZSTOP={'kfs':100}, SP=['LP0','HB','BP0'])
		r2r += run(lp,ISW=2, DS='-',ICP=['delta_P','kfs','B'], SP=['LP0','HB','BP0'],NMX=700)
		counter += 1
	return r2r+r2f

def par3cp():

	r0,r2f,r2r = fse()
	r2 = r2f + r2r
	cusps = r2('CP')
	cp3 = []
	counter = 0
	for cp in cusps:
		if counter > 0:
			break
		print(cp)
		wait()
		cp3 += run(cp,ISW=3, ICP=['delta_P','kfs','B']) + run(cp,ISW=3, DS='-',ICP=['delta_P','kfs','B'])
		counter += 1
	plot(cp3)
	return cp3


	#r2 = \
	#run(r0('LP1'),ISW=2, ICP=['delta_P','kfs'],UZSTOP={'delta_P':0,'kfs':0})+\
	#run(r0('LP1'),ISW=2, DS='-',ICP=['delta_P','kfs'],UZSTOP={'delta_P':0.1,'kfs':3})+\
	#run(r0('LP2'),ISW=2, ICP=['delta_P','kfs'],UZSTOP={'delta_P':0,'kfs':0})+\
	#run(r0('LP2'),ISW=2, DS='-',ICP=['delta_P','kfs'],UZSTOP={'delta_P':0.1,'kfs':3})+\
	#run(r0('LP3'),ISW=2, ICP=['delta_P','kfs'],UZSTOP={'delta_P':0.1,'kfs':0})+\
	#run(r0('LP3'),ISW=2, DS='-',ICP=['delta_P','kfs'],UZSTOP={'delta_P':-0.7,'kfs':3})+\
	#run(r0('LP4'),ISW=2, ICP=['delta_P','kfs'],UZSTOP={'delta_P':0.2,'kfs':0})+\
	#run(r0('LP4'),ISW=2, DS='-',ICP=['delta_P','kfs'],UZSTOP={'delta_P':-0.7,'kfs':3})

def lpcont(r0,par1,par2):
	lps = r0('LP')
	r2f = []
	r2r = []
	counter = 0
	for lp in lps:
		if counter > 3:
			break
		r2f += run(lp,ISW=2, ICP=[par1,par2], SP=['LP0','HB0','BP0'])
		r2r += run(lp,ISW=2, DS='-',ICP=[par1,par2], SP=['LP0','HB0','BP0'])
		counter += 1
	return r2f,r2r


# def changeB(b):
# 	return


#----------------------------for endpoint continuation


# 	epf = r2f
# 	epr = r2r

# 	#print(epr)
# 	#print(epf)
# 	#wait()

# 	i = 0
# 	n = 7

# 	while i < n:

# 		try:
# 			endpointf = epf('EP')
# 			epf = []
# 			for endf in endpointf:
# 				epf += run(endf,ISW=2, ICP=['delta_P','kfs'],PAR={8:b}, NMX=10000000, SP=['LP0','HB0','BP0'])
# 			r2f += epf
# 			# print(epf)
# 			# wait()
# 		except TypeError:
# 			print('probably MX in ALL brsnches or no branches?')
# 			wait()
# 			pass

		

# 		#plot(epf,bifurcation_y='kfs',top_title='forward EPs'+str(i))
# 		#wait()

# 		try:
# 			endpointr = epr('EP')
# 			epr = []
# 			for endr in endpointr:
# 				epr += run(endr,ISW=2, ICP=['delta_P','kfs'],PAR={8:b}, NMX=10000000, SP=['LP0','HB0','BP0']) 
# 				#literally fkn magic WHY NO DS='-'??????????????
# 			r2r += epr
# 			# print(epf)
# 			# wait()
# 		except TypeError:
# 			print('maybe MX in all branches or just no branches')
# 			wait()
# 			pass


# 		#plot(epr+epf,bifurcation_y='kfs',top_title='combo EPs'+str(i))
# 		#wait()

# 		i+=1


# 	return r2r+r2f+r0

# #plot(fse(2), bifurcation_y='Rac', top_title='all eps continued for B:2'

# 		#wait()



