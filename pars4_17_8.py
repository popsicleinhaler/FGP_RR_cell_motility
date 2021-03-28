from auto import *

#B = [2,4.26,10,20,25,30,40,50,100]

f = load('chem_Dawn_Kathy')

def fse(b):


	#r0 = run(f,ICP=['B'],NMX=9,PAR={8:b}) + run(f, DS='-', ICP=['B'],NMX=220)
	#\,PAR={8:b})
	#pl(r0, stability=True,bifurcation_y='Rac', top_title='1par in B')
	#wait()

	r0 = run(f,ICP=['delta_P'],NMX=170) + run(f, DS='-', ICP=['delta_P'],NMX=170)
	#\,PAR={8:b})
	pl(r0, stability=True,bifurcation_y='Rac', top_title='1par with B = '+str(b))
	wait()

	lps = r0('LP')

	r2f = []
	r2r = []

	counter = 0

	for lp in lps:

		if counter > 3:
			break

		r2f += run(lp,ISW=2, ICP=['delta_P','kfs'], UZSTOP={'kfs':100}, SP=['LP0','HB0','BP0'])
		#\,PAR={8:b})

	#plot(r2f,bifurcation_y='kfs',top_title='forward LPs in r2')
	#wait()

		r2r += run(lp,ISW=2, DS='-',ICP=['delta_P','kfs'], NMX=700, SP=['LP0','HB0','BP0'])
		#\,PAR={8:b})

		counter += 1

	#plot(r2r,bifurcation_y='kfs',top_title='bakcward LPs in r2')
	print(r2f+r2r)
	wait()
	plot(r2r+r2f,bifurcation_y='kfs',top_title='cont LPs in kfs and delta_P')
	return r2f+r2r+r0



def fse(b):


	#r0 = run(f,ICP=['B'],NMX=9,PAR={8:b}) + run(f, DS='-', ICP=['B'],NMX=220)
	#\,PAR={8:b})
	#pl(r0, stability=True,bifurcation_y='Rac', top_title='1par in B')
	#wait()

	r0 = run(f,ICP=['delta_P'],NMX=170) + run(f, DS='-', ICP=['delta_P'],NMX=170)
	#\,PAR={8:b})
	#pl(r0, stability=True,bifurcation_y='Rac', top_title='1par with B = '+str(b))
	#wait()

	lps = r0('LP')

	r2f = []
	r2r = []

	counter = 0

	for lp in lps:

		if counter > 3:
			break

		r2f += run(lp,ISW=2, ICP=['delta_P','kfs'], UZSTOP={'kfs':100}, SP=['LP0','HB0','BP0'])
		#\,PAR={8:b})

	#plot(r2f,bifurcation_y='kfs',top_title='forward LPs in r2')
	#wait()

		r2r += run(lp,ISW=2, DS='-',ICP=['delta_P','kfs'], NMX=700, SP=['LP0','HB0','BP0'])
		#\,PAR={8:b})

		counter += 1

	#plot(r2r,bifurcation_y='kfs',top_title='bakcward LPs in r2')
	print(r2f+r2r)
	wait()
	plot(r2r+r2f,bifurcation_y='kfs',top_title='cont LPs in kfs and delta_P')
	return r2f+r2r+r0







	#r2 = \
	#run(r0('LP1'),ISW=2, ICP=['delta_P','kfs'],UZSTOP={'delta_P':0,'kfs':0})+\
	#run(r0('LP1'),ISW=2, DS='-',ICP=['delta_P','kfs'],UZSTOP={'delta_P':0.1,'kfs':3})+\
	#run(r0('LP2'),ISW=2, ICP=['delta_P','kfs'],UZSTOP={'delta_P':0,'kfs':0})+\
	#run(r0('LP2'),ISW=2, DS='-',ICP=['delta_P','kfs'],UZSTOP={'delta_P':0.1,'kfs':3})+\
	#run(r0('LP3'),ISW=2, ICP=['delta_P','kfs'],UZSTOP={'delta_P':0.1,'kfs':0})+\
	#run(r0('LP3'),ISW=2, DS='-',ICP=['delta_P','kfs'],UZSTOP={'delta_P':-0.7,'kfs':3})+\
	#run(r0('LP4'),ISW=2, ICP=['delta_P','kfs'],UZSTOP={'delta_P':0.2,'kfs':0})+\
	#run(r0('LP4'),ISW=2, DS='-',ICP=['delta_P','kfs'],UZSTOP={'delta_P':-0.7,'kfs':3})



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



