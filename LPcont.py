from auto import *

def clps(par1,par2,nxf,nxr,B):

	f = load('chem_Dawn_Kathy')
	r0 = r(f,ICP=[par1],NMX=nxf)+r(f, DS='-', ICP=[par1],NMX=nxr)
	plot(r0, stability=True, bifurcation_y='Rac',top_title='1par')

	r2 = run(r0('LP1'),ISW=2, ICP=[par1,par2]) + run(r0('LP1'),ISW=2, DS='-',ICP=[par1,par2])\
		+run(r0('LP2'),ISW=2, ICP=[par1,par2]) + run(r0('LP2'),ISW=2, DS='-',ICP=[par1,par2])\
		+run(r0('LP4'),ISW=2, ICP=[par1,par2]) + run(r0('LP3'),ISW=2, DS='-',ICP=[par1,par2])\
		+run(r0('LP4'),ISW=2, ICP=[par1,par2]) + run(r0('LP4'),ISW=2, DS='-',ICP=[par1,par2])

	plot(r2, bifurcation_y=par2, bifurcation_z='Rac',top_title='4 LP continuation in:'+par1+' and '+par2, coloring_method= 'branch')

	plot(r0+r2, stability=True, top_tilte='overlaped LPs and 1par no B')

	r2=rl(r2)
	save(r2, '2par bifurcation in :'+par1+'and'+par2)

	

	# r0 = r(f,ICP=[par1],NMX=nxf,PAR={8:B})+r(f, DS='-', ICP=[par1],NMX=nxr,PAR={8:B})
	# plot(r0, stability=True, bifurcation_y='Rac')

	# r2 = run(r0('LP1'),ISW=2, ICP=[par1,par2],PAR={8:B}) + run(r0('LP1'),ISW=2, DS='-',ICP=[par1,par2],PAR={8:B})\
	# 	+run(r0('LP2'),ISW=2, ICP=[par1,par2],PAR={8:B}) + run(r0('LP2'),ISW=2, DS='-',ICP=[par1,par2],PAR={8:B})\
	# 	+run(r0('LP4'),ISW=2, ICP=[par1,par2],PAR={8:B}) + run(r0('LP3'),ISW=2, DS='-',ICP=[par1,par2],PAR={8:B})\
	# 	+run(r0('LP4'),ISW=2, ICP=[par1,par2],PAR={8:B}) + run(r0('LP4'),ISW=2, DS='-',ICP=[par1,par2],PAR={8:B})

	# plot(r2, bifurcation_y=par2, bifurcation_z='Rac',top_title='4 LP continuation in:'+par1+' and '+par2+'with B', coloring_method= 'branch')

	# plot(r0+r2, stability=True, top_tilte='overlaped LPs and 1par with B')

	