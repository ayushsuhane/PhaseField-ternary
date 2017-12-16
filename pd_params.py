import numpy as np
import math
from scipy.optimize import fsolve
# import pylab
# from matplotlib.font_manager import FontProperties
# 
# 
# 	
# def equations(p,data):
# 	T_eutectoid = 1.0
# 	x = p
# 	f = np.ones(len(x))
# 	print f,len(x)
# 	for i in range(len(x)-2):
# 		f[i] = x[0]*((data[i][0] - x[1])/float(x[1]))*(1.0 - data[i][1]) + data[i][1]*(x[2]*(data[i][0] - x[3])/float(x[3])) - x[4] - x[5]*(data[i][0] - T_eutectoid)/float(T_eutectoid) +(1.0 - data[i][1])*data[i][0]*math.log(1.0 - data[i][2]) + data[i][1]*data[i][0]*math.log(data[i][2])
# 	print i
# 	f[i+1] = x[0]*(data[i+1][0] - x[1])/float(x[1]) - data[i+1][0]*math.log((1.0 - data[i+1][1])/float((1.0 - data[i+1][2])))
# 	f[i+2] = x[2]*(data[i+1][0] - x[3])/float(x[3]) - data[i+1][0]*math.log((data[i+1][1])/float((data[i+1][2])))
# 
# 	return f
# 
# def alt_equations(p,data):
# 	x = p
# 	f = np.ones(len(x))
# 	for i in range(len(data)):	
# 		f[i] = (x[0]*(data[i][0] - x[1])/float(x[1])) - data[i][0]*math.log((1.0 - data[i][1])/float(1.0 - data[i][2]))
# 	return f
# 
# 	
# def Solve_function(flag,data):
# 	a = []
# 	print data
# 	if flag==1 : 
# 		init_value = np.ones(len(data)+1)*1.0
# 		a = fsolve(equations,init_value,args = data)
# 		return a	
# 	elif flag ==2:
# 		init_value = np.ones(len(data))*1.0
# 	#	print data
# 		b = [] 
# 		a = fsolve(alt_equations,init_value,args = data,xtol = 1.5e-5)
# 		data[:,-2:] =  1.0 - data[:,-2:]
# 		b = fsolve(alt_equations,init_value,args = data,xtol = 1.5e-5)
# 		return [a,b]
# 
# def pd_eq_below(p,temp):
# 	x = (p)
# 	f = (1.0 - c_0)*below[0]*(temp - below[1])/float(below[1]) + below[2]*c_0*(temp-below[3])/float(below[3]) + c_0*temp*math.log(abs(x[0]))+ temp*(1.0 - c_0)*math.log(abs(1.0 - x[0])) - below[4] - below[5]*(temp - T_eutectoid)/T_eutectoid
# 	return f
# 
# def pd_eq_above(p,temp):
# 	x = p
# 	f1 = above[0][0]*(temp - above[0][1])/float(above[0][1]) - temp*math.log(abs((1.0 - x[0])/(1.0 - x[1])))
# 	f2 = above[1][0]*(temp - above[1][1])/float(above[1][1]) - temp*math.log(abs((x[0])/(x[1])))
# 	return f1,f2
# 	 
# 
# def PD_Construction(temp):
# 	if temp<T_eutectoid:
# 		init_value = np.ones(1)*0.01
# 		c_alpha = fsolve(pd_eq_below,init_value,args = temp,xtol = 1.5e-5)
# 		return c_0,c_alpha[0]
# 	elif temp >= T_eutectoid:
# 		init_value = np.ones(2)*0.01
# 		c_alpha,c_gamma = fsolve(pd_eq_above,init_value,args = temp,xtol = 1.5e-5)
# 		return c_alpha,c_gamma
# 
# fp = open('thermo.inp','r')
# data = []
# for line in fp.readlines():
# 	data.append(line.split())
# for i in range(len(data)):
# 	data[i][0] = float(data[i][0])
# 	data[i][1] = wt2mole(float(data[i][1]))
# 	data[i][2] = wt2mole(float(data[i][2]))
# 
# c_0 = data[0][1]
# T_eutectoid = 1.0
# below = []
# above = []
# data_use = np.array(data)
# flag = 1
# below = Solve_function(flag,data_use[:-1,:])
# flag = 2
# above = Solve_function(flag,data_use[-2:,:])
# print below
# print "B_Fe,T_Fe,B_c,T_c",above
# print "Unknown Variables Calculated for fitting Phase Diagram"
# print "Phase Diagram Construction"
# sort_temp = [data[i][0] for i in range(len(data))] 
# T_steps = 1000
# #T_0 = sort_temp[0]
# #T_final = sort_temp[-1]
# T_0 = 0.7
# T_final = 1.17
# dT = T_final - T_0
# dT = dT/float(T_steps)
# delta = T_0
# print delta
# c1,c2 = 0.0,0.0
# 
# ##############Writing fitting parmeters to a file############
# fo = open('phasediagram.txt','w')
# fileout = open("pd_output.txt",'w')
# s = (str(below[0])+' '+str(below[1])+' '+str(below[2])+' '+str(below[3])+' '+str(below[4])+' '+str(below[5])+'\n')
# fileout.write(s)
# print s
# s = (str(above[0][0])+' '+str(above[0][1])+' '+str(above[1][0])+' '+str(above[1][1])+'\n')
# print s
# fileout.write(s)
# fileout.close()
# ###############################################################
# 
# for step in range(T_steps):
# 	c1,c2 = PD_Construction(delta)
# #	print delta,c1,c2
# 	if c1>0.0 and c2>0.0:
# 		s = str(delta)+' '+str(abs(c1))+' '+str(abs(c2))+'\n'
# 		fo.write(s)	
# 	delta += dT
# 
# #############PLotting G at eutectoid temperature##########
# 
# del_c = []
# for c_i in range(1,500):
# 	del_c.append(c_i/10000.0)
# G_alpha_below = [(below[0]*((1.0 - below[1])/(below[1]))*(1.0 - c) +below[2]*((1.0 - below[3])/(below[3]))*c + T_eutectoid*(c*math.log(c) + (1.0 - c)*(math.log(1.0 - c)))) for c in del_c] 
# G_alpha_above = [above[0][0]*((1.0 - above[0][1])/(above[0][1]))*(1.0 - c) +above[1][0]*((1.0 - above[1][1])/(above[1][1]))*c + T_eutectoid*(c*math.log(c) + (1.0 - c)*(math.log(1.0 - c))) for c in del_c]
# G_pearlite = below[4]
# G_gamma = [T_eutectoid*(c*math.log(c) +(1.0-c)*(math.log(1.0-c))) for c in del_c]
# 
# 
# 
# fontP = FontProperties()
# fontP.set_size('small')
# #legend([plot1], "title", prop = fontP)
# pylab.figure()
# pylab.plot(del_c,G_alpha_below,label="Alpha-below eutectoid")
# pylab.plot(del_c,G_alpha_above,label="Alpha-above eutectoid")
# pylab.plot(del_c,G_gamma,label="gamma")
# pylab.plot(c_0,below[4],'or',label="pearlite")
# pylab.legend(prop = fontP)
# pylab.savefig('Free_energy_Eutectoid.png') 

def wt2mole(Bulk_Comp,Mol_Weight):
	avg_comp = []
	C = 0
  	for i in range(0,len(Mol_Weight)):
  		C += (float(Bulk_Comp[i][1])/(float(Mol_Weight[i])))
 	for i in range(0,len(Mol_Weight)):
 		s = Bulk_Comp[i][0]+' '+str(float(Bulk_Comp[i][1])/((float)(Mol_Weight[i]))/(C))
 		avg_comp.append(s)
 	return avg_comp			 


def tielineestimation(Temperature):
	ref_data = []
	with open('ternary.inp','r') as fp:
		for line in fp:
			data = line.rstrip().split()
			if line.strip():
				if float(data[0]) >= (Temperature):
					T_ref = float(data[0])
					line = fp.next()
					ncounts = line.strip().split()[0]
					print T_ref,ncounts,Temperature
					for i in range(int(ncounts)):
						line = fp.next()
						if line.strip():
							ref_data.append(line.strip().split())
					break

	flag = 0	
	for num in range(len(ref_data)):
		tie_line_conf = []
		#print float(ref_data[num][1]),float(Bulk_Comp[0].split()[1]),float(ref_data[num][1]) - float(Bulk_Comp[0].split()[1])
		if float(ref_data[num][1]) - float(Bulk_Comp[0].split()[1]) >= 0.0:
			f_C = (float(Bulk_Comp[0].split()[1]) - float(ref_data[num][1]))/(float(ref_data[num][0]) - float(ref_data[num][1]))
			f_Mn =(float(Bulk_Comp[1].split()[1]) - float(ref_data[num][3]))/(float(ref_data[num][2]) - float(ref_data[num][3]))
			print ref_data[num]
			print f_C,f_Mn,"PLE"
			if(f_C > 0 and f_Mn>0 and f_Mn<=1):
				tie_line_conf.append(ref_data[num])  
				flag = 1
				break
			break
		if num == len(ref_data)-1 :
			flag = 1
			print "Still in austenite Phase at holding temperature"
			break

	for num in range(0,len(ref_data)):
		if float(ref_data[num][2]) - float(Bulk_Comp[1].split()[1]) <= 0.0 and flag == 0:
			f_C = (float(Bulk_Comp[0].split()[1]) - float(ref_data[num][1]))/(float(ref_data[num][0]) - float(ref_data[num][1]))
                	f_Mn =(float(Bulk_Comp[1].split()[1]) - float(ref_data[num][3]))/(float(ref_data[num][2]) - float(ref_data[num][3]))
			print ref_data[num]
                	print f_C,f_Mn,"NPLE"
			if(f_C > 0 and f_Mn>0 and f_Mn<=1):
                        	tie_line_conf.append(ref_data[num])
                        	flag = 1
                        	break
			break
	fp.close()
	return tie_line_conf


def equation_solve(p,c_data):
	x = p
	#print p
	f = np.ones(len(x))
	for i in range(0,len(c_data)):
		f[i] = x[0]*((c_data[i][2] - x[1])/(x[1]))  - c_data[i][2]*math.log(c_data[i][0]/c_data[i][1])
	#print f,p
	return f	



inp_data = []
for line in open('input.txt'):
	li = line.strip()
	if not li.startswith("#"):
		inp_data.append(li)
#print inp_data

Element = []
Mol_Weight = []
Bulk_Comp = []

current_line = 0

N_Components = (int)(inp_data[current_line])
current_line += 1
N_Phases = (int)(inp_data[current_line])
current_line += 1
for phase in range(0,N_Components):
	Element.append(inp_data[current_line].split(' ')[0])
	Mol_Weight.append(inp_data[current_line].split(' ')[1])
	current_line += 1

for phase in range(0,N_Components):
	Bulk_Comp.append(inp_data[current_line].split(' '))
	current_line += 1


##################Change into mole fraction####################
##############################################################
Bulk_Comp = wt2mole(Bulk_Comp,Mol_Weight)
print Bulk_Comp	

##############################################################
Heating_Rate,T_Initial,T_Hold = float(inp_data[current_line]), float(inp_data[current_line+1]), float(inp_data[current_line+2])
current_line += 3

if T_Hold == T_Initial:
	T_second = T_Initial - 20.0
else :
	T_second = T_Initial
#print T_Hold,T_Initial,Heating_Rate

###############################################################
####Read converted ternary data file#######################
tie_conf_1 = tielineestimation(T_Hold)
tie_conf_1.append(T_Hold)
tie_conf_2 = tielineestimation(T_second)
tie_conf_2.append(T_second)
############################################################################### 	 
print "Final Tie lines selected for two temperatures\n",tie_conf_1,tie_conf_2

###############################################################################
######Free Energy Construction###############################################

#For B_Fe,T_Fe#####
c_data = []
a,b = [],[]
a.append((1.0 - float(tie_conf_1[0][1])-float(tie_conf_1[0][3]))) 
a.append((1.0 - float(tie_conf_1[0][0])-float(tie_conf_1[0][2])))
a.append(float(tie_conf_1[1]))
b.append((1.0 - float(tie_conf_2[0][1])-float(tie_conf_2[0][3]))) 
b.append((1.0 - float(tie_conf_2[0][0])-float(tie_conf_2[0][2])))
b.append(float(tie_conf_2[1]))
c_data.append(a)
c_data.append(b)
print "Fe\n"
print "X_gamma \t X_alpha \t Temp\n",c_data
init_value = np.ones(2)*1.0
c_Fe = fsolve(equation_solve,init_value,args = c_data) 
print c_Fe

#For B_C,T_C

c_data = []
a,b = [],[]
a.append((float(tie_conf_1[0][1]))) 
a.append(float(tie_conf_1[0][0]))
a.append(float(tie_conf_1[1]))
b.append(float(tie_conf_2[0][1])) 
b.append(float(tie_conf_2[0][0]))
b.append(float(tie_conf_2[1]))
c_data.append(a)
c_data.append(b)
print "C\n"
print "X_gamma \t X_alpha \t Temp\n",c_data
init_value = np.array([-8000.0,1700])
c_C = fsolve(equation_solve,init_value,args = c_data,xtol = 1.5e-5) 
print c_C

#For B_Mn,T_Mn
c_data = []
a,b = [],[]
a.append(float(tie_conf_1[0][3])) 
a.append(float(tie_conf_1[0][2]))
a.append(float(tie_conf_1[1]))
b.append(float(tie_conf_2[0][3])) 
b.append(float(tie_conf_2[0][2]))
b.append(float(tie_conf_2[1]))
c_data.append(a)
c_data.append(b)
print "Mn\n"
print "X_gamma \t X_alpha \t Temp\n",c_data
init_value = np.array([-4000.0,1300])
c_Mn = fsolve(equation_solve,init_value,args = c_data) 
print c_Mn
###########################################################################


with open('constants.inp','w') as f:
	s = str(c_C[0])+' '+str(c_C[1])+'\n'
	f.write(s)
	s = str(c_Mn[0])+' '+str(c_Mn[1])+'\n'
	f.write(s)
	s = str(c_Fe[0])+' '+str(c_Fe[1])+'\n'
	f.write(s)
f.close()

with open('init_comp.inp','w') as fc:
	tie_line_init = tielineestimation(T_Initial)
	print tie_line_init
	s = str(tie_line_init[0][0])+' '+str(tie_line_init[0][1])+' '+str(tie_line_init[0][2])+' '+str(tie_line_init[0][3])+'\n'
	fc.write(s)
	tie_line_hold = tielineestimation(T_Hold)
	print tie_line_hold
	s = str(tie_line_hold[0][0])+' '+str(tie_line_hold[0][1])+' '+str(tie_line_hold[0][2])+' '+str(tie_line_hold[0][3])+'\n'
	fc.write(s)
fc.close()
	





