import sys, math
from copy import deepcopy

import numpy as np
from numpy import linalg

from DPpack.SetGlobals import *


epsilon = 1e-8

#######################################  functions  ######################################

def best_previous_point():
	
	min_energy = 0
	idx = 0
	for energy in internal['energy'][:-1]:
		if energy < min_energy or abs(energy - min_energy) < 1e-10:
			min_energy = energy
			min_idx = idx
		idx += 1
	
	return min_idx



def best_point():
	
	min_energy = 0
	idx = 0
	for energy in internal['energy']:
		if energy < min_energy or abs(energy - min_energy) < 1e-10:
			min_energy = energy
			min_idx = idx
		idx += 1
	
	return min_idx



def line_search(fh):
	
	X1 = internal['position'][-1]		# numpy array
	e1 = internal['energy'][-1]
	G1 = internal['gradient'][-1]		# numpy array
	
	idx = best_previous_point()
	X0 = internal['position'][idx]		# numpy array
	e0 = internal['energy'][idx]
	G0 = internal['gradient'][idx]		# numpy array
	
	# First try a quartic fit
	fh.write("Attempting a quartic fit.\n")
	success, y0 = quartic_fit(X0, X1, e0, e1, G0, G1, fh)
	if success and y0 > 0:
		if y0 < 1:
			new_point = X0 + y0*(X1 - X0)
			new_gradient = interpolate_gradient(G0, G1, y0)
			new_gradient = perpendicular_projection(new_gradient, X1 - X0)
			fh.write("Line search succeded.\n")
			return True, new_point, new_gradient
		else:
			idx = best_point()
			if idx == len(internal['energy']) - 1:
				new_point = X0 + y0*(X1 - X0)
				new_gradient = interpolate_gradient(G0, G1, y0)
				new_gradient = perpendicular_projection(new_gradient, X1 - X0)
				fh.write("Line search succeded.\n")
				return True, new_point, new_gradient
			else:
				fh.write("Quartic step is not acceptable. ")
	elif success:
		fh.write("Quartic step is not acceptable. ")
	
	# If no condition is met, then y0 is unacceptable. Try the cubic fit next
	fh.write("Attempting a cubic fit.\n")
	success, y0 = cubic_fit(X0, X1, e0, e1, G0, G1, fh)
	if success and y0 > 0:
		if y0 < 1:
			new_point = X0 + y0*(X1 - X0)
			new_gradient = interpolate_gradient(G0, G1, y0)
			new_gradient = perpendicular_projection(new_gradient, X1 - X0)
			fh.write("Line search succeded.\n")
			return True, new_point, new_gradient
		else:
			previous_step = X1 - internal['position'][-2]
			previous_step_size = linalg.norm(previous_step)
			new_point = X0 + y0*(X1 - X0)
			step = new_point - X1
			step_size = linalg.norm(step)
			if step_size < previous_step_size:
				new_gradient = interpolate_gradient(G0, G1, y0)
				new_gradient = perpendicular_projection(new_gradient, X1 - X0)
				fh.write("Line search succeded.\n")
				return True, new_point, new_gradient
			else:
				fh.write("Cubic step is not acceptable. ")
	elif success:
		fh.write("Cubic step is not acceptable. ")
	
	# If no condition is met again, then all fits fail.
	fh.write("All fits fail. ")
	
	# Then, if the latest point is not the best, use y0 = 0.5 (step to the midpoint)
	idx = best_point()
	if idx < len(internal['energy']) - 1:
		y0 = 0.5
		new_point = X0 + y0*(X1 - X0)
		new_gradient = interpolate_gradient(G0, G1, y0)
		new_gradient = perpendicular_projection(new_gradient, X1 - X0)
		fh.write("Moving to the midpoint.\n")
		return True, new_point, new_gradient
	
	# If the latest point is the best point, no linear search is done
	fh.write("No linear search will be used in this step.\n")
	
	return False, None, None



## For cubic and quartic fits, G0 and G1 are the gradient vectors

def cubic_fit(X0, X1, e0, e1, G0, G1, fh):
		
	line = X1 - X0
	line /= linalg.norm(line)
	
	g0 = np.dot(G0, line)
	g1 = np.dot(G1, line)
	
	De = e1 - e0
	
	fh.write("De = {:<18.15e}      g0 = {:<12.8f}      g1 = {:<12.8f}\n".format(De, g0, g1))
	
	alpha = g1 + g0 - 2*De
	if abs(alpha) < epsilon:
		fh.write("Cubic fit failed: alpha too small\n")
		return False, None
	
	beta = 3*De - 2*g0 - g1
	discriminant = 4 * (beta**2 - 3*alpha*g0)
	if discriminant < 0:
		fh.write("Cubic fit failed: no minimum found (negative Delta)\n")
		return False, None
	if abs(discriminant) < epsilon:
		fh.write("Cubic fit failed: no minimum found (null Delta)\n")
		return False, None
	
	y0 = (-beta + math.sqrt(discriminant/4)) / (3*alpha)
	fh.write("Minimum found with y0 = {:<8.4f}\n".format(y0))
	
	return True, y0



def quartic_fit(X0, X1, e0, e1, G0, G1, fh):
		
	line = X1 - X0
	line /= linalg.norm(line)
	
	g0 = np.dot(G0, line)
	g1 = np.dot(G1, line)
	
	De = e1 - e0
	Dg = g1 - g0
	
	fh.write("De = {:<18.15e}      g0 = {:<12.8f}      g1 = {:<12.8f}\n".format(De, g0, g1))
	
	if Dg < 0 or De - g0 < 0:
		fh.write("Quartic fit failed: negative alpha\n")
		return False, None
	if abs(Dg) < epsilon or abs(De - g0) < epsilon:
		fh.write("Quartic fit failed: alpha too small\n")
		return False, None
	
	discriminant = 16 * (Dg**2 - 3*(g1 + g0 - 2*De)**2)
	if discriminant < 0:
		fh.write("Quartic fit failed: no minimum found (negative Delta)\n")
		return False, None
	
	alpha1 = (Dg + math.sqrt(discriminant/16)) / 2
	alpha2 = (Dg - math.sqrt(discriminant/16)) / 2
	
	fh.write("alpha1 = {:<7.4e}      alpha2 = {:<7.4e}\n".format(alpha1, alpha2))
	
	alpha = alpha1
	beta = g1 + g0 - 2*De - 2*alpha
	gamma = De - g0 - alpha - beta
	
	y0 = (-1/(2*alpha)) * ((beta**3 - 4*alpha*beta*gamma + 8*g0*alpha**2)/4)**(1/3)
	fh.write("Minimum found with y0 = {:<8.4f}\n".format(y0))
	
	return True, y0



def rfo_step(gradient, hessian, type):
	
	dim = len(gradient)
	
	aug_hessian = []
	for i in range(dim):
		aug_hessian.extend(hessian[i,:].tolist())
		aug_hessian.append(gradient[i])
	
	aug_hessian.extend(gradient.tolist())
	aug_hessian.append(0)
	
	aug_hessian = np.array(aug_hessian).reshape(dim + 1, dim + 1)
	
	evals, evecs = linalg.eigh(aug_hessian)
	
	if type == "min":
		step = np.array(evecs[:-1,0])
	elif type == "ts":
		step = np.array(evecs[:-1,1])
	
	return step



def update_trust_radius():
	
	if internal['trust_radius'] == None:
		internal['trust_radius'] = player['maxstep']
	elif len(internal['energy']) > 1:
		X1 = internal['position'][-1]
		X0 = internal['position'][-2]
		Dx = X1 - X0
		displace = linalg.norm(Dx)
		e1 = internal['energy'][-1]
		e0 = internal['energy'][-2]
		De = e1 - e0
		g0 = internal['gradient'][-2]
		h0 = internal['hessian'][-2]
		
		rho = De / (np.dot(g0, Dx) + 0.5*np.dot(Dx, np.matmul(h0, Dx.T).T))
		
		if rho > 0.75 and displace > 0.8*internal['trust_radius']:
			internal['trust_radius'] = 2*internal['trust_radius']
		elif rho < 0.25:
			internal['trust_radius'] = 0.25*displace
	
	return



def interpolate_gradient(G0, G1, y0):
	
	DG = G1 - G0
	gradient = G0 + y0*DG
	
	return gradient



def perpendicular_projection(vector, line):
	
	direction = line / linalg.norm(line)
	projection = np.dot(vector, direction) * direction
	
	return vector - projection



