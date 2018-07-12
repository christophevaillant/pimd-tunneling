#!/usr/bin/python

# GDML Force Field
# Author: Stefan Chmiela (stefan@chmiela.com)
VERSION = 180118

import sys
import scipy.spatial.distance
import numpy as np
cimport numpy as np

np.import_array()

cdef np.float64_t[:,:] trainPts_D, r_d_desc_alpha, R_desc_perms, R_d_desc_alpha_perms
cdef np.int64_t[:] perms_lin
cdef np.int64_t[:,:] d_desc_mask_view
cdef double c, sig
cdef int n_train, n_atoms

cdef public GDML_load():
	global trainPts_D,r_d_desc_alpha,R_desc_perms,R_d_desc_alpha_perms,perms_lin,sig,n_perms,b_size, d_desc_mask,c,n_atoms

	model_filename = 'aspirin_ccsd.npz'

	# Load model.
	try:
		model = np.load(model_filename)
	except:
		sys.exit("ERROR: Reading model file failed.")

	trainPts_D = model['trainPts_D'].T
	r_d_desc_alpha = np.squeeze(model['r_d_desc_alpha'])

	sig = model['sig']

	c = model['c'] # self
	n_train = trainPts_D.shape[0] # self

	# Precompute batch permutations.
	if 'train_perms' in model:
		perms = model['train_perms']
	else:
		perms = np.arange(0,trainPts_D.shape[1]).reshape(1,-1)
	perm_offsets = np.arange(perms.shape[0])[:,None] * perms.shape[1]
	perms_lin = (perms + perm_offsets).ravel(order='F')

	# Precompute everything permutated training descriptors and its first derivatives multiplied with the coefficients (only needed for cached variant).
	n_perms = perm_offsets.shape[0]
	R_desc_perms = np.reshape(np.tile(trainPts_D, n_perms)[:,perms_lin], (n_train*n_perms,-1), order='F')
	R_d_desc_alpha_perms = np.reshape(np.tile(r_d_desc_alpha, n_perms)[:,perms_lin], (n_train*n_perms,-1), order='F')

	# Precompute indices for nonzero entries in desriptor derivatices.
	n_atoms = model['z'].shape[0]
	d_desc_mask = np.zeros((n_atoms,n_atoms-1), dtype=np.int64) # self
	for a in range(0,n_atoms): # for each partial deriavative
		rows,cols = np.tril_indices(n_atoms,-1)
		d_desc_mask[a,:] = np.concatenate([np.where( rows == a)[0], np.where( cols == a)[0]])

	d_desc_mask_view = d_desc_mask
	
cdef r_to_desc(r,pdist):
	n_atoms = r.shape[0]
	return 1 / pdist[np.tril_indices(n_atoms,-1)]

cdef r_to_d_desc(r,pdist):

	global d_desc_mask

	n_atoms = r.shape[0]
	d_dim = (n_atoms**2 - n_atoms)/2

	np.seterr(divide='ignore', invalid='ignore') # ignore division by zero below
	grad = np.zeros((d_dim,3*n_atoms))
	for a in range(0,n_atoms):

		d_dist = (r - r[a,:]) / (pdist[a,:]**3)[:,None]

		idx = d_desc_mask[a,:]
		grad[idx,(3*a):(3*a+3)] = np.delete(d_dist, a, axis=0)

	return grad
		
		
cdef public gdml_predict_py(double *r_c, double *grad, double *energy):

	global R_desc_perms,R_d_desc_alpha_perms,sig,n_perms,c
	
	r = np.empty(n_atoms*3)
	cdef Py_ssize_t i
	for i in range(n_atoms*3):
		r[i] = r_c[i]
	r = r.reshape(-1,3)
	
	pdist = scipy.spatial.distance.pdist(r,'euclidean')
	pdist = scipy.spatial.distance.squareform(pdist)
	
	# Generate input descriptor and its gradient.
	r_desc = r_to_desc(r,pdist)
	r_d_desc = r_to_d_desc(r,pdist)
	
	b_size = 1000

	cdef double [:] grad_mv = <double[:(r_d_desc.shape[1])]> grad
	cdef double* energy_mv = <double*> energy
	
	wkr_start = 0
	wkr_stop = 1000

	mat52_base_fact = 5/(3*sig**3)
	diag_scale_fact = 5/sig
	sqrt5 = np.sqrt(5)

	# Predict forces and energy.
	E = 0.0
	F = np.zeros((r_d_desc.shape[1],))

	wkr_start *= n_perms
	wkr_stop *= n_perms
	b_start = wkr_start
	for b_stop in range(wkr_start+b_size*n_perms,wkr_stop,b_size*n_perms) + [wkr_stop]:

		rj_desc_perms = np.asarray(R_desc_perms)[b_start:b_stop,:]
		rj_d_desc_alpha_perms = np.asarray(R_d_desc_alpha_perms)[b_start:b_stop,:]
		
		diff_ab_perms = r_desc - rj_desc_perms # b_len*6,36
		norm_ab_perms = sqrt5 * np.linalg.norm(diff_ab_perms, axis=1) # b_len*6,
		
		mat52_base = np.exp(-norm_ab_perms / sig) * mat52_base_fact # b_len*6,
		a_x2 = np.einsum('ji,ji->j', diff_ab_perms, rj_d_desc_alpha_perms) # b_len*6,

		F += np.einsum('ji,j->i', diff_ab_perms.dot(r_d_desc), a_x2 * mat52_base) * diag_scale_fact # 27,
		mat52_base *= (norm_ab_perms + sig) # b_len*6,

		F += np.einsum('ji,j', rj_d_desc_alpha_perms.dot(r_d_desc), -mat52_base)
		E += np.sum(a_x2 * mat52_base)
		
		b_start = b_stop
	
	energy_mv[0] = E + c
	for i from 0 <= i < 3*n_atoms:
		grad_mv[i] = -F[i]
	
	#cdef Py_ssize_t a,e
	#for a from 0 <= a < n_atoms:
	#	for e from 0 <= e < 3:
	#		grad_mv[a*3 + e] = -F[a + e*n_atoms]
