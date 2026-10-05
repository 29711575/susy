#!/usr/bin/env python3
"""Deterministic tests for off-grid finite-window structure factors."""
import unittest
import numpy as np

from vbs_window import (configuration_spatially_centered_F,
    connected_window_association_ratios, ensemble_connected_matrix,
    joint_phase_ratio_window_jackknife, largest_hermitian_eigenvalue,
    reciprocal_geometry, square_window_weights,
    translation_averaged_connected_matrix, window_factor, windowed_matrix)

Q = np.array([0.5804306232551663, 0.4195693767448336])


def translation_orbit_metrics(L=6, bond_length=1.0):
    nc = L * L
    xx, yy = np.meshgrid(np.arange(L), np.arange(L), indexing='ij')
    phase1 = 2*np.pi*(xx/L + 2*yy/L)
    phase2 = 2*np.pi*(2*xx/L + yy/L)
    c_values = np.array([[-.4,.1,.2],[.0,.3,-.1],[.5,-.2,.4]])
    avec = np.array([[.8,.2,-.3],[-.1,.5,.7]])
    bvec = np.array([[.1,-.6,.4],[.3,.2,.5]])
    base_fields = []
    for c in c_values:
        field = c[None,None,:] + np.cos(phase1)[...,None]*avec[0] + np.sin(phase1)[...,None]*avec[1]
        field += np.cos(phase2)[...,None]*bvec[0] + np.sin(phase2)[...,None]*bvec[1]
        base_fields.append(field)

    weights, _ = square_window_weights(L, Q, bond_length=bond_length)
    a1, a2, b1, b2 = reciprocal_geometry(bond_length=bond_length)
    qcart = Q[0]*b1 + Q[1]*b2
    rcart = xx[...,None]*a1 + yy[...,None]*a2
    qphase = np.exp(-1j*np.einsum('...i,i->...',rcart,qcart))
    raw_direct=[]; projected=[]; grid_matrices=[]; fqs=[]; fzeros=[]
    covariance_by_origin=np.zeros((L,L,L,L,3,3),dtype=float)
    ensemble_fields=[]
    for base in base_fields:
        for tx in range(L):
            for ty in range(L):
                field=np.roll(base,shift=(-tx,-ty),axis=(0,1))
                ensemble_fields.append(field)
                fq=np.einsum('xy,xyi->i',qphase,field)
                f0=field.sum(axis=(0,1)).astype(complex)
                fqs.append(fq);fzeros.append(f0)
                raw_direct.append(np.outer(fq,fq.conj())/nc)
                proj=configuration_spatially_centered_F(
                    fq,f0,window_factor(L,Q,bond_length=bond_length))
                projected.append(np.outer(proj,proj.conj())/nc)
                fk=np.empty((L,L,3),complex)
                for n1 in range(L):
                    for n2 in range(L):
                        ph=np.exp(-2j*np.pi*(n1*xx/L+n2*yy/L))
                        fk[n1,n2]=np.einsum('xy,xyi->i',ph,field)
                grid_matrices.append(fk[..., :,None]*fk[...,None,:].conj()/nc)
    raw_direct=np.mean(raw_direct,axis=0)
    projected=np.mean(projected,axis=0)
    grid_matrices=np.mean(grid_matrices,axis=0)
    conv=windowed_matrix(grid_matrices,weights)
    ensemble_fields=np.asarray(ensemble_fields)
    mean_field=ensemble_fields.mean(axis=0)
    for dx in range(L):
        for dy in range(L):
            shifted=np.roll(ensemble_fields,shift=(-dx,-dy),axis=(1,2))
            covariance_by_origin[dx,dy]=np.einsum('bxyi,bxyj->xyij',
                ensemble_fields-mean_field,shifted-mean_field)/len(ensemble_fields)
    cov_residual=float(np.max(np.ptp(covariance_by_origin,axis=(2,3))))
    mean_f=np.mean(fqs,axis=0)
    connected=ensemble_connected_matrix(raw_direct,mean_f,nc)
    translated_connected=translation_averaged_connected_matrix(
        conv,np.mean(fzeros,axis=0)/nc,window_factor(L,Q,bond_length=bond_length),nc)
    zero_contribution=float(weights[0,0])*grid_matrices[0,0]
    single=ensemble_fields[0]
    fk=np.empty((L,L,3),complex)
    for n1 in range(L):
        for n2 in range(L):
            ph=np.exp(-2j*np.pi*(n1*xx/L+n2*yy/L))
            fk[n1,n2]=np.einsum('xy,xyi->i',ph,single)
    single_grid=fk[..., :,None]*fk[...,None,:].conj()/nc
    single_q=np.einsum('xy,xyi->i',qphase,single)
    single_direct=np.outer(single_q,single_q.conj())/nc
    single_conv=windowed_matrix(single_grid,weights)
    return {"raw_direct_mean":raw_direct,"window_convolution":conv,
            "ensemble_connected":connected,"configuration_spatially_centered":projected,
            "translation_averaged_connected":translated_connected,
            "mean_FQ":mean_f,"mean_Fzero":np.mean(fzeros,axis=0),
            "zero_mode_contribution":zero_contribution,
            "matrix_identity_max_abs":float(np.max(np.abs(raw_direct-conv))),
            "covariance_translation_residual":cov_residual,
            "single_configuration_cross_term_residual":float(np.max(np.abs(single_direct-single_conv))),
            "single_to_windowed_matrix":single_conv}


class WindowStructureFactorTests(unittest.TestCase):
    def test_cartesian_and_fractional_phase_match_once(self):
        L=9; bond_length=1.0; p=np.array([.37,-.21])
        direct=window_factor(L,p,bond_length=bond_length)
        frac=sum(np.exp(-2j*np.pi*(p[0]*x+p[1]*y))
                 for x in range(L) for y in range(L))/(L*L)
        self.assertLess(abs(direct-frac),2e-14)
        a1,a2,b1,b2=reciprocal_geometry(bond_length=bond_length)
        np.testing.assert_allclose([np.dot(a1,b1),np.dot(a1,b2),np.dot(a2,b1),np.dot(a2,b2)],
                                   [2*np.pi,0,0,2*np.pi],atol=2e-14,rtol=0)
        x,y=3,5
        self.assertLess(abs(np.exp(-1j*np.dot(p[0]*b1+p[1]*b2,x*a1+y*a2))-
                            np.exp(-2j*np.pi*(p[0]*x+p[1]*y))),2e-14)

    def test_weight_sum_zero_bin_and_dominant_momenta(self):
        for L in (6,9,12):
            weights,rows=square_window_weights(L,Q,bond_length=1.0)
            self.assertAlmostEqual(float(weights.sum()),1.,places=13)
            self.assertAlmostEqual(float(weights[0,0]),
                abs(window_factor(L,Q,bond_length=1.0))**2,places=14)
            self.assertEqual(len(rows),L*L)
            self.assertGreater(max(row['weight'] for row in rows),
                               min(row['weight'] for row in rows))

    def test_translation_invariant_known_covariance_and_all_shifts(self):
        m=translation_orbit_metrics()
        self.assertLess(m['matrix_identity_max_abs'],2e-12)
        np.testing.assert_allclose(m['ensemble_connected'],
                                   m['translation_averaged_connected'],atol=2e-12,rtol=0)
        self.assertLess(m['covariance_translation_residual'],2e-13)
        self.assertGreater(m['single_configuration_cross_term_residual'],1e-5)
        self.assertGreater(np.linalg.norm(m['ensemble_connected']-
                                          m['configuration_spatially_centered']),1e-5)

    def test_matrix_is_weighted_before_largest_eigenvalue(self):
        weights=np.array([[.5,.5]])
        s=np.zeros((1,2,3,3),complex)
        s[0,0]=np.diag([4.,0.,0.]); s[0,1]=np.diag([0.,4.,0.])
        matrix=windowed_matrix(s,weights)
        self.assertAlmostEqual(largest_hermitian_eigenvalue(matrix),2.)
        self.assertAlmostEqual(sum(w*largest_hermitian_eigenvalue(x)
                                   for w,x in zip(weights.ravel(),s.reshape(-1,3,3))),4.)

    def test_joint_momentum_covariance_survives_window_jackknife(self):
        weights=np.array([[.5,.5]])
        phase=np.ones(4,complex)
        numer=np.zeros((4,1,2,3,3),complex)
        for b in range(4):
            numer[b,0,0]=np.diag([5+b,1,1+b])
            numer[b,0,1]=np.diag([1,7-b,2+2*b])
        out=joint_phase_ratio_window_jackknife(phase,numer,weights)
        cov=out['szz_real_jackknife_covariance']
        self.assertGreater(cov[0,1],0.)
        self.assertAlmostEqual(out['matrix'][0,0].real,
                               np.mean([(.5*(5+b+1)) for b in range(4)]))
        replicas=out['replica_leading_eigenvalues']
        expected=np.sqrt((len(replicas)-1)/len(replicas)*
                         np.sum((replicas-replicas.mean())**2))
        self.assertAlmostEqual(out['leading_eigenvalue_jackknife_error'],expected,places=14)

    def test_connected_window_jackknife_recomputes_first_moments(self):
        weights=np.array([[.25,.75]])
        amplitudes=np.array([[.5+0j,.5+0j]])
        phase=np.array([1.,.8,.9,1.1],complex)
        s=np.zeros((4,1,2,3,3),complex)
        d=np.zeros((4,1,2,3),complex)
        for b in range(4):
            s[b,0,0]=np.diag([1.,2.,3.+b])
            s[b,0,1]=np.diag([2.,1.,4.+2*b])
            d[b,0,0,2]=.03*(b+1)
            d[b,0,1,2]=-.02*b
        out=joint_phase_ratio_window_jackknife(phase,s,weights,d,amplitudes,Ncell=2)
        raw=np.einsum('xy,xyij->ij',weights,s.sum(axis=0)/phase.sum())
        mean_d=np.einsum('xy,xyi->i',amplitudes,d.sum(axis=0)/phase.sum())
        np.testing.assert_allclose(out['connected_matrix'],raw-2*np.outer(mean_d,mean_d.conj()))
        self.assertEqual(out['replica_connected_matrices'].shape,(4,3,3))
        self.assertEqual(out['connected_szz_real_jackknife_covariance'].shape,(2,2))
        for b in range(4):
            den=phase.sum()-phase[b]
            mean_db=np.einsum('xy,xyi->i',amplitudes,(d.sum(axis=0)-d[b])/den)
            rawb=np.einsum('xy,xyij->ij',weights,(s.sum(axis=0)-s[b])/den)
            np.testing.assert_allclose(out['replica_connected_matrices'][b],
                                       rawb-2*np.outer(mean_db,mean_db.conj()))

    def test_white_noise_covariance_and_projected_zero_mode_factor(self):
        # A balanced finite ensemble of signed site/channel impulses has
        # E[X_mu(r)X_nu(r')]=delta_mu,nu delta_r,r' exactly.
        L=4; nc=L*L; amp=np.sqrt(3.*nc)
        xx,yy=np.meshgrid(np.arange(L),np.arange(L),indexing='ij')
        Qs=[Q,Q+np.array([1/L,0]),Q-np.array([1/L,0]),
            Q+np.array([0,1/L]),Q-np.array([0,1/L])]
        fields=[]
        for x in range(L):
            for y in range(L):
                for mu in range(3):
                    for sign in (-1.,1.):
                        f=np.zeros((L,L,3),float);f[x,y,mu]=sign*amp
                        fields.append(f)
        fields=np.asarray(fields)
        for q in Qs:
            qphase=np.exp(-2j*np.pi*(q[0]*xx+q[1]*yy))
            w=window_factor(L,q)
            direct=[];projected=[];fqs=[]
            for f in fields:
                fq=np.einsum('xy,xyi->i',qphase,f)
                f0=f.sum(axis=(0,1)).astype(complex)
                fqs.append(fq)
                direct.append(np.outer(fq,fq.conj())/nc)
                fp=configuration_spatially_centered_F(fq,f0,w)
                projected.append(np.outer(fp,fp.conj())/nc)
            raw=np.mean(direct,axis=0)
            ensemble=ensemble_connected_matrix(raw,np.mean(fqs,axis=0),nc)
            sprj=np.mean(projected,axis=0)
            np.testing.assert_allclose(raw,np.eye(3),atol=2e-14)
            np.testing.assert_allclose(ensemble,np.eye(3),atol=2e-14)
            np.testing.assert_allclose(sprj,(1.-abs(w)**2)*np.eye(3),atol=2e-14)
            gridmat=np.empty((L,L,3,3),complex)
            for n1 in range(L):
                for n2 in range(L):
                    k=np.array([n1/L,n2/L])
                    kp=np.exp(-2j*np.pi*(k[0]*xx+k[1]*yy))
                    matrices=[]
                    for f in fields:
                        fk=np.einsum('xy,xyi->i',kp,f)
                        matrices.append(np.outer(fk,fk.conj())/nc)
                    gridmat[n1,n2]=np.mean(matrices,axis=0)
            weights,_=square_window_weights(L,q)
            np.testing.assert_allclose(windowed_matrix(gridmat,weights),np.eye(3),atol=2e-14)

    def test_window_association_keeps_negative_values_and_block_floor(self):
        phase=np.ones(8,complex)
        s=np.zeros((8,5,3,3),complex)
        d=np.zeros((8,5,3),complex)
        for b in range(8):s[b,:,2,2]=[1.,4.,4.,3.,3.]
        ids={'center':0,'plus_b1':1,'minus_b1':2,'plus_b2':3,'minus_b2':4}
        result=connected_window_association_ratios(phase,s,d,ids,L=6)
        self.assertTrue(result['ratios']['b1']['available'])
        self.assertAlmostEqual(result['ratios']['b1']['estimate'][0],-3.)
        self.assertLess(result['ratios']['b1']['estimate'][0],0.)
        self.assertEqual(result['connected_szz_real_jackknife_covariance_across_selected_q'].shape,(5,5))
        short=connected_window_association_ratios(phase[:2],s[:2],d[:2],ids,L=6)
        self.assertFalse(short['ratios']['b1']['available'])
        self.assertEqual(short['ratios']['b1']['unavailable_reason'],'too_few_equal_blocks')
        unresolved=connected_window_association_ratios(phase,np.zeros_like(s),d,ids,L=6)
        self.assertFalse(unresolved['ratios']['b1']['available'])
        self.assertEqual(unresolved['ratios']['b1']['unavailable_reason'],
                         'connected_center_denominator_unresolved')
        self.assertIsNone(unresolved['ratios']['b1']['estimate'])


if __name__=='__main__':
    unittest.main(verbosity=2)
