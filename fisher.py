import pandas as pd
import pyccl as ccl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from copy import deepcopy
from math import fsum
from functools import lru_cache
from numpy.linalg import inv
import sys, scipy
import numdifftools as nd
from scipy.stats import uniform, norm
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from copy import deepcopy
from numbers import Number
from itertools import chain
import seaborn as sns
from sklearn.neighbors import KernelDensity
import pickle

def FoM(matrix):
    return np.sqrt(np.linalg.det(matrix))

def plot_contours(matrix, sigmas, fid, **kwargs):
    prefactor = {1:1.52, 2:2.48}
    prefactor = prefactor[sigmas]
    matrix = np.linalg.inv(matrix)
    s00, s01, s11 = matrix[0][0], matrix[0][1], matrix[1][1]
    a = np.sqrt(
        0.5*(s00 + s11) + np.sqrt(s01**2 + 0.25*(s00-s11)**2)
    )
    b = np.sqrt(
        0.5*(s00 + s11) - np.sqrt(s01**2 + 0.25*(s00-s11)**2)
    )
    b *= prefactor
    a *= prefactor
    theta = np.arctan(2*s01/(s00-s11))/2
    eig = np.linalg.eig(matrix)
    maxeig = eig[1][np.argmax(eig[0])]
    theta = np.arctan2(maxeig[1], maxeig[0])
    el = matplotlib.patches.Ellipse(fid, 2*a, 2*b, angle=-np.degrees(theta), alpha=0.3, **kwargs)
    return el, ((fid[0]-3*a*np.cos(theta), fid[0]+3*a*np.cos(theta)), (fid[1]-3*a*np.sin(theta), fid[1]+3*a*np.sin(theta)))

def marginalize(fisher_matrix, i, j):
    return np.linalg.inv(np.linalg.inv(fisher_matrix)[np.ix_([i,j], [i,j])]) 


def get_cov_matrix(l, data_order, cl_vals, orderings, fsky):
    prefac = 1/(2*l + 1)/fsky
    out_cov = np.zeros((len(data_order), len(data_order)))
    for z in range(out_cov.shape[0]):
        for y in range(out_cov.shape[1]):
            i, j = (data_order[z][_] for _ in [0, 1])
            k, l = (data_order[y][_] for _ in [0, 1])

            try:
                cl_ik = cl_vals[orderings.index([i, k])]
            except ValueError:
                #this happens if only autocorrelation is queried
                cl_ik = 0

            try:
                cl_jl = cl_vals[orderings.index([j, l])]
            except ValueError:
                cl_jl = 0

            try:
                cl_il = cl_vals[orderings.index([i, l])]
            except ValueError:
                cl_il = 0

            try:
                cl_jk = cl_vals[orderings.index([j, k])]
            except ValueError:
                cl_jk = 0

            out_cov[z, y] = prefac * (cl_ik*cl_jl + cl_il*cl_jk)
    return out_cov

def multi_bin_cov(fsky, Clbins, Cl_ordering, num_dens): 
        #first column --> l values

    combination_cl = Clbins[:, 0]
    orderings = []
    for i in range(1, Clbins.shape[1]):
        #Cl_orderings has no l entry
        index_comb = Cl_ordering[i-1, :].tolist()
        if index_comb[0] != index_comb[1]:
            #add the entry plus the
            #reversed entry --> would NOT be valid for cross correlations since different
            #bins!!!
            combination_cl = np.column_stack((combination_cl, Clbins[:, i]))
            combination_cl = np.column_stack((combination_cl, Clbins[:, i]))
            orderings.append(index_comb)
            orderings.append([index_comb[1], index_comb[0]])
        else:
            #autocorrelations
            combination_cl = np.column_stack((combination_cl, Clbins[:, i]))
            orderings.append(index_comb)

    #remove the first column from combination_cl --> makes it easier since now
    #it corresponds to orderings vector

    combination_cl = combination_cl[:, 1:]

    assert len(orderings) == combination_cl.shape[1]

    #add the shotnoise to each of the cl combinations:
    # TODO: make this an option for clustering vs lensing
    shotnoise = []
    for el in orderings:
        shotnoise.append(1.0/num_dens[int(el[0] - 1)])  # because ordering starts with 1
    shotnoise = np.array(shotnoise)
    assert len(shotnoise) == combination_cl.shape[1]

    for i in range(combination_cl.shape[1]):
        #only the autocorrelation is affected by shot noise
        if orderings[i][0] == orderings[i][1]:
            combination_cl[:, i] += shotnoise[i]

    #now calculate the covariance matrice for each of the Cl_orderings
    out_mat = []
    for i in range(Clbins.shape[0]):
        curr_l = Clbins[i, 0]
        matrix_out = get_cov_matrix(curr_l, Cl_ordering, combination_cl[i, :], orderings, fsky)
        out_mat.append(matrix_out)
    return out_mat

def one_bin_cov(fsky, Clbins, num_dens): 
    added_shotnoise = (Clbins[:, 1] + 1.0/num_dens)**2
    prefactor = 2.0/((2.0 * Clbins[:, 0] + 1.)*fsky) 
    covariance = prefactor * added_shotnoise 
    cov_matrix = np.diag(covariance)
    return cov_matrix
    

class PhotoZ_core(object):
    """Defines p(zp | z) of the core distribution (i.e. no outliers)
    """

    def __init__(self, zp_support, zbias, sigma_z):
        self.zp_support = zp_support
        self.zbias = zbias
        self.sigma_z = sigma_z

    def get_pdf_given(self, z: Number):
        rv = norm()
        scale = self.sigma_z * (1 + z)
        loc = z - self.zbias
        return rv.pdf((np.array(self.zp_support) - loc) / scale) / scale

class Core(object):
    """Defines the core and outlier population
    """
    def __init__(self, zp_support, zbias, sigma_z):
        self.zp_support = zp_support
        self.core = PhotoZ_core(zp_support, zbias, sigma_z)

    def get_pdf_given(self, z: Number):
        core_pdf = self.core.get_pdf_given(z)
        return core_pdf/np.trapz(core_pdf, self.zp_support)


class SmailZ(object):
    """Define the True photometric redshift distribution
    that follows a Smail Type distribution
    """
    def __init__(self, z_support, pdf_values_evaluated_at_zsupport):
        self.z_support = z_support
        self.pdf = pdf_values_evaluated_at_zsupport/np.trapz(pdf_values_evaluated_at_zsupport, z_support)

    def get_pdf(self):
        return self.pdf

    def get_pdf_convoled(self, filter_list):
        output_tomo_list = np.array([el * self.pdf for el in filter_list]).T
        output_tomo_list = np.column_stack((self.z_support, output_tomo_list))
        return output_tomo_list


class PhotozModel(object):
    """Convolve the joint distribution p(z_s, z_p) with
    a set of filter functions (e.g. gaussian or tophat)

    The class function get_pdf produces an array of tomographic
    bins

    """
    def __init__(self, pdf_z, pdf_zphot_given_z, filters):
        self.pdf_zphot_given_z = pdf_zphot_given_z
        self.pdf_z = pdf_z
        self.filter_list = filters


    def get_pdf(self):
        return np.array([self.get_pdf_tomo(el) for el in self.filter_list])

    def get_pdf_tomo(self, filter):
        z_support = self.pdf_z.z_support
        zp_support = self.pdf_zphot_given_z.zp_support
        z_pdf = self.pdf_z.get_pdf()
        pdf_joint = np.zeros((len(zp_support), len(z_support)))
        for i in range(len(z_support)):
            pdf_joint[:, i] = self.pdf_zphot_given_z.get_pdf_given(z_support[i]) * z_pdf[i] * filter[i]

        pdf_zp = np.zeros((len(zp_support),))
        for i in range(len(zp_support)):
            pdf_zp[i] = np.trapz(pdf_joint[i, :], z_support)

        return pdf_zp

def centroid_shift(unbiased, biased):
    cl_unbiased = np.array(unbiased.ccl_cls['C_ell']).reshape(15, 15).T 
    cl_biased = np.array(biased.ccl_cls['C_ell']).reshape(15, 15).T 
    
    diff_cl = cl_biased - cl_unbiased
    bias_vec = []
    for i, param in enumerate(unbiased.param_order):
        bias_vec.append(sum(diff_cl[idx].dot(
                np.array(unbiased.invcov_list[idx]).dot(unbiased.derivs_sig[param][idx])
                ) for idx in range(len(unbiased.ell))))
    bias_vec = np.array(bias_vec)
    para_bias = np.linalg.inv(fisher).dot(bias_vec) 
    para_bias = {param_order[i]: para_bias[i] for i in range(7)}
    return para_bias
    


def FullPlot(params, *args, labels='LongString'):
    colors = sns.color_palette(palette='colorblind', n_colors=len(args))
    es = []
    for i, arg in enumerate(args):
        if not arg.has_run:
            raise ValueError('This class has not been processed yet.')
        matrix = arg.fisher
        e, (xlim, ylim) = plot_contours(
            marginalize(matrix, 
                        arg.param_order.index(params[0]), 
                        arg.param_order.index(params[1])),
            sigmas=2,
            fid=(arg.vals[params[0]], arg.vals[params[1]]))
        e.set_facecolor(colors[i])
        e.set_label(labels[i])
        e.set_alpha(1/(1+len(args)))
        es.append(e)

    fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
    for i, e in enumerate(es):
        ax.add_artist(e)

    _ = plt.xlim((xlim))
    _ = plt.ylim((ylim))
    #_ = plt.legend(handles=es, loc=(0.7, 0.7))
    _ = plt.xlabel(arg.param_labels[arg.param_order.index(params[0])])
    _ = plt.ylabel(arg.param_labels[arg.param_order.index(params[1])])
    _ = plt.figure()



class Fisher:
    def __init__(self, cosmo, zvariance=[1,1,1,1,1], zbias=[0,0,0,0,0], outliers=[1,1,1,1,1], calcCov=False, plot_label=None, end=None, masked=None):
        self.zvariance = zvariance
        self.zbias = zbias
        self.outliers = outliers
        self.has_run = False
        self.A0 = 5
        self.intCache = {}
        self.beta = 1
        self.etal = 0
        self.etah = 0
        self.gbias = [1.376695, 1.451179, 1.528404, 1.607983, 1.689579, 1.772899, 1.857700, 1.943754, 2.030887, 2.118943] 
        self.IA_interp = pickle.load(open('IA_interp.p', 'rb'))
        self.calcCov = calcCov
        self.plot_label = plot_label
        self.cosmo = cosmo
        df = pd.read_csv('nzdist.txt', sep=' ') 
        self.zmid = list(df['zmid'])
        self.dneff = df['dneff']
        self.z = [0] + self.zmid
        self.funcs_ss = {
            'sigma_8': self.getC_ellOfSigma8_ss,
            'omega_b': self.getC_ellOfOmegab_ss,
            'h': self.getC_ellOfh_ss,
            'n_s': self.getC_ellOfn_s_ss,
            'omega_m': self.getC_ellOfOmegam_ss,
            'w_0': self.getC_ellOfw0_ss,
            'w_a': self.getC_ellOfwa_ss,
            'zbias1': self.getC_ellOfzbias1_ss,
            'zbias2': self.getC_ellOfzbias2_ss,
            'zbias3': self.getC_ellOfzbias3_ss,
            'zbias4': self.getC_ellOfzbias4_ss,
            'zbias5': self.getC_ellOfzbias5_ss,
            'zvariance1': self.getC_ellOfzvariance1_ss,
            'zvariance2': self.getC_ellOfzvariance2_ss,
            'zvariance3': self.getC_ellOfzvariance3_ss,
            'zvariance4': self.getC_ellOfzvariance4_ss,
            'zvariance5': self.getC_ellOfzvariance5_ss,
            'zoutlier1': self.getC_ellOfzoutlier1_ss,
            'zoutlier2': self.getC_ellOfzoutlier2_ss,
            'zoutlier3': self.getC_ellOfzoutlier3_ss,
            'zoutlier4': self.getC_ellOfzoutlier4_ss,
            'zoutlier5': self.getC_ellOfzoutlier5_ss,
            'A0': self.getC_ellOfIA_amp_ss,
            'beta': self.getC_ellOfIA_beta_ss,
            'etal': self.getC_ellOfIA_lowz_ss,
            'etah': self.getC_ellOfIA_highz_ss
        }
        self.funcs_sp = {
            'sigma_8': self.getC_ellOfSigma8_sp,
            'omega_b': self.getC_ellOfOmegab_sp,
            'h': self.getC_ellOfh_sp,
            'n_s': self.getC_ellOfn_s_sp,
            'omega_m': self.getC_ellOfOmegam_sp,
            'w_0': self.getC_ellOfw0_sp,
            'w_a': self.getC_ellOfwa_sp,
            'zbias1': self.getC_ellOfzbias1_sp,
            'zbias2': self.getC_ellOfzbias2_sp,
            'zbias3': self.getC_ellOfzbias3_sp,
            'zbias4': self.getC_ellOfzbias4_sp,
            'zbias5': self.getC_ellOfzbias5_sp,
            'zvariance1': self.getC_ellOfzvariance1_sp,
            'zvariance2': self.getC_ellOfzvariance2_sp,
            'zvariance3': self.getC_ellOfzvariance3_sp,
            'zvariance4': self.getC_ellOfzvariance4_sp,
            'zvariance5': self.getC_ellOfzvariance5_sp,
            'zoutlier1': self.getC_ellOfzoutlier1_sp,
            'zoutlier2': self.getC_ellOfzoutlier2_sp,
            'zoutlier3': self.getC_ellOfzoutlier3_sp,
            'zoutlier4': self.getC_ellOfzoutlier4_sp,
            'zoutlier5': self.getC_ellOfzoutlier5_sp,
            'A0': self.getC_ellOfIA_amp_sp,
            'beta': self.getC_ellOfIA_beta_sp,
            'etal': self.getC_ellOfIA_lowz_sp,
            'etah': self.getC_ellOfIA_highz_sp
        }
        self.funcs_pp = {
            'sigma_8': self.getC_ellOfSigma8_pp,
            'omega_b': self.getC_ellOfOmegab_pp,
            'h': self.getC_ellOfh_pp,
            'n_s': self.getC_ellOfn_s_pp,
            'omega_m': self.getC_ellOfOmegam_pp,
            'w_0': self.getC_ellOfw0_pp,
            'w_a': self.getC_ellOfwa_pp,
        }

        self.vals = {
            'sigma_8': 0.831, 
            'omega_b': 0.049, 
            'h': 0.6727, 
            'n_s': 0.9645, 
            'omega_m': 0.3156,
            'w_0': -1,
            'w_a': 0,
            'A0': 5,
            'etal':0,
            'etah': 0,
            'beta': 0
        }
        for i in range(5):
            self.vals['zbias'+str(i+1)] = 0
            self.vals['zvariance'+str(i+1)] = 1
            self.vals['zoutlier'+str(i+1)] = 1
        for i in range(1,11):
            self.funcs_pp['gbias'+str(i)] = self.getC_ellOfgbias_ss
            self.vals['gbias'+str(i)] = self.gbias[i-1]
        self.priors = {
            'sigma_8': 1/0.14**2, 
            'omega_b': 1/0.006**2, 
            'h': 1/0.063**2,
            'n_s': 1/0.08**2, 
            'omega_m': 1/0.2**2,
            'w_0': 1/0.8**2,
            'w_a': 1/2**2,
            'zbias1': 1/2**2,
            'zbias2': 1/2**2,
            'zbias3': 1/2**2,
            'zbias4': 1/2**2,
            'zbias5': 1/2**2,
            'zvariance1': 1/2**2,
            'zvariance2': 1/2**2,
            'zvariance3': 1/2**2,
            'zvariance4': 1/2**2,
            'zvariance5': 1/2**2,
            'zoutlier1': 1/2**2,
            'zoutlier2': 1/2**2,
            'zoutlier3': 1/2**2,
            'zoutlier4': 1/2**2,
            'zoutlier5': 1/2**2,
            'A0': 1/2**2,
            'beta': 1/2**2,
            'etal': 1/2**2,
            'etah': 1/2**2
        }
        for i in range(10):
            string = 'gbias'+str(i+1)
            self.priors[string] = 1/0.9**2
        self.param_order = ['omega_m', 'sigma_8', 'n_s', 'w_0', 'w_a', 'omega_b', 'h', 'A0', 'beta', 'etal', 'etah'] + ['zbias'+str(i) for i in range(1, 6)] + ['zvariance'+str(i) for i in range(1,6)] + ['zoutlier'+str(i) for i in range(1,6)] + ['gbias'+str(i) for i in range(1,11)]
        self.param_labels = [r'$\Omega_m$', r'$\sigma_8$', r'$n_s$', r'$w_0$', r'$w_a$', r'$\Omega_b$', r'$h$', r'$A_0$', r'$\beta$', r'$\eta_l$', r'$\eta_h$'] + [r'$z_{bias}$'+str(i) for i in range(1, 6)] + [r'$\std{z}$'+str(i) for i in range(1,6)] + [r'$out$'+str(i) for i in range(1,6)] + [rf'$b_g^{i}' for i in range(1,11)]
        if end:
            self.end = end
        else:
            self.end = len(self.vals)
        if masked:
            self.masked = None
        else:
            self.masked = masked
        
    def __repr__(self):
        return f'Run status: {self.has_run} \
                with cosmology: {self.cosmo} \
                with Photo-z error model:  \
                - bias:  {self.zbias} \
                - variance {[v*0.05 for v in self.zvariance]} \
                - outliers: {[o*0.03 for o in self.outliers]}'
        
    def _makeLensPZ(self):
        
        bins = np.linspace(0.2, 1.2, 11)
        self.gbias_dict = {}
        bin_centers = [.5*fsum([bins[i]+bins[i+1]]) for i in range(len(bins[:-1]))]
        self.pdf_z = SmailZ(self.zmid, np.array(self.dneff))
        dNdz_dict_lens = {}
        qs = []
        for index, (x,x2) in enumerate(zip(bins[:-1], bins[1:])):
            core = Core(self.zmid, zbias=0, sigma_z=0.03)
            tomofilter = uniform.pdf(self.zmid, loc=x, scale=x2-x)
            photoz_model = PhotozModel(self.pdf_z, core, [tomofilter])
            dNdz_dict_lens[bin_centers[index]] = photoz_model.get_pdf()[0]
            self.gbias_dict[bin_centers[index]] = self.gbias[index]
        self.dNdz_dict_lens = dNdz_dict_lens
        
        
    def _NormalizePZ(self, qs, dNdz_dict_source, m=1):
        for q, k in zip(qs, dNdz_dict_source.keys()):
            dNdz_dict_source[k] = dNdz_dict_source[k]*sum(qs)/q
            f = CubicSpline(self.zmid, dNdz_dict_source[k])
            d = quad(f, 0, 4)[0]
            for k in dNdz_dict_source.keys():
                dNdz_dict_source[k] /= (d*m)
            return dNdz_dict_source
        
    
    def _makeSourcePZ(self, implement_outliers=True):
        # print('Making Source Photo-z')
        n = len(self.zmid)
        datapts = ([list(np.ones(int(self.dneff[i]/min(self.dneff)))*self.zmid[i]) for i in range(n)])
        datapts = list(chain.from_iterable(datapts)) # flatten
        bins = datapts[0::int(len(datapts)/5)]
        self.bins = bins
        bin_centers = [.5*fsum([bins[i]+bins[i+1]]) for i in range(len(bins[:-1]))]
        self.bin_centers = bin_centers
        self.pdf_z = SmailZ(self.zmid, np.array(self.dneff))
        dNdz_dict_source = {}
        qs = []
        for index, (x,x2) in enumerate(zip(bins[:-1], bins[1:])):
            bias = self.zbias[index]
            variance = self.zvariance[index]
            core = Core(self.zmid, zbias=0.1*bias, sigma_z=variance*0.05)
            tomofilter = uniform.pdf(self.zmid, loc=x, scale=x2-x)
            photoz_model = PhotozModel(self.pdf_z, core, [tomofilter])
            dNdz_dict_source[bin_centers[index]] = photoz_model.get_pdf()[0]
            f = CubicSpline(self.zmid, dNdz_dict_source[bin_centers[index]])
            qs.append(quad(f, 0, 4)[0])
        norm = 1
        
        if implement_outliers:
            # inbin = [0.009780926739896716, 0.0022554556214481247, 0.020947481173020983, 0.04831089597165704, 0.07238322930340861]
            inbin = [0.03]*5
            inbin = [inbin[i]*self.outliers[i] for i in range(len(inbin))]
            norm = 1 - sum(inbin)

        dNdz_dict_source = self._NormalizePZ(qs, dNdz_dict_source, 5/norm)
        self.dNdz_dict_source = dNdz_dict_source
        
        if implement_outliers:
            self.scores = []
            self.KDEs = pickle.load(open('KDEs.p', 'rb'))
            for i, b in enumerate(list(sorted(self.dNdz_dict_source.keys()))):
                kde = self.KDEs[i]
                self.scores.append(np.exp(kde.score_samples(np.array(self.zmid).reshape(-1, 1))))
                self.dNdz_dict_source[b] += self.scores[i]*inbin[i]
            
    def getElls(self, file='ell-values.txt'):
        print('Getting Ells')
        ell = pd.read_csv(file, names=['ell'])
        self.ell = list(ell.to_dict()['ell'].values())
        
    def A_h(self, z, etah):
        if z>2:
            return (1+z)**etah
        return 1

        
    def A_l(self, z, etal):
        return (1+z)**etal
        
    
    def makeShearShearCells(self, cosmo=None):
        if not cosmo:
            cosmo = self.cosmo
        ccl_cls = pd.DataFrame()
        zbin = 0
        j = 0
        ia0 = self.A0 * np.array([self.A_l(zi, self.etal) for zi in self.zmid]) * np.array([self.A_h(zi, self.etah) for zi in self.zmid])
        lst = list(self.dNdz_dict_source.keys())
        for i, key in enumerate(lst):
            ia = self.getAi(self.beta, cosmo, dNdz=tuple(self.dNdz_dict_source[key])) * ia0
            lens1 = ccl.WeakLensingTracer(cosmo, dndz=(self.zmid, self.dNdz_dict_source[key]), ia_bias=(self.zmid, ia))
            for keyj in lst[i:]:
                ia = self.getAi(self.beta, cosmo, dNdz=tuple(self.dNdz_dict_source[keyj])) * ia0
                lens2 = ccl.WeakLensingTracer(cosmo, dndz=(self.zmid, self.dNdz_dict_source[keyj]), ia_bias=(self.zmid, ia))
                cls = ccl.angular_cl(cosmo, lens1, lens2, self.ell)
                newdf = pd.DataFrame({'zbin': [int(k) for k in j*np.ones(len(cls))],
                                      'ell': self.ell,
                                      'C_ell': cls})
                ccl_cls = pd.concat((ccl_cls, newdf))
                j += 1


        self.shearshear_cls = ccl_cls.reset_index()
        
        C_ells = []
        for i in set(ccl_cls['zbin']):
            C_ells.append(list(ccl_cls[ccl_cls['zbin']==i]['C_ell']))
        return C_ells
    
    def makePosPosCells(self, cosmo=None):
        if not cosmo:
            cosmo = self.cosmo
        ccl_cls = pd.DataFrame()
        zbin = 0
        j = 0
        lst = list(self.dNdz_dict_lens.keys())
        for i, key in enumerate(lst):
            pos1 = ccl.NumberCountsTracer(cosmo, dndz=(self.zmid, self.dNdz_dict_lens[key]), has_rsd=False, bias=(self.zmid, self.gbias[i]*np.ones_like(self.zmid)))
            cls = ccl.angular_cl(cosmo, pos1, pos1, self.ell)
            newdf = pd.DataFrame({'zbin': [int(k) for k in i*np.ones(len(cls))],
                                  'ell': self.ell,
                                  'C_ell': cls})
            ccl_cls = pd.concat((ccl_cls, newdf))


        self.pospos_cls = ccl_cls.reset_index()
        
        C_ells = []
        for i in set(ccl_cls['zbin']):
            C_ells.append(list(ccl_cls[ccl_cls['zbin']==i]['C_ell']))
        return C_ells


    def makePosShearCells(self, cosmo=None):
        if not cosmo:
            cosmo = self.cosmo
        ccl_cls = pd.DataFrame()
        zbin = 0
        j = 0
        llst = list(self.dNdz_dict_lens.keys())
        slst = list(self.dNdz_dict_source.keys())
        self.accept = {(0,1), (0,2), (0,3), (0,4),
                  (1,1), (1,2), (1,3), (1,4),
                  (2,2), (2,3), (2,4),
                  (3,2), (3,3), (3,4),
                  (4,2), (4,3), (4,4),
                  (5,3), (5,4),
                  (6,3), (6,4),
                  (7,3), (7,4),
                  (8,4), 
                  (9,4)
                  }
        for l, key in enumerate(llst):
            pos = ccl.NumberCountsTracer(cosmo, dndz=(self.zmid, self.dNdz_dict_lens[key]), has_rsd=False, bias=(self.zmid, self.gbias[l]*np.ones_like(self.zmid)))
            for s, keyj in enumerate(slst):
                if (l, s) in self.accept:
                    ia0 = self.A0 * np.array([self.A_l(zi, self.etal) for zi in self.zmid]) * np.array([self.A_h(zi, self.etah) for zi in self.zmid])
                    ia = self.getAi(self.beta, cosmo, dNdz=self.dNdz_dict_source[keyj]) * ia0
                    shear = ccl.WeakLensingTracer(cosmo, dndz=(self.zmid, self.dNdz_dict_source[keyj]), ia_bias=(self.zmid, ia))
                    cls = ccl.angular_cl(cosmo, pos, shear, self.ell)                

                    newdf = pd.DataFrame({'zbin': [int(o) for o in j*np.ones(len(cls))],
                                          'ell': self.ell,
                                          'C_ell': cls})
                    ccl_cls = pd.concat((ccl_cls, newdf))
                    j += 1


        self.posshear_cls = ccl_cls.reset_index()
        
        C_ells = []
        for i in set(ccl_cls['zbin']):
            C_ells.append(list(ccl_cls[ccl_cls['zbin']==i]['C_ell']))
        return C_ells

    def makeEmceeCells(self, cosmo):
        C_ells = []
        lst = list(self.dNdz_dict_source.keys())
        ia0 = self.A0 * np.array([self.A_l(zi, self.etal) for zi in self.zmid]) * np.array([self.A_h(zi, self.etah) for zi in self.zmid])
        for i, key in enumerate(lst):
            ia = self.getAi(self.beta, cosmo, dNdz=tuple(self.dNdz_dict_source[key])) * ia0
            lens1 = ccl.WeakLensingTracer(cosmo, dndz=(self.zmid, self.dNdz_dict_source[key]), ia_bias=(self.zmid, ia))
            for keyj in lst[i:]:
                ia = self.getAi(self.beta, cosmo, dNdz=tuple(self.dNdz_dict_source[keyj])) * ia0
                lens2 = ccl.WeakLensingTracer(cosmo, dndz=(self.zmid, self.dNdz_dict_source[keyj]), ia_bias=(self.zmid, ia))
                cls = ccl.angular_cl(cosmo, lens1, lens2, self.ell)
                C_ells.append(cls)

        return C_ells
    
    def buildCovMatrix(self):
        print('Getting covariance matrix')
        invcov_SRD = pd.read_csv('Y10_3x2pt_inv.txt', names=['a','b'], delimiter=' ')
        self.invcov = np.array(invcov_SRD['b']).reshape(1000, 1000)
        
        
    def getDerivs(self):
        print('Getting derivatives')

        self.derivs_sig = {}
        for var in self.param_order[:self.end]:
            print(var)
            zbin = 0
            j = 0
            derivs1, derivs2, derivs3 = [], [], []
            slst = list(self.dNdz_dict_source.keys())
            llst = list(self.dNdz_dict_lens.keys())
            if var not in self.funcs_ss.keys():
                derivs1 = list(np.zeros_like(self.ShearShearFid))
            else:
                for i, key in enumerate(slst):    
                    for keyj in slst[i:]:
                        self.key = key
                        self.keyj = keyj
                        f = nd.Derivative(self.funcs_ss[var], full_output=True, step=0.01)
                        val, info = f(self.vals[var])
                        derivs1.append(val)
                derivs1 = np.array(derivs1)
            if var not in self.funcs_sp.keys():
                derivs2 = list(np.zeros_like(self.PosShearFid))
            else:
                for l, keyl in enumerate(llst):
                    for s, keys in enumerate(slst):
                        if (l, s) in self.accept:
                            self.keyl = keyl
                            self.keys = keys 
                            f = nd.Derivative(self.funcs_sp[var], full_output=True, step=0.01)
                            val, info = f(self.vals[var])
                            derivs2.append(val)
                derivs2 = np.array(derivs2) 
            if var not in self.funcs_pp.keys():
                derivs3 = list(np.zeros_like(self.PosPosFid))
            else:
                for i, key in enumerate(llst):
                    self.key = key
                    f = nd.Derivative(self.funcs_pp[var], full_output=True, step=0.01)
                    val, info = f(self.vals[var])
                    derivs3.append(val)
            derivs3 = np.array(derivs3)
            self.derivs_sig[var] = np.vstack((derivs1, derivs2, derivs3))
            
    def getFisher(self):
        print('Building fisher matrix')
        fisher = np.zeros((self.end, self.end))
        derivs = deepcopy(self.derivs_sig)
        for var in self.param_order[:self.end]:
            if self.masked is not None:
                if 'ss' in self.masked:
                    derivs[var][:15] = 0
                if 'sl' in self.masked:
                    derivs[var][15:40] = 0
                if 'll' in self.masked:
                    derivs[var][40:] = 0
        for i, var1 in enumerate(self.param_order[:self.end]):
            for j, var2 in enumerate(self.param_order[:self.end]):
                fisher[i][j] = derivs[var1].flatten().T @ self.invcov @ derivs[var2].flatten()
                if self.calcCov:
                    fisher[i][j] *= len(self.ell)**2
        for i in range(self.end):
            fisher[i][i] += self.priors[self.param_order[i]]
        return fisher
        
    def process(self):
        self._makeSourcePZ()
        self._makeLensPZ()
        self.getElls()
        self.ShearShearFid = self.makeShearShearCells()
        self.PosShearFid = self.makePosShearCells()
        self.PosPosFid = self.makePosPosCells()
        self.ccl_cls = np.vstack((self.ShearShearFid, self.PosShearFid, self.PosPosFid))
        self.buildCovMatrix()
        self.getDerivs()
        self.fisher = self.getFisher()
        self.has_run = True
        print('Done')
        
        
#    @lru_cache
    def getAi(self, beta, cosmo, dNdz, Mr_s=-20.70, Q=1.23, alpha_lum=-1.23, phi_0=0.0094, P=-0.3, mlim=25.3, Lp= 1.):
        """ Get the amplitude of the 2-halo part of w_{l+}
        A0 is the amplitude, beta is the power law exponent (see Krause et al. 2016) 
        cosmo is a CCL cosmology object 
        Lp is a pivot luminosity (default = 1)
        dndz is (z, dNdz) - vector of z and dNdz for the galaxies
        (does not need to be normalised) 
        """
        z_input = self.zmid
        dNdz = list(dNdz)
        z = np.array(z_input)

        # Get the luminosity function
        (L, phi_normed) = self.get_phi(z, cosmo, Mr_s, Q, alpha_lum, phi_0, P, mlim)
        # Pivot luminosity:
        Lp = 1.

        # Get Ai as a function of lens redshift.
        if beta == 1:
            #print('no integration')
            dl = ccl.luminosity_distance(cosmo, 1./(1.+z))
            Ai_ofzl = self.IA_interp(dl)
        else:
            #print('integrating')
            Ai_ofzl = np.array([scipy.integrate.simps(np.asarray(phi_normed[zi]) * (np.asarray(L[zi]) / Lp)**(beta), np.asarray(L[zi])) for zi in range(len(z))])

        # Integrate over dNdz
        Ai = scipy.integrate.simps(Ai_ofzl * dNdz, z) / scipy.integrate.simps(dNdz, z)

        return Ai

    def get_phi(self, z, cosmo, Mr_s, Q, alpha_lum, phi_0, P, mlim, Mp=-22.):

        """ This function outputs the Schechter luminosity function with parameters fit in Loveday 2012, following the same procedure as Krause et al. 2015, as a function of z and L 
        The output is L[z][l], list of vectors of luminosity values in z, different at each z due to the different lower luminosity limit, and phi[z][l], a list of luminosity functions at these luminosity vectors, at each z
        cosmo is a CCL cosmology object
        mlim is the magnitude limit of the survey
        Mp is the pivot absolute magnitude.
        other parameteres are the parameters of the luminosity function that are different for different samples, e.g. red vs all. lumparams = [Mr_s, Q, alpha_lum, phi_0, P]
        Note that the luminosity function is output normalized (appropriate for getting Ai)."""

        # Get the amplitude of the Schechter luminosity function as a function of redshift.
        phi_s = phi_0 * 10.**(0.4 * P * z)

        # Get M_* (magnitude), then convert to L_*
        Ms = Mr_s - Q * (z - 0.1)
        Ls = 10**(-0.4 * (Ms - Mp))

        # Import the kcorr and ecorr correction from Poggianti (assumes elliptical galaxies)
        # No data for sources beyon z = 3, so we keep the same value at higher z as z=3
        (z_k, kcorr, x,x,x) = np.loadtxt('kcorr.dat', unpack=True)
        (z_e, ecorr, x,x,x) = np.loadtxt('ecorr.dat', unpack=True)
        kcorr_interp = CubicSpline(z_k, kcorr)
        ecorr_interp = CubicSpline(z_e, ecorr)
        kcorr = kcorr_interp(z)
        ecorr = ecorr_interp(z)

        # Get the absolute magnitude and luminosity corresponding to limiting apparent magntiude (as a function of z)
        dl = ccl.luminosity_distance(cosmo, 1./(1.+z))
        Mlim = mlim - (5. * np.log10(dl) + 25. + kcorr + ecorr)
        Llim = 10.**(-0.4 * (Mlim-Mp))


        L = [scipy.logspace(np.log10(Llim[zi]), 2., 1000) for zi in range(len(z))]

        # Now get phi(L,z), where this exists for each z because the lenghts of the L vectors are different.
        phi_func = [0]*len(z)
        for zi in range(len(z)):
            phi_func[zi] = [phi_s[zi] * (L[zi][li] / Ls[zi]) ** (alpha_lum) * np.exp(- L[zi][li] / Ls[zi]) for li in range(len(L[zi]))]

        norm = np.zeros(len(z))
        phi_func_normed = [0]*len(z)
        for zi in range(len(z)):
            norm[zi] = scipy.integrate.simps(phi_func[zi], L[zi])
            phi_func_normed[zi] = phi_func[zi] / norm[zi]

        return (L, phi_func_normed)
        
    def _outlier_helper(self, idx, zoutlier):
        dNdz_dict_source = {}
        qs = []
        for index, (x,x2) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            core = Core(self.zmid, zbias=0.1*self.zbias[index], sigma_z=self.zvariance[index]*0.05)
            tomofilter = uniform.pdf(self.zmid, loc=x, scale=x2-x)
            photoz_model = PhotozModel(self.pdf_z, core, [tomofilter])
            dNdz_dict_source[self.bin_centers[index]] = photoz_model.get_pdf()[0]
            f = CubicSpline(self.zmid, dNdz_dict_source[self.bin_centers[index]])
            qs.append(quad(f, 0, 4)[0])

        inbin = [0.03]*5
        for index in range(len(self.bins)-1):
            if index==idx: 
                inbin[index] = inbin[index]*self.outliers[index]*zoutlier
            else:
                inbin[index] = inbin[index]*self.outliers[index]
        norm = 1 - sum(inbin)

        dNdz_dict_source = self._NormalizePZ(qs, dNdz_dict_source, 5/norm)

        for i, b in enumerate(list(sorted(dNdz_dict_source.keys()))):
            dNdz_dict_source[b] += self.scores[i]*inbin[i]
        return dNdz_dict_source

    def getC_ellOfzoutlier1_ss(self, zoutlier):
        index = 0
        dNdz_dict_source = self._outlier_helper(index, zoutlier)
        return self._helper_ss(self.cosmo, dNdz_dict_source)
    
    def getC_ellOfzoutlier2_ss(self, zoutlier):
        index = 1
        dNdz_dict_source = self._outlier_helper(index, zoutlier)
        return self._helper_ss(self.cosmo, dNdz_dict_source)

    def getC_ellOfzoutlier3_ss(self, zoutlier):
        index = 2
        dNdz_dict_source = self._outlier_helper(index, zoutlier)
        return self._helper_ss(self.cosmo, dNdz_dict_source)

    def getC_ellOfzoutlier4_ss(self, zoutlier):
        index = 3
        dNdz_dict_source = self._outlier_helper(index, zoutlier)
        return self._helper_ss(self.cosmo, dNdz_dict_source)

    def getC_ellOfzoutlier5_ss(self, zoutlier):
        index = 4
        dNdz_dict_source = self._outlier_helper(index, zoutlier)
        return self._helper_ss(self.cosmo, dNdz_dict_source)
    
    def getC_ellOfzoutlier1_sp(self, zoutlier):
        index = 0
        dNdz_dict_source = self._outlier_helper(index, zoutlier)
        return self._helper_sp(self.cosmo, dNdz_dict_source)
    
    def getC_ellOfzoutlier2_sp(self, zoutlier):
        index = 1
        dNdz_dict_source = self._outlier_helper(index, zoutlier)
        return self._helper_sp(self.cosmo, dNdz_dict_source)

    def getC_ellOfzoutlier3_sp(self, zoutlier):
        index = 2
        dNdz_dict_source = self._outlier_helper(index, zoutlier)
        return self._helper_sp(self.cosmo, dNdz_dict_source)

    def getC_ellOfzoutlier4_sp(self, zoutlier):
        index = 3
        dNdz_dict_source = self._outlier_helper(index, zoutlier)
        return self._helper_sp(self.cosmo, dNdz_dict_source)

    def getC_ellOfzoutlier5_sp(self, zoutlier):
        index = 4
        dNdz_dict_source = self._outlier_helper(index, zoutlier)
        return self._helper_sp(self.cosmo, dNdz_dict_source)
    
    def _bias_helper(self, idx, zbias):
        dNdz_dict_source = {}
        qs = []
        for index, (x,x2) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            if index==idx:
                core = Core(self.zmid, zbias=self.zbias[index]+zbias, sigma_z=0.05*self.zvariance[index])
            else:
                core = Core(self.zmid, zbias=self.zbias[index], sigma_z=0.05*self.zvariance[index])
            tomofilter = uniform.pdf(self.zmid, loc=x, scale=x2-x)
            photoz_model = PhotozModel(self.pdf_z, core, [tomofilter])
            dNdz_dict_source[self.bin_centers[index]] = photoz_model.get_pdf()[0]
            f = CubicSpline(self.zmid, dNdz_dict_source[self.bin_centers[index]])
            qs.append(quad(f, 0, 4)[0])

        inbin = [0.03]*5
        for index in range(len(self.bins)-1):
            inbin[index] = inbin[index]*self.outliers[index]
        norm = 1 - sum(inbin)

        dNdz_dict_source = self._NormalizePZ(qs, dNdz_dict_source, 5/norm)
        
        for i, b in enumerate(list(sorted(self.dNdz_dict_source.keys()))):
            dNdz_dict_source[b] += self.scores[i]*inbin[i]
        return dNdz_dict_source
    
    
    def getC_ellOfIA_amp_ss(self, A0):
        return self._helper_ss(self.cosmo, self.dNdz_dict_source, A0=A0)
        
    def getC_ellOfIA_highz_ss(self, etah):
        return self._helper_ss(self.cosmo, self.dNdz_dict_source, etah=etah)
        
    def getC_ellOfIA_lowz_ss(self, etal):
        return self._helper_ss(self.cosmo, self.dNdz_dict_source, etal=etal)
        
    def getC_ellOfIA_beta_ss(self, beta):
        return self._helper_ss(self.cosmo, self.dNdz_dict_source, beta=beta)


    def getC_ellOfIA_amp_sp(self, A0):
        return self._helper_sp(self.cosmo, self.dNdz_dict_source, A0=A0)
        
    def getC_ellOfIA_highz_sp(self, etah):
        return self._helper_sp(self.cosmo, self.dNdz_dict_source, etah=etah)
        
    def getC_ellOfIA_lowz_sp(self, etal):
        return self._helper_sp(self.cosmo, self.dNdz_dict_source, etal=etal)
        
    def getC_ellOfIA_beta_sp(self, beta):
        return self._helper_sp(self.cosmo, self.dNdz_dict_source, beta=beta)

    def getC_ellOfzbias1_ss(self, zbias):
        index = 0
        dNdz_dict_source = self._bias_helper(index, zbias)
        return self._helper_ss(self.cosmo, dNdz_dict_source)
    
    def getC_ellOfzbias2_ss(self, zbias):
        index = 1
        dNdz_dict_source = self._bias_helper(index, zbias)
        return self._helper_ss(self.cosmo, dNdz_dict_source)

    def getC_ellOfzbias3_ss(self, zbias):
        index = 2
        dNdz_dict_source = self._bias_helper(index, zbias)
        return self._helper_ss(self.cosmo, dNdz_dict_source)

    def getC_ellOfzbias4_ss(self, zbias):
        index = 3
        dNdz_dict_source = self._bias_helper(index, zbias)
        return self._helper_ss(self.cosmo, dNdz_dict_source)

    def getC_ellOfzbias5_ss(self, zbias):
        index = 4
        dNdz_dict_source = self._bias_helper(index, zbias)
        return self._helper_ss(self.cosmo, dNdz_dict_source)

    def getC_ellOfzbias1_sp(self, zbias):
        index = 0
        dNdz_dict_source = self._bias_helper(index, zbias)
        return self._helper_sp(self.cosmo, dNdz_dict_source)
    
    def getC_ellOfzbias2_sp(self, zbias):
        index = 1
        dNdz_dict_source = self._bias_helper(index, zbias)
        return self._helper_sp(self.cosmo, dNdz_dict_source)

    def getC_ellOfzbias3_sp(self, zbias):
        index = 2
        dNdz_dict_source = self._bias_helper(index, zbias)
        return self._helper_sp(self.cosmo, dNdz_dict_source)

    def getC_ellOfzbias4_sp(self, zbias):
        index = 3
        dNdz_dict_source = self._bias_helper(index, zbias)
        return self._helper_sp(self.cosmo, dNdz_dict_source)

    def getC_ellOfzbias5_sp(self, zbias):
        index = 4
        dNdz_dict_source = self._bias_helper(index, zbias)
        return self._helper_sp(self.cosmo, dNdz_dict_source)

    def _variance_helper(self, idx, zvar):
        dNdz_dict_source = {}
        qs = []
        for index, (x,x2) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            if index==idx:
                core = Core(self.zmid, zbias=self.zbias[index], sigma_z=0.05*self.zvariance[index]*zvar)
            else:
                core = Core(self.zmid, zbias=self.zbias[index], sigma_z=0.05*self.zvariance[index])
            tomofilter = uniform.pdf(self.zmid, loc=x, scale=x2-x)
            photoz_model = PhotozModel(self.pdf_z, core, [tomofilter])
            dNdz_dict_source[self.bin_centers[index]] = photoz_model.get_pdf()[0]
            f = CubicSpline(self.zmid, dNdz_dict_source[self.bin_centers[index]])
            qs.append(quad(f, 0, 4)[0])

        inbin = [0.03]*5
        for index in range(len(self.bins)-1):
            inbin[index] = inbin[index]*self.outliers[index]
        norm = 1 - sum(inbin)

        dNdz_dict_source = self._NormalizePZ(qs, dNdz_dict_source, 5/norm)

        for i, b in enumerate(list(sorted(dNdz_dict_source.keys()))):
            dNdz_dict_source[b] += self.scores[i]*inbin[i]
        return dNdz_dict_source

    def getC_ellOfzvariance1_ss(self, zvariance):
        index = 0
        dNdz_dict_source = self._variance_helper(index, zvariance)
        return self._helper_ss(self.cosmo, dNdz_dict_source)

    def getC_ellOfzvariance2_ss(self, zvariance):
        index = 1
        dNdz_dict_source = self._variance_helper(index, zvariance)
        return self._helper_ss(self.cosmo, dNdz_dict_source)

    def getC_ellOfzvariance3_ss(self, zvariance):
        index = 2
        dNdz_dict_source = self._variance_helper(index, zvariance)
        return self._helper_ss(self.cosmo, dNdz_dict_source)

    def getC_ellOfzvariance4_ss(self, zvariance):
        index = 3
        dNdz_dict_source = self._variance_helper(index, zvariance)
        return self._helper_ss(self.cosmo, dNdz_dict_source)

    def getC_ellOfzvariance5_ss(self, zvariance):
        index = 4
        dNdz_dict_source = self._variance_helper(index, zvariance)
        return self._helper_ss(self.cosmo, dNdz_dict_source)

    def getC_ellOfzvariance1_sp(self, zvariance):
        index = 0
        dNdz_dict_source = self._variance_helper(index, zvariance)
        return self._helper_sp(self.cosmo, dNdz_dict_source)

    def getC_ellOfzvariance2_sp(self, zvariance):
        index = 1
        dNdz_dict_source = self._variance_helper(index, zvariance)
        return self._helper_sp(self.cosmo, dNdz_dict_source)

    def getC_ellOfzvariance3_sp(self, zvariance):
        index = 2
        dNdz_dict_source = self._variance_helper(index, zvariance)
        return self._helper_sp(self.cosmo, dNdz_dict_source)

    def getC_ellOfzvariance4_sp(self, zvariance):
        index = 3
        dNdz_dict_source = self._variance_helper(index, zvariance)
        return self._helper_sp(self.cosmo, dNdz_dict_source)

    def getC_ellOfzvariance5_sp(self, zvariance):
        index = 4
        dNdz_dict_source = self._variance_helper(index, zvariance)
        return self._helper_sp(self.cosmo, dNdz_dict_source)

    def _helper_ss(self, cosmo, dNdz_dict_source, A0=None, beta=None, etal=None, etah=None):
        if not beta:
            beta = self.beta
        if not etal:
            etal = self.etal
        if not etah:
            etah = self.etah
        if not A0:
            A0 = self.A0
        ia0 =  A0 * np.array([self.A_l(zi, etal) for zi in self.zmid]) * np.array([self.A_h(zi, etah) for zi in self.zmid])
        ia = self.getAi(beta, cosmo, dNdz=tuple(dNdz_dict_source[self.key])) * ia0
        lens1 = ccl.WeakLensingTracer(cosmo, dndz=(self.zmid, dNdz_dict_source[self.key]), ia_bias=(self.zmid, ia))
        ia = self.getAi(beta, cosmo, dNdz=tuple(dNdz_dict_source[self.keyj])) * ia0
        lens2 = ccl.WeakLensingTracer(cosmo, dndz=(self.zmid, dNdz_dict_source[self.keyj]), ia_bias=(self.zmid, ia))
        return ccl.angular_cl(cosmo, lens1, lens2, self.ell)

    def _helper_sp(self, cosmo, dNdz_dict_source, A0=None, beta=None, etal=None, etah=None):
        if not beta:
            beta = self.beta
        if not etal:
            etal = self.etal
        if not etah:
            etah = self.etah
        if not A0:
            A0 = self.A0
        pos = ccl.NumberCountsTracer(self.cosmo, dndz=(self.zmid, self.dNdz_dict_lens[self.keyl]), has_rsd=False, bias=(self.zmid, self.gbias_dict[self.keyl]*np.ones_like(self.zmid)))
        ia0 =  A0 * np.array([self.A_l(zi, etal) for zi in self.zmid]) * np.array([self.A_h(zi, etah) for zi in self.zmid])
        ia = self.getAi(beta, cosmo, dNdz=tuple(dNdz_dict_source[self.keys])) * ia0
        lens = ccl.WeakLensingTracer(cosmo, dndz=(self.zmid, dNdz_dict_source[self.keys]), ia_bias=(self.zmid, ia))
        return ccl.angular_cl(cosmo, pos, lens, self.ell)


    def getC_ellOfgbias_ss(self, gbias):
        pos1 = ccl.NumberCountsTracer(self.cosmo, dndz=(self.zmid, self.dNdz_dict_lens[self.key]), has_rsd=False, bias=(self.zmid, gbias*np.ones_like(self.zmid)))
        return ccl.angular_cl(self.cosmo, pos1, pos1, self.ell)


    def _helper_pp(self, cosmo):
        pos1 = ccl.NumberCountsTracer(cosmo, dndz=(self.zmid, self.dNdz_dict_lens[self.key]), has_rsd=False, bias=(self.zmid, self.gbias_dict[self.key]*np.ones_like(self.zmid)))
        return ccl.angular_cl(cosmo, pos1, pos1, self.ell)

    def getC_ellOfSigma8_ss(self, sigma8):
        cosmo = ccl.Cosmology(Omega_c=0.2666, Omega_b=0.049, h=0.6727, sigma8=sigma8, n_s=0.9645, 
            transfer_function='eisenstein_hu')
        return self._helper_ss(cosmo, dNdz_dict_source=self.dNdz_dict_source)

    def getC_ellOfOmegab_ss(self, Omega_b):
        cosmo = ccl.Cosmology(Omega_c=0.2666, Omega_b=Omega_b, h=0.6727, sigma8=0.831, n_s=0.9645, 
            transfer_function='eisenstein_hu')
        return self._helper_ss(cosmo, dNdz_dict_source=self.dNdz_dict_source)

    def getC_ellOfh_ss(self, h):
        cosmo = ccl.Cosmology(Omega_c=0.2666, Omega_b=0.049, h=h, sigma8=0.831, n_s=0.9645, 
            transfer_function='eisenstein_hu')
        return self._helper_ss(cosmo, dNdz_dict_source=self.dNdz_dict_source)

    def getC_ellOfn_s_ss(self, n_s):
        cosmo = ccl.Cosmology(Omega_c=0.2666, Omega_b=0.049, h=0.6727, sigma8=0.831, n_s=n_s, 
            transfer_function='eisenstein_hu')
        return self._helper_ss(cosmo, dNdz_dict_source=self.dNdz_dict_source)


    def getC_ellOfOmegam_ss(self, Omega_m):
        cosmo = ccl.Cosmology(Omega_c=Omega_m-0.049, Omega_b=0.049, h=0.6727, sigma8=0.831, n_s=0.9645, 
            transfer_function='eisenstein_hu')
        return self._helper_ss(cosmo, dNdz_dict_source=self.dNdz_dict_source)


    def getC_ellOfw0_ss(self, w_0):
        cosmo = ccl.Cosmology(Omega_c=0.2666, Omega_b=0.049, h=0.6727, sigma8=0.831, n_s=0.9645, w0=w_0, 
            transfer_function='eisenstein_hu')
        return self._helper_ss(cosmo, dNdz_dict_source=self.dNdz_dict_source)


    def getC_ellOfwa_ss(self, w_a):
        cosmo = ccl.Cosmology(Omega_c=0.2666, Omega_b=0.049, h=0.6727, sigma8=0.831, n_s=0.9645, wa=w_a, 
            transfer_function='eisenstein_hu')
        return self._helper_ss(cosmo, dNdz_dict_source=self.dNdz_dict_source)


    def getC_ellOfSigma8_sp(self, sigma8):
        cosmo = ccl.Cosmology(Omega_c=0.2666, Omega_b=0.049, h=0.6727, sigma8=sigma8, n_s=0.9645, 
            transfer_function='eisenstein_hu')
        return self._helper_sp(cosmo, dNdz_dict_source=self.dNdz_dict_source)

    def getC_ellOfOmegab_sp(self, Omega_b):
        cosmo = ccl.Cosmology(Omega_c=0.2666, Omega_b=Omega_b, h=0.6727, sigma8=0.831, n_s=0.9645, 
            transfer_function='eisenstein_hu')
        return self._helper_sp(cosmo, dNdz_dict_source=self.dNdz_dict_source)

    def getC_ellOfh_sp(self, h):
        cosmo = ccl.Cosmology(Omega_c=0.2666, Omega_b=0.049, h=h, sigma8=0.831, n_s=0.9645, 
            transfer_function='eisenstein_hu')
        return self._helper_sp(cosmo, dNdz_dict_source=self.dNdz_dict_source)

    def getC_ellOfn_s_sp(self, n_s):
        cosmo = ccl.Cosmology(Omega_c=0.2666, Omega_b=0.049, h=0.6727, sigma8=0.831, n_s=n_s, 
            transfer_function='eisenstein_hu')
        return self._helper_sp(cosmo, dNdz_dict_source=self.dNdz_dict_source)


    def getC_ellOfOmegam_sp(self, Omega_m):
        cosmo = ccl.Cosmology(Omega_c=Omega_m-0.049, Omega_b=0.049, h=0.6727, sigma8=0.831, n_s=0.9645, 
            transfer_function='eisenstein_hu')
        return self._helper_sp(cosmo, dNdz_dict_source=self.dNdz_dict_source)


    def getC_ellOfw0_sp(self, w_0):
        cosmo = ccl.Cosmology(Omega_c=0.2666, Omega_b=0.049, h=0.6727, sigma8=0.831, n_s=0.9645, w0=w_0, 
            transfer_function='eisenstein_hu')
        return self._helper_sp(cosmo, dNdz_dict_source=self.dNdz_dict_source)


    def getC_ellOfwa_sp(self, w_a):
        cosmo = ccl.Cosmology(Omega_c=0.2666, Omega_b=0.049, h=0.6727, sigma8=0.831, n_s=0.9645, wa=w_a, 
            transfer_function='eisenstein_hu')
        return self._helper_sp(cosmo, dNdz_dict_source=self.dNdz_dict_source)


    def getC_ellOfSigma8_pp(self, sigma8):
        cosmo = ccl.Cosmology(Omega_c=0.2666, Omega_b=0.049, h=0.6727, sigma8=sigma8, n_s=0.9645, 
            transfer_function='eisenstein_hu')
        return self._helper_pp(cosmo)

    def getC_ellOfOmegab_pp(self, Omega_b):
        cosmo = ccl.Cosmology(Omega_c=0.2666, Omega_b=Omega_b, h=0.6727, sigma8=0.831, n_s=0.9645, 
            transfer_function='eisenstein_hu')
        return self._helper_pp(cosmo)

    def getC_ellOfh_pp(self, h):
        cosmo = ccl.Cosmology(Omega_c=0.2666, Omega_b=0.049, h=h, sigma8=0.831, n_s=0.9645, 
            transfer_function='eisenstein_hu')
        return self._helper_pp(cosmo)

    def getC_ellOfn_s_pp(self, n_s):
        cosmo = ccl.Cosmology(Omega_c=0.2666, Omega_b=0.049, h=0.6727, sigma8=0.831, n_s=n_s, 
            transfer_function='eisenstein_hu')
        return self._helper_pp(cosmo)


    def getC_ellOfOmegam_pp(self, Omega_m):
        cosmo = ccl.Cosmology(Omega_c=Omega_m-0.049, Omega_b=0.049, h=0.6727, sigma8=0.831, n_s=0.9645, 
            transfer_function='eisenstein_hu')
        return self._helper_pp(cosmo)


    def getC_ellOfw0_pp(self, w_0):
        cosmo = ccl.Cosmology(Omega_c=0.2666, Omega_b=0.049, h=0.6727, sigma8=0.831, n_s=0.9645, w0=w_0, 
            transfer_function='eisenstein_hu')
        return self._helper_pp(cosmo)


    def getC_ellOfwa_pp(self, w_a):
        cosmo = ccl.Cosmology(Omega_c=0.2666, Omega_b=0.049, h=0.6727, sigma8=0.831, n_s=0.9645, wa=w_a, 
            transfer_function='eisenstein_hu')
        return self._helper_pp(cosmo)