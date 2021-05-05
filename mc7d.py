import sys
sys.path = ['/global/u2/h/husni/PZ_Project',
 '/global/u2/h/husni/PZ_Project',
 '',
 '/global/homes/h/husni/.local/lib/python3.8/site-packages',
 '/global/homes/h/husni/.local/bin',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/sims_almanac/master-g021b69e146+616205b9df/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/sims_downtimeModel/master-g55f72efa65+3e384ed16a/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/sims_cloudModel/master-ge3724df529+632df0f48d/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/sims_seeingModel/master-ga4bf72ea44+6e62bc95d5/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/sims_featureScheduler/master-g5f42dc1c76+c0a10aa3f3/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/sims_skybrightness_pre/2.13.0.sims-11-g9ab127b+5f00cb8bf4/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/sims_skybrightness_data/2017.05.05-2-gf1b2499/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/sims_skybrightness/2.13.0.sims-7-g57b020e+ba01f1af57/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/python_extinction/master-g3504403c44+15e725ac18/lib/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/astropy_helpers/3.0.2.lsst-1-g9c76ab2/lib/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/astropy_helpers/3.0.2.lsst-1-g9c76ab2',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/sncosmo/1.3.0.lsst2-5-gea39e10+3c1685c4f7/lib/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/sims_catUtils/2.13.0.sims-16-g9f88571c+4cf1f4b6fb/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/sims_survey_fields/2.13.0.sims-1-g3f5255d+616205b9df/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/scarlet_extensions/lsst-dev-g9d18589735+5de2e29f4b/lib/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/sims_coordUtils/2.13.0.sims-7-gc9b9334+d5a977bfcb/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/sims_photUtils/2.13.0.sims-5-ge14a266+aba1a31c29/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/pymssql/2.1.1-2-g168068f+ad9f4941ed/lib/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/sims_catalogs/2.13.0.sims-9-g8bb9136+47a81b2cae/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/sims_utils/2.13.0.sims-46-g370dc29+872406ca88/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/palpy/1.8.1-1-ga780397/lib/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/ephem/3.7.6.0-1-g6182364/lib/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/sims_maf/2.13.0.sims-94-g33e25da1+7297bd3ce7/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/dustmaps_cachedata/master-gd5d683fb32+04719a4bac/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/faro/master-ge470644580+765c7fdc98/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/fgcmcal/21.0.0-14-g8629b32+ada1ffd932/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/validate_drp/21.0.0-3-g107f81d+8475e67332/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/sdm_schemas/21.0.0-1-ge52cb69+04719a4bac/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/cp_pipe/21.0.0-16-g434b3eb+a783293e9c/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/display_matplotlib/21.0.0-2-gfc62afb+2a0303cc17/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/display_firefly/21.0.0-6-gc675373+2a0303cc17/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/ap_verify/21.0.0-12-g786902c+5e4455bbc1/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/jointcal/21.0.0-13-g3a573fe+7df2a27b45/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/cbp/21.0.0-2-g143869c+2a0303cc17/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/obs_subaru/21.0.0-29-ge89e7495+5cf2761b2b/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/obs_lsst/21.0.0-33-g03ef56e+8cd588cbf3/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/obs_cfht/21.0.0-11-g29e88f6+d9d31afa03/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/ctrl_pool/21.0.0-2-ga326454+7857418d36/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/pipe_drivers/21.0.0-5-gd00fb1e+9635e151a8/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/meas_extensions_shapeHSM/21.0.0-2-gf484ad7+2dd4ab1e4d/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/meas_extensions_photometryKron/21.0.0-3-ge02ed75+2dd4ab1e4d/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/ctrl_bps/21.0.0-13-g0c9ff3e+5a28569b12/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/ctrl_mpexec/21.0.0-25-g2bc75c1+ab0f4fc97f/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/ctrl_orca/21.0.0-2-gfca2c14+fcd967797e/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/ctrl_execute/21.0.0-2-g3d7e797+38373445df/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/meas_extensions_piff/master-gac4afde19b+2dd4ab1e4d/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/display_ds9/21.0.0-3-g4a4ce7f+2a0303cc17/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/meas_extensions_simpleShape/21.0.0-2-g103fe59+e59df2cecd/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/obs_decam/21.0.0-15-g5a7caf0+0ea4a5585d/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/alert_packet/master-gbd80f2fea6/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/ap_association/21.0.0-9-g4828449+7bc8d59392/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/ap_pipe/21.0.0-12-g7d1023a+7582d51a30/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/verify_metrics/21.0.0-2-gfc76737+04719a4bac/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/dax_apdb/21.0.0-2-g5242d73+2a0303cc17/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/verify/21.0.0-8-ga0979ac+3aab7553ee/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/ip_diffim/21.0.0-12-g5009899+33d7f927f3/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/ip_isr/21.0.0-14-g07a2928+2dd4ab1e4d/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/meas_astrom/21.0.0-4-g591bb35+2dd4ab1e4d/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/kht/lsst-dev-g14ffe67dc2+fe93bf5141/cpp/build',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/psfex/21.0.0-2-g03166ea+04719a4bac/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/meas_extensions_psfex/21.0.0-3-g357aad2+5f7b813397/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/proxmin/lsst-dev-g79c0498783/lib/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/scarlet/lsst-dev-g965bb5fbbf+f31336177f/lib/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/meas_extensions_scarlet/21.0.0-12-g771b9a2+644d7c70fa/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/pipe_tasks/21.0.0-65-g9ea87ca1+ca40f17e88/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/shapelet/21.0.0-3-g4be5c26+2a0303cc17/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/meas_modelfit/21.0.0-2-gecfae73+c0704111f2/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/astro_metadata_translator/0.1.0-25-gffa9b84+f5e6047307/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/obs_base/21.0.0-48-g6f8f15d+42060c315e/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/obs_test/21.0.0-4-ge8a399c+7d01a1857e/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/pipe_base/21.0.0-17-g8839fe5+c812bf64e9/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/coadd_utils/21.0.0-2-g7f82c8f+7857418d36/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/daf_butler/21.0.0-69-g3ff10e88+06509c8b61/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/skymap/21.0.0-4-g65b4814+c812bf64e9/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/meas_base/21.0.0-8-g72414d6+5fcc31e360/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/meas_algorithms/21.0.0-20-g8cd22d88+ae1e48c0d5/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/astshim/21.0.0-2-g45278ab+04719a4bac/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/sphgeom/21.0.0+04719a4bac/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/geom/21.0.0-4-g05c9afb+06509c8b61/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/pex_config/21.0.0-1-ga51b5d4+f5e6047307/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/log/21.0.0-3-g7d9da8d+616205b9df/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/daf_persistence/21.0.0-8-g5674e7b+d1bd76f71f/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/pex_exceptions/21.0.0-2-gde069b7+5e4aea9c2f/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/sconsUtils/21.0.0-7-gd4ecef8/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/base/21.0.0-7-gdf92d54+04719a4bac/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/utils/21.0.0-4-gccdca77+0de219a2bc/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/daf_base/21.0.0-2-g8faa9b5+616205b9df/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/afw/21.0.0-27-g37a8c363b+1e633d4884/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/meas_deblender/21.0.0-3-g65f322c+6c73791a99/python',
 '/opt/lsst/software/stack/conda/miniconda3-py38_4.9.2/envs/lsst-scipipe-0.4.3/eups/python',
 '/opt/mods/lib/python3.6/site-packages',
 '/opt/ovis/lib/python3.6/site-packages',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/meas_extensions_convolved/21.0.0-2-g7713827+138fd9ef6c/python',
 '/opt/lsst/software/stack/stack/miniconda3-py38_4.9.2-0.4.3/Linux64/fgcm/lsst-dev-g14fa2066fb/lib/python',
 '/opt/lsst/software/stack/supreme',
 '/opt/lsst/software/stack/conda/miniconda3-py38_4.9.2/envs/lsst-scipipe-0.4.3/lib/python38.zip',
 '/opt/lsst/software/stack/conda/miniconda3-py38_4.9.2/envs/lsst-scipipe-0.4.3/lib/python3.8',
 '/opt/lsst/software/stack/conda/miniconda3-py38_4.9.2/envs/lsst-scipipe-0.4.3/lib/python3.8/lib-dynload',
 '/opt/lsst/software/stack/conda/miniconda3-py38_4.9.2/envs/lsst-scipipe-0.4.3/lib/python3.8/site-packages',
 '/opt/lsst/software/stack/conda/miniconda3-py38_4.9.2/envs/lsst-scipipe-0.4.3/lib/python3.8/site-packages/IPython/extensions',
 '/global/u2/h/husni/.ipython']

import numpy as np
import pandas as pd;
from collections import namedtuple
from pprint import pprint;
import sys;
from copy import deepcopy;
import pickle;
import pyccl as ccl
import emcee
from fisher import Fisher, FullPlot, marginalize, plot_contours
from multiprocessing import Pool


obj = pickle.load(open('/global/u2/h/husni/PZ_Project/zbiasedFisher.p','rb'))

data = np.array(obj.ccl_cls['C_ell']).reshape(15, 15).T
def lnprob(theta):
    cosmo = ccl.Cosmology(Omega_c=theta[1]-0.049, 
           Omega_b=0.049, 
           h=0.6727, 
           sigma8=theta[0], 
           n_s=0.9645, 
           transfer_function='eisenstein_hu')
    obj.bias = [theta[i] for i in range(2, 7)]
    obj._makeSourcePZ()
    theory = np.array(obj.makeEmceeCells(cosmo)).T;
    diff = np.array(data) - np.array(theory)
    return -sum(np.dot(diff[l], np.dot(obj.invcov_list[l], diff[l]))/2
               for l in range(15)) + lnprior(theta)

def lnprior(theta):
    X = {'sigma_8':theta[0], 'omega_m':theta[1]}
    mu = {'sigma_8':0.831, 'omega_m':0.3156}
    sigma = {'sigma_8':0.14, 'omega_m':0.2}
    for i in range(5):
        mu['zbias'+str(i)] = 0
        X['zbias'+str(i)] = theta[i+2]
        sigma['zbias'+str(i)] = 2
    return sum([
        np.log(1/(np.sqrt(2*np.pi)*sigma[i])) - (X[i]-mu[i])**2/(2*(sigma[i]**2)) for i in X.keys()
    ])


ndim = 7
nwalkers = 64
nsteps = 100000
fid = np.array([0.831, 0.3156, 0, 0, 0, 0, 0])
p00 = 0.1 * np.random.randn(nwalkers, ndim) + fid
filename = "/global/u2/h/husni/PZ_Project/mc7d.h5"
backend = emcee.backends.HDFBackend(filename)
backend.reset(nwalkers, ndim)
with Pool() as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool, backend=backend)
    sampler.run_mcmc(p00, nsteps)

