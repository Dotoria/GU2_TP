import numpy as np
import astropy.io.fits as fits
from astropy.cosmology import FlatLambdaCDM
from data_processing import data_processing as dat

class SAMI :
  # SAMI catalog data
  dtype_cat = [('CATID','<f8'),('type','<f8'),('z','<f8'),('color','<f8'),('mu','<f8'),('Re','<f8'),('mass','<f8'),('Mr','<f8'),('ra','<f8'),('dec','<f8')]
  GAMAdata = np.genfromtxt("GU2_HW4_GAMA_region.csv",dtype=dtype_cat,delimiter=',',skip_header=1,usecols=(1,2,3,4,9,10,7,8,5,6))
  clustdata = np.genfromtxt("GU2_HW4_cluster_region.csv",dtype=dtype_cat,delimiter=',',skip_header=1,usecols=(1,2,3,4,9,10,7,8,5,6))
   # Data processing
  GAMAdata = dat.removeNan2(GAMAdata['type'],GAMAdata)
  GAMAdata = dat.removeNan2(GAMAdata['z'],GAMAdata)
  GAMAdata = dat.removeNan2(GAMAdata['color'],GAMAdata)
  GAMAdata = dat.removeNan2(GAMAdata['mu'],GAMAdata)
  GAMAdata = dat.removeNan2(GAMAdata['Re'],GAMAdata)
  GAMAdata = dat.removeNan2(GAMAdata['mass'],GAMAdata)
  GAMAdata = dat.removeNan2(GAMAdata['Mr'],GAMAdata)
  clustdata = dat.removeNan2(clustdata['type'],clustdata)
  clustdata = dat.removeNan2(clustdata['z'],clustdata)
  clustdata = dat.removeNan2(clustdata['color'],clustdata)
  clustdata = dat.removeNan2(clustdata['mu'],clustdata)
  clustdata = dat.removeNan2(clustdata['Re'],clustdata)
  clustdata = dat.removeNan2(clustdata['mass'],clustdata)
  clustdata = dat.removeNan2(clustdata['Mr'],clustdata)

  # dtype_cat = [('CATID', '<f8'), ('type', '<f8'), ('z', '<f8'), ('color', '<f8'), ('ra', '<f8'), ('dec', '<f8'), ('Re', '<f8'), ('R_R200', '<f8')]
  # clustdata = np.genfromtxt("GU2_HW6_cluster.csv",dtype=dtype_cat,delimiter=',',skip_header=1,usecols=(1,2,3,4,5,6,10,11))
  # clustdata = dat.removeNan2(clustdata['type'], clustdata)
  # clustdata = dat.removeNan2(clustdata['z'],clustdata)
  # clustdata = dat.removeNan2(clustdata['color'],clustdata)
  # clustdata = dat.removeNan2(clustdata['ra'],clustdata)
  # clustdata = dat.removeNan2(clustdata['dec'],clustdata)
  # clustdata = dat.removeNan2(clustdata['Re'],clustdata)
  # clustdata = dat.removeNan2(clustdata['R_R200'],clustdata)

  def SAMI_psf(catid):
    # SAMI psf data
    dtype_psf = [('ID', '<f8'), ('Re', '<f8'), ('psf', '<f8')]
    psfdata = np.genfromtxt("psf.csv", dtype=dtype_psf, delimiter=',', skip_header=1, usecols=(1, 12, 17))
    # Data processing
    # NaN
    Renanind = np.isnan(psfdata['Re'])
    psfdata = psfdata[~Renanind]
    psfnanind = np.isnan(psfdata['psf'])
    psfdata = psfdata[~psfnanind]
    # Re > 1/2 psf
    effective_ind = psfdata['Re'] > (1 / 2) * psfdata['psf']
    effective_ID = psfdata[effective_ind]['ID']
     # Data loading
    mypath = './sami/dr3/ifs/'
    vel = '_A_stellar-velocity_default_two-moment.fits'
    disp = '_A_stellar-velocity-dispersion_default_two-moment.fits'
    veldata = fits.open(mypath+'{}'.format(catid)+'/'+'{}'.format(catid)+vel)
    dispdata = fits.open(mypath+'{}'.format(catid)+'/'+'{}'.format(catid)+disp)

    velocity = veldata[0].data
    vel_error = veldata[1].data
    vel_snr = veldata[4].data
    dispersion = dispdata[0].data
    disp_error = dispdata[1].data
    disp_snr = dispdata[4].data
     # Data processing
      # SNR
    vel_snr_cri = vel_snr > 3
    disp_snr_cri = disp_snr > 3
      # error
    vel_error_cri = vel_error < 30
    disp_error_cri = disp_error < (dispersion*0.1+25)
      # disp
    dispersion_cri = dispersion > 35

    velocity = velocity[vel_error_cri * vel_snr_cri]
    vel_error = vel_error[vel_error_cri * vel_snr_cri]
    vel_snr = vel_snr[vel_error_cri * vel_snr_cri]
    dispersion = dispersion[dispersion_cri * disp_error_cri * disp_snr_cri]
    disp_error = disp_error[dispersion_cri * disp_error_cri * disp_snr_cri]
    disp_snr = disp_snr[dispersion_cri * disp_error_cri * disp_snr_cri]
      # NaN
    velocity = dat.removeNan(velocity)
    vel_error = dat.removeNan(vel_error)
    vel_snr = dat.removeNan(vel_snr)
    dispersion = dat.removeNan(dispersion)
    disp_error = dat.removeNan(disp_error)
    disp_snr = dat.removeNan(disp_snr)
    return velocity,vel_error,vel_snr,dispersion,disp_error,disp_snr
