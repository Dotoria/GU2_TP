# %%
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15 as cosmo
import astropy.units as u
from tqdm import tqdm

# %% [markdown]
# specObjID,z,z_err,ra,dec,petroMag_g,petroMag_r,petroMag_i,h_alpha_flux,h_alpha_flux_err,h_beta_flux,h_beta_flux_err,oiii_5007_flux,oiii_5007_flux_err,oi_6300_flux,oi_6300_flux_err,nii_6584_flux,nii_6584_flux_err,sii_6717_flux,sii_6717_flux_err,sii_6731_flux,sii_6731_flux_err,h_delta_flux,h_delta_flux_err

# %%
data=np.loadtxt('GU2_TP_third_BeomChan_Koh.csv',skiprows=1,delimiter=',')

# %%
data=data.T

# %%
ID=data[0]
z=data[1]
z_err=data[2]
RA=data[3]
DEC=data[4]
gmag=data[5]
rmag=data[6]
imag=data[7]
halpha=data[8]
halpha_err=data[9]
hbeta=data[10]
hbeta_err=data[11]
oiii=data[12]
oiii_errr=data[13]
oi=data[14]
oi_err=data[15]
nii=data[16]
nii_err=data[17]
sii6717=data[18]
sii6717_err=data[19]
sii6731=data[20]
sii6731_err=data[21]
hdelta=data[22]
hdelta_err=data[23]


# %%
Mr= rmag - 5.0*np.log10(cosmo.luminosity_distance(z).to('pc').value/10.0)

# %%
rcrit=17.77 -5.0*np.log10(cosmo.luminosity_distance(0.17).value)-25


# %%
volsam=np.where((Mr<rcrit)&(rmag<17.7)&(rmag>14.5))[0]
volsamemgal=np.where((Mr<rcrit)&(rmag<17.7)&(rmag>14.5)&(halpha>0)&(hbeta>0)&(oiii>0)&(nii>0))[0]
len(volsam),len(volsamemgal)

# %%
niiHal=np.log10(nii[volsamemgal]/halpha[volsamemgal])
oiiiHbeta=np.log10(oiii[volsamemgal]/hbeta[volsamemgal])

sigma5=np.zeros(len(volsam))
for i, idx in enumerate(tqdm(volsam)):
    rai=RA[idx]
    deci=DEC[idx]
    onedeg=(cosmo.angular_diameter_distance(z[idx])/u.rad).to('Mpc/degree')
    samp2=np.where((RA<rai+2)&(RA>rai-2)&(DEC<deci+2)&(DEC>deci-2)&(z<z[idx]+0.003)&(z<z[idx]-0.003))[0]
#    print(len(samp2))
    if len(samp2)<5:
        sigma5[i]=0.000000001
    else:
        distance=np.zeros(len(samp2))
        for j, idx2 in enumerate(samp2):
            raj=RA[idx2]
            decj=DEC[idx2]
            radist=(abs(rai-raj)*u.degree*onedeg).value
            decdist=(abs(deci-decj)*u.degree*onedeg).value
            distance[j]=np.sqrt(radist**2+decdist**2)
        top5dist=np.sort(distance)[4]
        sigma5[i]=5/(np.pi*top5dist**2)
np.save('sigma5',sigma5)

