import math
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP9 as cosmo
from data_processing import data_processing as dat
from SAMI import SAMI
import astropy.io.fits as fits
import spectral as sp
import spectral.io.aviris as aviris
from PyAstronomy import pyasl

GAMAid = list(SAMI.GAMAdata['CATID'])
clustid = list(SAMI.clustdata['CATID'])
CATID = GAMAid + clustid
CATID.sort()
#print(len(CATID))

def SAMI_spectra(catid):
    mypath = './sami/dr3/ifs/'
    vel = '_A_stellar-velocity_default_two-moment.fits'
    disp = '_A_stellar-velocity-dispersion_default_two-moment.fits'
    Halpha = '_A_Halpha_default_recom-comp.fits'
    Hbeta = '_A_Hbeta_default_recom-comp.fits'
    N2 = '_A_NII6583_default_recom-comp.fits'
    O1 = '_A_OI6300_default_recom-comp.fits'
    O2 = '_A_OII3728_default_recom-comp.fits'
    O3 = '_A_OIII5007_default_recom-comp.fits'
    S6716 = '_A_SII6716_default_recom-comp.fits'
    S6731 = '_A_SII6731_default_recom-comp.fits'
    #veldata = fits.open(mypath+'{}'.format(catid)+'/'+'{}'.format(catid)+vel)
    #dispdata = fits.open(mypath+'{}'.format(catid)+'/'+'{}'.format(catid)+disp)
    Hadata = fits.open(mypath+'{}'.format(catid)+'/'+'{}'.format(catid) + Halpha)[0].data
    Hbdata = fits.open(mypath+'{}'.format(catid)+'/'+'{}'.format(catid) + Hbeta)[0].data
    N2data = fits.open(mypath+'{}'.format(catid)+'/'+'{}'.format(catid) + N2)[0].data
    O1data = fits.open(mypath+'{}'.format(catid)+'/'+'{}'.format(catid) + O1)[0].data
    O2data = fits.open(mypath+'{}'.format(catid)+'/'+'{}'.format(catid) + O2)[0].data
    O3data = fits.open(mypath+'{}'.format(catid)+'/'+'{}'.format(catid) + O3)[0].data
    S6731data = fits.open(mypath+'{}'.format(catid)+'/'+'{}'.format(catid) + S6731)[0].data
    S6716data = fits.open(mypath+'{}'.format(catid)+'/'+'{}'.format(catid) + S6716)[0].data
    #print(Hadata.info(),Hbdata.info(),N2data.info(),O1data.info(),O2data.info(),O3data.info(),S6716data.info(),S6731data.info())
    Hadata = dat.removeNan(Hadata)
    Hbdata = dat.removeNan(Hbdata)
    N2data = dat.removeNan(N2data)
    O1data = dat.removeNan(O1data)
    O2data = dat.removeNan(O2data)
    O3data = dat.removeNan(O3data)
    S6731data = dat.removeNan(S6731data)
    S6716data = dat.removeNan(S6716data)
    return Hadata, Hbdata, N2data, O1data, O2data, O3data, S6731data, S6716data
HaList = []
HbList = []
O3List = []
N2List = []
LINER = []
Seyfert = []
AGN = []
for id in CATID:
    try :
        Hadata, Hbdata, N2data, O1data, O2data, O3data, S6731data, S6716data = SAMI_spectra(int(id))

        Ha = np.nanmax(Hadata)
        Hb = np.nanmax(Hbdata)
        O3 = np.nanmax(O3data)
        N2 = np.nanmax(N2data)

        HaList.append(Ha)
        HbList.append(Hb)
        O3List.append(O3)
        N2List.append(N2)
        if np.log10(O3/Hb) > 0.61/(np.log10(N2/Ha)-0.05) + 1.3:

            AGN.append(int(id))
        elif O3/Hb < 3 and N2/Ha > 0.6 :
            LINER.append(int(id))
        elif O3/Hb > 3 and N2/Ha > 0.6 :
            Seyfert.append(int(id))
    except(ValueError) : continue

def view(catid):
        mypath = './sami/dr3/ifs/'
        Halpha = '_A_Halpha_default_recom-comp.fits'
        Hbeta = '_A_Hbeta_default_recom-comp.fits'
        N2 = '_A_NII6583_default_recom-comp.fits'
        O1 = '_A_OI6300_default_recom-comp.fits'
        O2 = '_A_OII3728_default_recom-comp.fits'
        O3 = '_A_OIII5007_default_recom-comp.fits'
        S6716 = '_A_SII6716_default_recom-comp.fits'
        S6731 = '_A_SII6731_default_recom-comp.fits'
        Hadata = fits.open(mypath + '{}'.format(catid) + '/' + '{}'.format(catid) + Halpha)[0].data
        Hbdata = fits.open(mypath + '{}'.format(catid) + '/' + '{}'.format(catid) + Hbeta)[0].data
        N2data = fits.open(mypath + '{}'.format(catid) + '/' + '{}'.format(catid) + N2)[0].data
        O1data = fits.open(mypath + '{}'.format(catid) + '/' + '{}'.format(catid) + O1)[0].data
        O2data = fits.open(mypath + '{}'.format(catid) + '/' + '{}'.format(catid) + O2)[0].data
        O3data = fits.open(mypath + '{}'.format(catid) + '/' + '{}'.format(catid) + O3)[0].data
        S6731data = fits.open(mypath + '{}'.format(catid) + '/' + '{}'.format(catid) + S6731)[0].data
        S6716data = fits.open(mypath + '{}'.format(catid) + '/' + '{}'.format(catid) + S6716)[0].data

        plt.subplot(3 ,3, 1)
        plt.imshow(Hadata[0])
        plt.title('Ha')

        plt.subplot(3, 3, 2)
        plt.imshow(Hbdata)
        plt.title('Hb')

        plt.subplot(3, 3, 3)
        plt.imshow(N2data)
        plt.title('N2')

        plt.subplot(3, 3, 4)
        plt.imshow(O1data)
        plt.title('O1')

        plt.subplot(3, 3, 5)
        plt.imshow(O2data)
        plt.title('O2')

        plt.subplot(3, 3, 6)
        plt.imshow(O3data)
        plt.title('O3')

        plt.subplot(3, 3, 7)
        plt.imshow(S6731data)
        plt.title('S6731')

        plt.subplot(3, 3, 8)
        plt.imshow(S6716data)
        plt.title('S6716')
        plt.legend(str(catid))
        plt.show()


print(len(AGN))

for id in AGN:
    view(id)
    # ra=SAMI.GAMAdata['ra']['CATID'==id]
    # dec=SAMI.GAMAdata['dec']['CATID' == id]
    # print(id, ra, dec)

"""plt.scatter(np.log10(np.array(N2List)/np.array(HaList)), np.log10(np.array(O3List)/np.array(HbList)),s=0.5,c='black')
plt.axhline(np.log10(3))
plt.axvline(np.log10(0.6))
plt.show()
print(AGN,len(AGN))"""
