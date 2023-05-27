import numpy as np
import matplotlib.pyplot as plt
from data_processing import data_processing as dat
from astropy.cosmology import WMAP9 as cosmo
from astropy import constants as const

dtype = [('id' ,'<f8') ,('z' ,'<f8') ,('z_err' ,'<f8'),('ra' ,'<f8') ,('dec' ,'<f8') ,('g' ,'<f8') ,('r' ,'<f8')
             ,('i' ,'<f8') ,('Ha_flux' ,'<f8') ,('Ha_err' ,'<f8') ,('Hb_flux' ,'<f8'),('Hb_err' ,'<f8'),('O3_flux' ,'<f8'),('O3_err' ,'<f8'),('N2_flux' ,'<f8'),('N2_err' ,'<f8')]
data0 = np.genfromtxt("GU2_TP_Test2.csv" ,dtype=dtype ,delimiter=',' ,skip_header=1
                         ,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,16,17))

# mag cutting
cri_mr = data0['r'] < 17.77
cri_mr2 = data0['r'] > 14.5 # Kauffmann et al. 2003
# cri_z1 = data0['z'] < 0.1
# cri_z2 = data0['z'] > 0.08
# cri_snr_Ha = data0['Ha_err'] < 3
# cri_snr_Hb = data0['Hb_err'] < 3
# cri_snr_O3 = data0['O3_err'] < 3
# cri_snr_N2 = data0['N2_err'] < 3
z0 = data0['z']
dL0 = cosmo.luminosity_distance(z0).value
Mr0 = data0['r'] - 5 * np.log10(dL0) + 25
Mr_limit = 17.77 - 5 * np.log10(max(dL0)) + 25
cri_Mr = Mr0 < Mr_limit
cri_dL = dL0 > 0 # Volume limited sampling
cri = cri_mr & cri_Mr & cri_dL & cri_mr2
data = data0[cri]
#print('Original data :',data0.shape,'Sampled data :', data.shape)

# BPT
Ha_list = data['Ha_flux']
Hb_list = data['Hb_flux']
O3_list = data['O3_flux']
N2_list = data['N2_flux']
"""
plt.scatter(np.log10(N2_list/Ha_list),np.log(O3_list/Hb_list),s=0.2,c='black',alpha=0.5)
plt.axhline(np.log10(3))
plt.axvline(np.log10(0.6))
x = np.linspace(-1.5,-0.1,100)
y = 0.61 / (x-0.05) + 1.3
plt.plot(x,y,'r--')
plt.xlim(-1.5,1)
plt.ylim(-1.5,1.5)
plt.xlabel(r'[NII] / $h_{\alpha}$')
plt.ylabel(r'[OIII] / $h_{\beta}$')
plt.show()
"""
# Classifying
AGN = []
LINER = []
Seyfert = []
for i in range(len(Ha_list)):
    """
    Append object id into AGN, LINER, Seyfert lists
    """
    Ha = Ha_list[i]
    Hb = Hb_list[i]
    O3 = O3_list[i]
    N2 = N2_list[i]
    if (O3 / Hb) < 3 and (N2 / Ha) > 0.6 and np.log10(O3 / Hb) > (0.61 / (np.log10(N2 / Ha) - 0.05) + 1.3):
        LINER.append(int(data['id'][i]))
    elif (O3 / Hb) > 3 and (N2 / Ha) > 0.6 and np.log10(O3 / Hb) > (0.61 / (np.log10(N2 / Ha) - 0.05) + 1.3):
        Seyfert.append(int(data['id'][i]))
    elif np.log10(O3 / Hb) > (0.61 / (np.log10(N2 / Ha) - 0.05) + 1.3):
        AGN.append(int(data['id'][i]))
"""Use list name of AGN, Seyfert and LINER"""

# Weak & Strong
"""
For L[O3] > 1e7 * L_solar, AGN fraction no longer depends on z/z_max : Strong AGNs
"""
Strong_Sey = []
Weak_Sey = []
#print(const.L_sun.value) # 3.828e+26 W
"""
sdss flux의 단위 1e-17 erg/s/cm^2
"""
solar_lumi_erg = const.L_sun.value * 1e7
for id in Seyfert:
    O3_lumi = 4 * np.pi * ((cosmo.luminosity_distance(data[data['id']==id]['z']).value * 3.086e+22) ** 2) * data[data['id']==id]['O3_flux'] * 1e-13
    if np.log10(O3_lumi/solar_lumi_erg) > 7:
        Strong_Sey.append(id)
    else : Weak_Sey.append(id)
Strong_L = []
Weak_L = []
for id in LINER:
    O3_lumi = 4 * np.pi * (cosmo.luminosity_distance(data[data['id']==id]['z']).value * 3.086e+22) ** 2 * data[data['id']==id]['O3_flux'] * 1e-13
    if np.log10(O3_lumi/solar_lumi_erg) > 7:
        Strong_L.append(id)
    else : Weak_L.append(id)
Strong_AGN = []
Weak_AGN = []
for id in AGN:
    O3_lumi = 4 * np.pi * (cosmo.luminosity_distance(data[data['id']==id]['z']).value * 3.086e+22) ** 2 * data[data['id']==id]['O3_flux'] * 1e-13
    if np.log10(O3_lumi/solar_lumi_erg) > 7:
        Strong_AGN.append(id)
    else : Weak_AGN.append(id)
print('Strong :',len(Strong_AGN)+len(Strong_L)+len(Strong_Sey))
print('Weak : ',len(Weak_AGN)+len(Weak_L)+len(Weak_Sey))

"""
추가할것 : 
- A, L, S 색 표시한 BPT
- logL[O3] 범위별로 분리한 BPT
"""
class classifyAGN :
    AGN = AGN
    Seyfert = Seyfert
    LINER = LINER
    Strong_Sey = Strong_Sey
    Weak_Sey = Weak_Sey
    Strong_L = Strong_L
    Weak_L = Weak_L
    Strong_AGN = Strong_AGN
    Weak_AGN = Weak_AGN