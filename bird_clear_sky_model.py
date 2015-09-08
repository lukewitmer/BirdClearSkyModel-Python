from __future__ import division  # ensures no rounding errors from division involving integers

__author__ = 'Luke Witmer'

from math import * # enables use of pi, trig functions, and more.
import pandas as pd # gives us the dataframe concept
pd.options.display.max_columns = 50
pd.options.display.max_rows = 9

### User Inputs

phi = 40.72
longitude = -77.93
tz = -5
P_mb = 970
Ozone_cm = 0.3
H20_cm = 1.5
AOD500nm = 0.193
AOD380nm = 0.298
Taua = 0.15
Ba = 0.85
albedo = 0.2

G_sc = 1367 # W/m^2
std_mer = longitude-longitude%15+15  # This Standard Meridian calculation is only a guide!!
                                     # Please double check this value for your location!

### Day of the Year Column

n = range(1, 366)  # julian day of the year
n_hrly = list(pd.Series(n).repeat(24))  # julian day numbers repeated hourly to create 8760 datapoints in dataset

ds = pd.DataFrame(n_hrly, columns=['DOY'])  # create dataframe with julian days

### Hr of the Day Column

ds['HR'] = [(hr)%24 for hr in ds.index.tolist()]  # append dataframe with hr of the day for each day

### Extraterrestrial Radiation

def etr(n):
    return G_sc*(1.00011+0.034221*cos(2*pi*(n-1)/365)+0.00128*sin(2*pi*(n-1)/365)+0.000719*cos(2*(2*pi*(n-1)/365))+0.000077*sin(2*(2*pi*(n-1)/365)))

ds['ETR'] = [etr(n) for n in ds['DOY']] # append dataframe with etr for day

### Intermediate Parameters

ds['Dangle'] = [2*pi*(n-1)/365 for n in ds['DOY']]

def decl(Dangle):
    return (0.006918-0.399912*cos(Dangle)+0.070257*sin(Dangle)-0.006758*cos(2*Dangle)+0.000907*sin(2*Dangle)-0.002697*cos(3*Dangle)+0.00148*sin(3*Dangle))*(180/pi)
ds['DEC'] = [decl(Dangle) for Dangle in ds['Dangle']]

def eqtime(Dangle):
    return (0.0000075+0.001868*cos(Dangle)-0.032077*sin(Dangle)-0.014615*cos(2*Dangle)-0.040849*sin(2*Dangle))*229.18
ds['EQT'] = [eqtime(Dangle) for Dangle in ds['Dangle']]

def omega(hr, eqt):
    return 15*(hr-12.5) + longitude - tz*15 + eqt/4
ds['Hour Angle'] = [omega(hr, eqt) for hr, eqt in zip(ds['HR'],ds['EQT'])]

def zen(dec, hr_ang):
    return acos(cos(dec/(180/pi))*cos(phi/(180/pi))*cos(hr_ang/(180/pi))+sin(dec/(180/pi))*sin(phi/(180/pi)))*(180/pi)
ds['Zenith Ang'] = [zen(dec, hr_ang) for dec, hr_ang in zip(ds['DEC'],ds['Hour Angle'])]

def airmass(zenang):
    if zenang < 89:
        return 1/(cos(zenang/(180/pi))+0.15/(93.885-zenang)**1.25)
    else:
        return 0
ds['Air Mass'] = [airmass(zenang) for zenang in ds['Zenith Ang']]

### Intermediate Results

def T_rayleigh(airmass):
    if airmass > 0:
        return exp(-0.0903*(P_mb*airmass/1013)**0.84*(1+P_mb*airmass/1013-(P_mb*airmass/1013)**1.01))
    else:
        return 0
ds['T rayleigh'] = [T_rayleigh(airmass) for airmass in ds['Air Mass']]

def T_ozone(airmass):
    if airmass > 0:
        return 1-0.1611*(Ozone_cm*airmass)*(1+139.48*(Ozone_cm*airmass))**-0.3034-0.002715*(Ozone_cm*airmass)/(1+0.044*(Ozone_cm*airmass)+0.0003*(Ozone_cm*airmass)**2)
    else:
        return 0
ds['T ozone'] = [T_ozone(airmass) for airmass in ds['Air Mass']]

def T_gasses(airmass):
    if airmass > 0:
        return exp(-0.0127*(airmass*P_mb/1013)**0.26)
    else:
        return 0
ds['T gases'] = [T_gasses(airmass) for airmass in ds['Air Mass']]

def T_water(airmass):
    if airmass > 0:
        return 1-2.4959*airmass*H20_cm/((1+79.034*H20_cm*airmass)**0.6828+6.385*H20_cm*airmass)
    else:
        return 0
ds['T water'] = [T_water(airmass) for airmass in ds['Air Mass']]

def T_aerosol(airmass):
    if airmass > 0:
        return exp(-(Taua**0.873)*(1+Taua-Taua**0.7088)*airmass**0.9108)
    else:
        return 0
ds['T aerosol'] = [T_aerosol(airmass) for airmass in ds['Air Mass']]

def taa(airmass, taerosol):
    if airmass > 0:
        return 1-0.1*(1-airmass+airmass**1.06)*(1-taerosol)
    else:
        return 0
ds['TAA'] = [taa(airmass, taerosol) for airmass, taerosol in zip(ds['Air Mass'],ds['T aerosol'])]

def rs(airmass, taerosol, taa):
    if airmass > 0:
        return 0.0685+(1-Ba)*(1-taerosol/taa)
    else:
        return 0
ds['rs'] = [rs(airmass, taerosol, taa) for airmass, taerosol, taa in zip(ds['Air Mass'],ds['T aerosol'],ds['TAA'])]

def Id(airmass, etr, taerosol, twater, tgases, tozone, trayleigh):
    if airmass > 0:
        return 0.9662*etr*taerosol*twater*tgases*tozone*trayleigh
    else:
        return 0
ds['Id'] = [Id(airmass, etr, taerosol, twater, tgases, tozone, trayleigh) for airmass, etr, taerosol, twater, tgases, tozone, trayleigh in zip(ds['Air Mass'],ds['ETR'],ds['T aerosol'],ds['T water'],ds['T gases'],ds['T ozone'],ds['T rayleigh'])]

def idnh(zenang, Id):
    if zenang < 90:
        return Id*cos(zenang/(180/pi))
    else:
        return 0
ds['IdnH'] = [idnh(zenang, Id) for zenang, Id in zip(ds['Zenith Ang'],ds['Id'])]

def ias(airmass, etr, zenang, tozone, tgases, twater, taa, trayleigh, taerosol):
    if airmass > 0:
        return etr*cos(zenang/(180/pi))*0.79*tozone*tgases*twater*taa*(0.5*(1-trayleigh)+Ba*(1-(taerosol/taa)))/(1-airmass+(airmass)**1.02)
    else:
        return 0
ds['Ias'] = [ias(airmass, etr, zenang, tozone, tgases, twater, taa, trayleigh, taerosol) for airmass, etr, zenang, tozone, tgases, twater, taa, trayleigh, taerosol in zip(ds['Air Mass'],ds['ETR'],ds['Zenith Ang'],ds['T ozone'],ds['T gases'],ds['T water'],ds['TAA'],ds['T rayleigh'],ds['T aerosol'])]

def gh(airmass, idnh, ias, rs):
    if airmass > 0:
        return (idnh+ias)/(1-albedo*rs)
    else:
        return 0
ds['GH'] = [gh(airmass, idnh, ias, rs) for airmass, idnh, ias, rs in zip(ds['Air Mass'],ds['IdnH'],ds['Ias'],ds['rs'])]

### Decimal Time

def dectime(doy, hr):
    return doy+(hr-0.5)/24
ds['Decimal Time'] = [dectime(doy, hr) for doy, hr in zip(ds['DOY'],ds['HR'])]


### Model Results (W/m^2)

ds['Direct Beam'] = ds['Id']

ds['Direct Hz'] = ds['IdnH']

ds['Global Hz'] = ds['GH']

ds['Dif Hz'] = ds['Global Hz']-ds['Direct Hz']

# ds

# get_ipython().magic(u'matplotlib inline')
# import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
pylab.rcParams['figure.figsize'] = 16, 6  # this sets the default image size for this session

ax = ds[ds['DOY']==212].plot('HR',['Global Hz','Direct Hz','Dif Hz'],title='Bird Clear Sky Model Results')
ax.set_xlabel("Hour")
ax.set_ylabel("Irradiance W/m^2")
majorx = ax.set_xticks(range(0,25,1))
majory = ax.set_yticks(range(0,1001,200))


### SURFRAD Data

initial_cols_list = list(pd.read_table('surfrad_col_headers.txt', sep="\t", skipinitialspace=True, header=None)[0])
qc_flag_cols = ['qc'+str(n) for n in range(1,21)]
col_headers = initial_cols_list[:8] + [item for sublist in zip(initial_cols_list[8:],qc_flag_cols) for item in sublist]

surfrad = pd.read_table('psu07212.dat', skiprows=[0,1], sep=" ", skipinitialspace=True, header=None, names=col_headers)

### Calculate local time and add to surfrad dataframe

surfrad['jday_dt'] = [jday+dt/24 for jday, dt in zip(surfrad['jday'],surfrad['dt'])]

surfrad['Dangle'] = [2*pi*(n-1)/365 for n in surfrad['jday_dt']]

surfrad['DEC'] = [decl(Dangle) for Dangle in surfrad['Dangle']]

surfrad['EQT'] = [eqtime(Dangle) for Dangle in surfrad['Dangle']]

# surfrad['t_local'] = [(dt+tz)+(4*(longitude-std_mer)+eqt)/60 for dt,eqt in zip(surfrad['dt'],surfrad['EQT'])]
surfrad['t_local'] = [(dt+tz+1)+(4*(longitude-std_mer)+eqt)/60 for dt,eqt in zip(surfrad['dt'],surfrad['EQT'])]
# surfrad['t_local'] = [(dt+tz+1) for dt in surfrad['dt']]

# surfrad

ay = surfrad.plot('t_local',['dw_solar','netsolar','diffuse'],title='SURFRAD Data - Penn State')
ay.set_xlabel("Hour")
ay.set_ylabel("Irradiance W/m^2")
majorx = ay.set_xticks(range(0,25,1))
ay.set_xlim([0,24])
ay.set_ylim(0)
majory = ay.set_yticks(range(0,1001,200))


### Combine the two plots above

ax = ds[ds['DOY']==212].plot('HR',['Global Hz','Direct Hz','Dif Hz'],title='Bird Clear Sky Model with SURFRAD Data')
ay = surfrad.plot('t_local',['dw_solar','netsolar','diffuse'], ax=ax)
ay.set_xlabel("Hour")
ay.set_ylabel("Irradiance W/m^2")
majorx = ay.set_xticks(range(0,25,1))
ay.set_xlim([0,24])
ay.set_ylim(0)
majory = ay.set_yticks(range(0,1001,200))
