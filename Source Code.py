
# Import numpy and matplotlib. For the latter, use a nicer set of plot
# parameters and set up support for plotting/converting quantities.
import numpy as np
import matplotlib.pyplot as mpl
from astropy.visualization import astropy_mpl_style, quantity_support
mpl.style.use(astropy_mpl_style)
quantity_support()

# Import the packages necessary for finding coordinates and making
# coordinate transformations
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

# Get the coordinates of deep sky object:
# In this case, Orion Nebula (M42)
m42 = SkyCoord.from_name('M42')

# Use `astropy.coordinates.EarthLocation` to provide the location of the
# observer and set the time to 11pm GMT+8 on 11 February 2021:
# For example, the location is in Bagan Datoh
lat = 3.0 
lon = 100.0
Bagan_Datuk = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=4*u.m)
utcoffset = +8*u.hour  # UTC +8
time = Time('2021-2-11 23:00:00') - utcoffset

# Use `astropy.coordinates` to find the Alt, Az coordinates of M42 at as
# observed from Bagan Datuk at 11pm on 11 February 2021.
m42altaz = m42.transform_to(AltAz(obstime=time,location=Bagan_Datuk))
print("-----------------------------------------")
print("The Altitude/Azimuth Position for M42")
print("The Coordinate of the Observer: ("+str(lat)+","+str(lon)+") ")
print("The Time and Date: ",time)
print("-----------------------------------------")
print("M42 Altitude = {0.alt:.10}".format(m42altaz))
print("M42 Azimuth = {0.az:.10}".format(m42altaz))

midnight = Time('2021-2-11 00:00:00') - utcoffset
delta_midnight = np.linspace(-2, 10, 100)*u.hour
frame_Feb11Night = AltAz(obstime=midnight+delta_midnight,
                          location=Bagan_Datuk)
m42altazs_Feb11Night = m42.transform_to(frame_Feb11Night)

# Use  `~astropy.coordinates.get_sun` to find the location of the Sun at 1000
# evenly spaced times between noon on 11 February and noon on 12 February:
from astropy.coordinates import get_sun
delta_midnight = np.linspace(-12, 12, 1000)*u.hour
times_Feb11_to_12 = midnight + delta_midnight
frame_Feb11_to_12 = AltAz(obstime=times_Feb11_to_12, location=Bagan_Datuk)
sunaltazs_Feb11_to_12 = get_sun(times_Feb11_to_12).transform_to(frame_Feb11_to_12)

# Do the same with `~astropy.coordinates.get_moon` to find when the moon is.
from astropy.coordinates import get_moon
moon_Feb11_to_12 = get_moon(times_Feb11_to_12)
moonaltazs_Feb11_to_12 = moon_Feb11_to_12.transform_to(frame_Feb11_to_12)

# Find the alt,az coordinates of M42 at those same times:
m42altazs_Feb11_to_12 = m42.transform_to(frame_Feb11_to_12)

#Plotting
mpl.plot(delta_midnight, sunaltazs_Feb11_to_12.alt, color='r', label='Sun')
mpl.plot(delta_midnight, moonaltazs_Feb11_to_12.alt, color=[0.75]*3, ls='--', label='Moon')
mpl.scatter(delta_midnight, m42altazs_Feb11_to_12.alt,
            c=m42altazs_Feb11_to_12.az, label='M42 Orion Nebula', lw=0, s=8,
            cmap='viridis')
mpl.fill_between(delta_midnight, 0*u.deg, 90*u.deg,
                 sunaltazs_Feb11_to_12.alt < -0*u.deg, color='0.5', zorder=0)
mpl.fill_between(delta_midnight, 0*u.deg, 90*u.deg,
                 sunaltazs_Feb11_to_12.alt < -18*u.deg, color='k', zorder=0)
mpl.colorbar().set_label('Azimuth [deg]')
mpl.legend(loc='upper left')
mpl.xlim(-12*u.hour, 12*u.hour)
mpl.xticks((np.arange(13)*2-12)*u.hour)
mpl.ylim(0*u.deg, 90*u.deg)
mpl.xlabel('Hours from midnight at given coordinate')
mpl.ylabel('Altitude [deg]')
mpl.title("Altitude and Azimuth position of M42")
mpl.show()