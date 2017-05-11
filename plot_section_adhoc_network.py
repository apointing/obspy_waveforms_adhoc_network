import os
import obspy
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from matplotlib.transforms import blended_transform_factory
from obspy import read, Stream, UTCDateTime
from obspy.geodetics import gps2dist_azimuth, locations2degrees
from obspy.taup.taup import travelTimePlot
from obspy.core import AttribDict
from obspy.geodetics import kilometer2degrees
import numpy as np
import warnings
from obspy.taup.taup import getTravelTimes
from obspy.clients.fdsn import Client
from obspy.core import read
from obspy.core.trace import Trace

# User defined data and parameters
# Specify the location of the mseed files
path = '/Applications/swarm/obspy/station_location/data/'
# Earthquakes' epicenter
p_onset = obspy.UTCDateTime("2017-04-29T22:58:50.4")
eq_lat = 58.19
eq_lon = 6.81
depth=16
# define record length of waveform in seconds
recl=240
# Define the calculated travel times of phases - use ttall for all phases
phases = ['P', 'S', 'PP']
model='iasp91'

c = Client("ORFEUS")

class lat(object): pass
class lon(object): pass
a=lat()
b=lon()

# User defined location data
# Define the a-lattude and b-longitude for those stations not in the FDSN network
setattr(a, "DKMUD", 56.455)
setattr(b, "DKMUD", 9.173)
setattr(a, "AMRD4A6", 55.76)
setattr(b, "AMRD4A6", 12.4)

xaxispos = np.zeros(100)
trmax = np.zeros(100)

st = obspy.core.stream.Stream()
for filename in os.listdir(path):
    try:
        if filename.endswith(".DS_Store"):
           print("skipping")
        else:
            st += read(path + filename, starttime=p_onset, endtime=p_onset + 240)
    except Exception as e:
        raise e
        print "No files found here!"
print(st.__str__(extended=True))
# User defined selection of channels
st2 = st.select(id="*.00.BHZ")
st2 += st.select(id="*..BHZ")
st2 += st.select(id="*.00.HHZ")
st2 += st.select(id="*.SHZ")
print(st2.__str__(extended=True))
pos = 0
for tr in st2:
    try:
        inventory = c.get_stations(network=tr.stats.network, station=tr.stats.station, level="station")
        for network in inventory:
            for station in network:
                tr.stats.distance = gps2dist_azimuth(station.latitude,station.longitude,
                                                     eq_lat, eq_lon)[0]
                xaxispos[pos] = kilometer2degrees(tr.stats.distance / 1e3)
                tr.stats.coordinates = \
                   AttribDict(dict(latitude=station.latitude,
                                   longitude=station.longitude))
                pos += 1
    except Exception as e:
#        raise e
        print("station not found in FDSN, switching to user defined locations", tr.stats.station)
        latstr=''.join((tr.stats.network,tr.stats.station))
        lonstr=''.join((tr.stats.network,tr.stats.station))
        tr.stats.distance = gps2dist_azimuth(getattr(a,latstr), getattr(b,lonstr),
                                             eq_lat, eq_lon)[0]
        xaxispos[pos] = kilometer2degrees(tr.stats.distance / 1e3)
        tr.stats.coordinates = \
           AttribDict(dict(latitude=getattr(a,latstr),
                           longitude=getattr(b,lonstr)))
        pos += 1

# User defined section - use if yuo need it
#st2.filter('bandpass', freqmin=0.5, freqmax=4)

# Do the section plot..
# If no customization is done after the section plot command, figure
# initialization can be left out and also option ".., show=False, fig=fig)" can
# be omitted, and figure is shown automatically
fig = plt.figure()
st2.plot(type='section', dist_degree=True, ev_coord=(eq_lat,eq_lon), plot_dx=20e3, recordlength=recl,
        time_down=True, linewidth=.25, grid_linewidth=.25, show=False, fig=fig)

ax = fig.axes[0]
xmin, xmax = ax.get_xlim()

data = {}

min_degree=round(xmin)
max_degree=round(xmax)
npoints = max_degree - min_degree + 1

degrees = np.linspace(min_degree, max_degree, npoints)
# Loop over all degrees.
for degree in degrees:
    with warnings.catch_warnings(record=True):
        warnings.simplefilter('always')
        tt = getTravelTimes(degree, depth, model, phase_list=phases)
        # Mirror if necessary.
    if degree > 180:
        degree = 180 - (degree - 180)
    for item in tt:
        phase = item['phase_name']
        if phase not in data:
            data[phase] = [[], []]
        data[phase][1].append(item['time'] )
        data[phase][0].append(degree)
# Plot and some formatting.
for key, value in data.items():
    plt.plot(value[0], value[1], '.', label=key)
plt.grid()
plt.xlabel('Distance (degrees)')
plt.ylabel('Time (minutes)')
ax.set_xticks(degrees)
if max_degree <= 180:
    plt.xlim(min_degree, max_degree)
else:
    plt.xlim(min_degree, 180)
plt.legend(numpoints=1)

# Plot customization: Add station labels to offset axis
transform = blended_transform_factory(ax.transData, ax.transAxes)
pos = 0
for tr in st2:
    ax.text(xaxispos[pos], 1.0, tr.stats.station, rotation=270,
            va="bottom", ha="center", transform=transform, 
            zorder=10, fontsize=8)
    pos += 1
ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
plt.tight_layout(pad=3)
plt.show()
