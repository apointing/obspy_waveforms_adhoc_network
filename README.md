# obspy_waveforms_adhoc_network
Use obspy to plot waveforms from stations downloaded from adhoc networks (including raspberry shakes)

Use obspy to plot waveforms for stations downloaded from swarm or other adhoc networks and overlay with calculated travel times

You need the obspy python package in order to execute this python script, you can download it from obspy.org

This script has been tested on Apple Mac

From a terminal command line, type:

python plot_section_adhoc_network.py

For example, if you have found an event that is picked up by a number of raspberry shakes, then export the waveforms as mseed files 

from Swarm to a local directory. You can also export mseed waveforms from other adhoc networks e.g. in Europe the ORFEUS network, which 

allows you to download waveforms. Collect the waveforms into a directory. 

Then edit the python script and fill in the user data and parameters shown in the script. There is example data filled in, in this

script for an event in Norway. The example data is uploaded to the repository here, download it and rename the path in the script.

There is an issue at the moment with the plotting routine that doesn't scale the HHZ channels for stations, for some reason.
