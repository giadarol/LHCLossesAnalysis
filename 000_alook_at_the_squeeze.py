import sys, os
BIN = os.path.expanduser("../LHC_fullRun2_analysis_scripts/")
sys.path.append(BIN)

import LHCMeasurementTools.mystyle as ms

import pickle
import numpy as np

from matplotlib.ticker import MaxNLocator

filln = 6639
filln = 6060

import pytimber
ldb = pytimber.LoggingDB(source='ldb')

# fills_pkl_name = '../LHC_2018_followup/fills_and_bmodes.pkl'
# with open(fills_pkl_name, 'rb') as fid:
#     dict_fill_bmodes = pickle.load(fid)
# t_start = dict_fill_bmodes[filln]['t_start_FLATTOP']
# t_stop = dict_fill_bmodes[filln]['t_stop_SQUEEZE']

fillinfo = ldb.getLHCFillData(filln)
bmodes = fillinfo['beamModes']
squeeze_found = False
for bm in bmodes:
    if bm['mode'] == 'FLATTOP':
        squeeze_found = True
        t_start = bm['startTime']
    if bm['mode'] == 'STABLE':
        squeeze_found = True
        t_stop = bm['startTime']+10*60.
        

data = {}

data.update(ldb.get([
			'LHC.BCTFR.A6R4.B1:BUNCH_INTENSITY',
				'LHC.BCTFR.A6R4.B2:BUNCH_INTENSITY'], t_start, t_stop))
print 'Downloaded Intensity'

data.update(ldb.get([
			'HX:BETASTAR_IP1',
            'CMS:LUMI_TOT_INST',
            'ATLAS:LUMI_TOT_INST',
				], t_start, t_stop))
print 'Downloaded Beta'

import matplotlib.pyplot as plt
import matplotlib.cm as cm
plt.close('all')
ms.mystyle_arial(fontsz=14, dist_tick_lab=5)

axt = None
axbun = None
axbet = None
for beam in [1,2]:

    bint = data['LHC.BCTFR.A6R4.B%d:BUNCH_INTENSITY'%beam][1]
    t_stamps = data['LHC.BCTFR.A6R4.B%d:BUNCH_INTENSITY'%beam][0]

    bint[bint<0.5e11]=1.

    slots = np.arange(len(bint[0, :]))

    bint_norm = bint*0

    for ii in xrange(len(t_stamps)):
    	bint_norm[ii, :] = (1.-bint[ii, :]/bint[0, :])*100.

    t_mat = np.dot(np.ones((len(slots), 1)), np.atleast_2d(t_stamps)).T
    lifet_h = -1/(np.diff(bint, axis=0)/np.diff(t_mat, axis=0)/bint[:-1,:])/3600.
    lifet_h[bint[:-1,:]<0.8e11] = 0.
    lifet_h[lifet_h<0]= 100

    t_ref = t_stamps[0]
    t_fbct_minutes = (t_stamps-t_ref)/60.
    

    fig = plt.figure(beam, figsize=(8*1.5, 6))
    fig.set_facecolor('w')
    sploss = plt.subplot2grid(shape=(1, 4), loc=(0, 1),colspan=3, sharey=axt, sharex = axbun)
    axt = sploss
    axbun = sploss
    plt.pcolormesh(slots, t_fbct_minutes, bint_norm)
    plt.colorbar()

    t_temp = data['HX:BETASTAR_IP1'][0]
    t_betast_minutes = (t_temp-t_ref)/60.

    spbet = plt.subplot2grid(shape=(1, 4), loc=(0, 0),colspan=1, sharey=axt, sharex=axbet)
    axbet = spbet
    plt.plot(data['HX:BETASTAR_IP1'][1], t_betast_minutes)
    spbet.set_ylim(0.25, 110.)
    splumi = spbet.twiny()
    splumi.xaxis.set_major_locator(MaxNLocator(5))
    var = 'ATLAS:LUMI_TOT_INST'; splumi.plot(data[var][1]/1e4, (data[var][0]-t_ref)/60., color='r')
    plt.grid('on')

    figl = plt.figure(100+beam, figsize=(8*1.5, 6))
    figl.set_facecolor('w')
    splifet = plt.subplot2grid(shape=(1, 4), loc=(0, 1),colspan=3, sharey=axt, sharex = axbun)
    plt.pcolormesh(slots, t_fbct_minutes[1:], lifet_h, vmin=0, vmax = 60, cmap=cm.jet_r)
    splifet.set_xlim(0, 3500)
    plt.colorbar()

    spbet = plt.subplot2grid(shape=(1, 4), loc=(0, 0),colspan=1, sharey=axt, sharex=axbet)
    axbet = spbet
    plt.plot(data['HX:BETASTAR_IP1'][1], t_betast_minutes)
    spbet.set_xlim(25., 110.)
    spbet.xaxis.set_major_locator(MaxNLocator(5))

    splumi = spbet.twiny()
    var = 'ATLAS:LUMI_TOT_INST'; splumi.plot(data[var][1]/1e4, (data[var][0]-t_ref)/60., color='r')
    splumi.xaxis.set_major_locator(MaxNLocator(5))
    plt.grid('on')
    splumi.set_xlabel('Lumi')
    splumi.set_ylim(bottom=(t_start-t_ref)/60.,top=(t_stop-t_ref)/60.)
    splumi.xaxis.label.set_color('r')
    splumi.tick_params(axis='x', colors='red')

    figl.suptitle('Fill %d Beam %d'%(filln, beam))

    figl.subplots_adjust(left=.05, right=1.)

plt.figure(3)
spl = plt.subplot(2,1,1)
plt.plot(bint[-1,:])
plt.subplot(2,1,2, sharex=spl)
plt.plot(lifet_h[-1,:])

plt.show()