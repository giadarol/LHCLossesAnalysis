import sys, os
BIN = os.path.expanduser("../LHC_fullRun2_analysis_scripts/")
sys.path.append(BIN)

import LHCMeasurementTools.mystyle as ms

import pickle
import numpy as np

from matplotlib.ticker import MaxNLocator

filln = 6639
filln = 6646
# filln = 6060
filln = 6645
filln=6640
filln=6641
filln=6650
filln=6662
filln = 6666
filln = 6674 
filln = 6672
filln = 6860

import pytimber
ldb = pytimber.LoggingDB(source='ldb')
mdb = pytimber.LoggingDB(source='mdb')
# fills_pkl_name = '../LHC_2018_followup/fills_and_bmodes.pkl'
# with open(fills_pkl_name, 'rb') as fid:
#     dict_fill_bmodes = pickle.load(fid)
# t_start = dict_fill_bmodes[filln]['t_start_FLATTOP']
# t_stop = dict_fill_bmodes[filln]['t_stop_SQUEEZE']

fillinfo = mdb.getLHCFillData(filln)
bmodes = fillinfo['beamModes']
for bm in bmodes:
    if bm['mode'] == 'INJPHYS':
        t_start = bm['endTime']
    if bm['mode'] == 'RAMP':
        t_stop = bm['endTime']+10*60.
        break
        

data = {}

data.update(ldb.get([
			'LHC.BCTFR.A6R4.B1:BUNCH_INTENSITY',
			'LHC.BCTFR.A6R4.B2:BUNCH_INTENSITY'], t_start, t_stop))
print 'Downloaded Intensity'

data.update(ldb.get([
			'HX:BETASTAR_IP1',
            'CMS:LUMI_TOT_INST',
            'ATLAS:LUMI_TOT_INST',
            'LHC.BSRA.US45.B1:ABORT_GAP_ENERGY',
				], t_start, t_stop))
print 'Downloaded Beta'

import matplotlib.pyplot as plt
import matplotlib.cm as cm
plt.close('all')
ms.mystyle_arial(fontsz=14, dist_tick_lab=5)

axt = None
axbun = None
axbet = None

figlist = []

for beam in [1,2]:
    spbet_list = []

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
    lifet_h[lifet_h<0]= 200

    t_ref = t_stamps[0]
    t_fbct_minutes = (t_stamps-t_ref)/60.
    

    fig = plt.figure(beam, figsize=(8*1.5, 6))
    fig.set_facecolor('w')
    sploss = plt.subplot2grid(shape=(1, 4), loc=(0, 1),colspan=3, sharey=axt, sharex = axbun)
    axt = sploss
    axbun = sploss
    spbet = plt.subplot2grid(shape=(1, 4), loc=(0, 0),colspan=1, sharey=axt, sharex=axbet)
    spbet_list.append(spbet)

    cc = sploss.pcolormesh(slots, t_fbct_minutes, bint_norm)
    fig.colorbar(cc, ax=sploss, label='Accumulated losses [%]')
    sploss.set_xlabel('25ns slot')


    t_temp = data['HX:BETASTAR_IP1'][0]
    t_betast_minutes = (t_temp-t_ref)/60.

    figl = plt.figure(100+beam, figsize=(8*1.5, 6))
    figl.set_facecolor('w')
    splifet = plt.subplot2grid(shape=(1, 4), loc=(0, 1),colspan=3, sharey=axt, sharex = axbun)
    spbet = plt.subplot2grid(shape=(1, 4), loc=(0, 0),colspan=1, sharey=axt, sharex=axbet)
    spbet_list.append(spbet)

    cc = splifet.pcolormesh(slots, t_fbct_minutes[1:], lifet_h, vmin=0, vmax = 120, cmap=cm.jet_r)
    splifet.set_xlim(0, 3500)
    splifet.set_xlabel('25ns slot')
    figl.colorbar(cc, ax=splifet, label='Lifetime [h]')

    
    for spbet in spbet_list:
        spbet.plot(data['HX:BETASTAR_IP1'][1], t_betast_minutes, lw=2)
        #spbet.set_xlim(25., 110.)
        spbet.set_xlabel('Beta* [cm]')
        spbet.set_ylabel('Time [min]')
        spbet.xaxis.set_major_locator(MaxNLocator(5))
        spbet.xaxis.label.set_color('b')
        spbet.tick_params(axis='x', colors='b')

        splumi = spbet.twiny()
        var = 'LHC.BSRA.US45.B1:ABORT_GAP_ENERGY'; splumi.plot(data[var][1], (data[var][0]-t_ref)/60., color='r', lw=2)
        splumi.xaxis.set_major_locator(MaxNLocator(5))
        spbet.grid('on')
        splumi.set_xlabel('Energy [GeV]')
        splumi.set_ylim(bottom=(t_start-t_ref)/60.,top=(t_stop-t_ref)/60.)
        splumi.xaxis.label.set_color('r')
        splumi.tick_params(axis='x', colors='red')

        fig.suptitle('Fill %d Losses Beam %d'%(filln, beam))
        fig.subplots_adjust(left=.05, right=1.)

        figl.suptitle('Fill %d Lifetime Beam %d'%(filln, beam))
        figl.subplots_adjust(left=.05, right=1.)

        figlist += [fig, figl]

# plt.figure(3)
# spl = plt.subplot(2,1,1)
# plt.plot(bint[-1,:])
# plt.subplot(2,1,2, sharex=spl)
# plt.plot(lifet_h[-1,:])

# To save
# for fg in figlist: fg.savefig(ff+'/'+fg._suptitle.get_text().replace(' ','_')+'_zoom.png', dpi=200)


plt.show()
