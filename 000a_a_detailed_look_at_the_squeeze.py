import sys, os
BIN = os.path.expanduser("../LHC_fullRun2_analysis_scripts/")
sys.path.append(BIN)

import LHCMeasurementTools.mystyle as ms

import pickle
import numpy as np
import scipy.io as sio

from matplotlib.ticker import MaxNLocator

filln = 6639
#~ filln = 6646
#~ # filln = 6060
#~ filln = 6645
#~ filln=6640
#~ filln=6641
#~ filln=6650
#~ filln=6662
#~ filln = 6666
#~ filln = 6654
#~ #filln = 6674
#~ filln=6677
#~ filln=6681
#~ filln=6714

ff = './results_squeeze/'

if len(sys.argv)>1:
   filln = int(sys.argv[1])

beta_obs_cm = 33.
#beta_obs_cm = 43.



import pytimber
ldb = pytimber.LoggingDB(source='ldb')
mdb = pytimber.LoggingDB(source='mdb')

fillinfo = mdb.getLHCFillData(filln)
bmodes = fillinfo['beamModes']
for bm in bmodes:
    if bm['mode'] == 'FLATTOP':
        t_start = bm['startTime']
    if bm['mode'] == 'ADJUST':
        t_stop = bm['endTime']+5*60.
        break

#hack
t_stop = t_start+20*60.
        

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

data.update(ldb.getScaled(['RPMBB.RR17.ROD.A12B1:I_MEAS', 'RPMBB.RR17.ROD.A12B2:I_MEAS'],
                    t_start, t_stop, scaleAlgorithm='AVG', scaleInterval='MINUTE',scaleSize='1'))
print 'Downloaded Octupole current'

import matplotlib.pyplot as plt
import matplotlib.cm as cm
plt.close('all')
ms.mystyle_arial(fontsz=16, dist_tick_lab=5)

axt = None
axbun = None
axbet = None
ax_lifet_obs = None

figlist = []

fig_obs = plt.figure(200, figsize=(8*1.5, 6))
fig_obs.set_facecolor('w')

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

    cc = splifet.pcolormesh(slots, t_fbct_minutes[1:], lifet_h, vmin=0, vmax=120, cmap=cm.jet_r)
    splifet.set_xlim(0, 3500)
    splifet.set_xlabel('25ns slot')
    figl.colorbar(cc, ax=splifet, label='Lifetime [h]')

    
    for spbet in spbet_list:
        spbet.plot(data['HX:BETASTAR_IP1'][1], t_betast_minutes, lw=2)
        spbet.set_xlim(25., 110.)
        spbet.set_xlabel('Beta* [cm]')
        spbet.set_ylabel('Time [min]')
        spbet.xaxis.set_major_locator(MaxNLocator(5))
        spbet.xaxis.label.set_color('b')
        spbet.tick_params(axis='x', colors='b')

        splumi = spbet.twiny()
        var = 'ATLAS:LUMI_TOT_INST'; splumi.plot(data[var][1]/1e4*100, (data[var][0]-t_ref)/60., color='r', lw=2)
        
        var = 'RPMBB.RR17.ROD.A12B%d:I_MEAS'%beam; splumi.plot(-data[var][1], (data[var][0]-t_ref)/60., color='k', lw=2)
        
        splumi.xaxis.set_major_locator(MaxNLocator(5))
        spbet.grid('on')
        splumi.set_xlabel('Lumi [1e32 cm^-2.s^-1], Oct. curr [A]')
        splumi.set_ylim(bottom=0.,top=(t_stop-t_ref)/60.)
        #~ splumi.xaxis.label.set_color('r')
        #~ splumi.tick_params(axis='x', colors='red')

        fig.suptitle('Fill %d Losses Beam %d'%(filln, beam))
        fig.subplots_adjust(left=.07, right=1., top=.85)

        figl.suptitle('Fill %d Lifetime Beam %d'%(filln, beam))
        figl.subplots_adjust(left=.07, right=1., top=.85)

        figlist += [figl]
        
    sortbeta = np.argsort(data['HX:BETASTAR_IP1'][1])
    t_obs = np.interp(beta_obs_cm, 
                        np.take(data['HX:BETASTAR_IP1'][1], sortbeta),
                        np.take(data['HX:BETASTAR_IP1'][0], sortbeta))
                        
    i_obs = np.argmin(np.abs(t_stamps-t_obs))
    t_fbct_obs = t_stamps[i_obs]
    
    
    lifet_obs = np.mean(lifet_h[i_obs-1:i_obs+2, :], axis=0)
    avg_lifet = 1./np.mean(1./lifet_obs[lifet_obs!=0.])
    
    
    lr_obs = 1./lifet_obs*100.
    lr_obs[lifet_obs==0.] = -1
    
    var = 'RPMBB.RR17.ROD.A12B%d:I_MEAS'%beam
    oct_obs = -np.interp(t_obs, data[var][0], data[var][1])
    
    if beam==1:
        ax_lifet_obs = fig_obs.add_subplot(1,1,1, sharex=axbun, sharey=ax_lifet_obs)
    ax_lifet_obs.plot(lr_obs, color={1:'b', 2:'r'}[beam], lw=2, label='B%d, tau=%.1fh'%(beam, avg_lifet))
    ax_lifet_obs.legend(loc='upper left', prop={'size':16})
                        
    ax_lifet_obs.set_xlim(0, 3500)
    ax_lifet_obs.set_ylim(bottom=0., top=30)
    ax_lifet_obs.grid('on')
    
    with open('avglifet/avg_lifet_fill%d_beam%d_betast%.1fcm.txt'%(filln, beam, beta_obs_cm), 'w') as fid:
        fid.write('%.3f %.3f'%(oct_obs, avg_lifet))

    fnamemat = ff+'/lossrate_fill%d_beam%d_betast%.1fcm'%(filln, beam, beta_obs_cm)
    sio.savemat(fnamemat, {'lossrate':lr_obs}, oned_as='row')

    
fig_obs.suptitle('Fill %d Loss rate at betast %.1fcm\nIoct=%.1fA'%(filln, beta_obs_cm, oct_obs))
fig_obs.subplots_adjust(left=.07, right=.95)
ax_lifet_obs.set_xlabel('25 ns slot')
ax_lifet_obs.set_ylabel('Loss rate [%/h]')

figlist.append(fig_obs)



[fg.savefig(ff+'/'+fg._suptitle.get_text().split('\n')[0].replace(' ','_')+'.png', dpi=200) for fg in figlist]
ax_lifet_obs.set_xlim(770, 1270);
[fg.savefig(ff+'/'+fg._suptitle.get_text().split('\n')[0].replace(' ','_')+'_zoom.png', dpi=200) for fg in figlist]


#~ plt.show()
