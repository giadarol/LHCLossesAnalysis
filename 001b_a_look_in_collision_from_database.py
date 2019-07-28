import sys, os
import numpy as np
import matplotlib.pyplot as plt

BIN = os.path.expanduser("../LHC_fullRun2_analysis_scripts/")
sys.path.append(BIN)

filln = 7174
T_observ_h = 5
t_obs_minutes = 30
beam = 1

filln = 7145
T_observ_h = 5
t_obs_minutes = 50 
beam = 1

import LHCMeasurementTools.mystyle as ms
import LHCMeasurementTools.TimestampHelpers as th

ms.mystyle_arial(fontsz=16, dist_tick_lab=5)
import pytimber
ldb = pytimber.LoggingDB(source='ldb')

fillinfo = ldb.getLHCFillData(filln)
bmodes = fillinfo['beamModes']
for bm in bmodes:
    if bm['mode'] == 'FLATTOP':
        t_start = bm['startTime']

t_stop = t_start + T_observ_h*3600

data = {}

data.update(ldb.get([
            'LHC.BCTFR.A6R4.B1:BUNCH_INTENSITY',
            'LHC.BCTFR.A6R4.B2:BUNCH_INTENSITY'], t_start, t_stop))
print 'Downloaded Intensity'
data.update(ldb.get([
            'CMS:LUMI_TOT_INST',
            'ATLAS:LUMI_TOT_INST'], t_start, t_stop))
print 'Downloaded luminosity'



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

axt = None
axbun = None
axbet = None

fig = plt.figure(beam, figsize=(8*1.5, 6))
fig.set_facecolor('w')
sploss = plt.subplot2grid(shape=(1, 4), loc=(0, 1),colspan=3,
        sharey=axt, sharex = axbun)
axt = sploss
axbun = sploss
spbet = plt.subplot2grid(shape=(1, 4), loc=(0, 0),colspan=1,
        sharey=axt, sharex=axbet)

cc = sploss.pcolormesh(slots, t_fbct_minutes, bint_norm, 
        vmin=0, vmax=10)
fig.colorbar(cc, ax=sploss, label='Accumulated losses [%]')
sploss.set_xlabel('25ns slot')
for ee in ['ATLAS', 'CMS']:
    dd = data['%s:LUMI_TOT_INST'%ee]
    tt = (dd[0]-t_start)/60.
    spbet.plot(dd[1]*1e-3, tt)

t_obs = t_start + t_obs_minutes*60

i_obs = np.argmin(np.abs(t_stamps-t_obs))
lifet_obs = np.mean(lifet_h[i_obs-1:i_obs+2, :], axis=0)
avg_lifet = 1./np.mean(1./lifet_obs[lifet_obs!=0.])
    
t_fbct_obs = t_stamps[i_obs]
lr_obs = 1./lifet_obs*100.
lr_obs[lifet_obs==0.] = -1

fig_obs = plt.figure(200, figsize=(8*1.5, 6))
fig_obs.set_facecolor('w')
ax_lifet_obs = fig_obs.add_subplot(1,1,1, sharex=axbun)
ax_lifet_obs.plot(lr_obs, color={1:'b', 2:'r'}[beam], lw=2,
         label='B%d, tau=%.1fh'%(beam, avg_lifet))
ax_lifet_obs.legend(loc='upper left', prop={'size':16})
                    
ax_lifet_obs.set_xlim(0, 3500)
ax_lifet_obs.set_ylim(bottom=0., top=30)
ax_lifet_obs.grid('on')

ax_lifet_obs.set_xlabel('25 ns slot')
ax_lifet_obs.set_ylabel('Loss rate [%/h]')

fig300 = plt.figure(300)
plt.plot(bint_norm[i_obs, :])
plt.show()
