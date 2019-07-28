import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

plt.close('all')

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

sigma_m2 = 80e-3*1e-28

import LHCMeasurementTools.mystyle as ms
import LHCMeasurementTools.TimestampHelpers as th

ms.mystyle_arial(fontsz=16, dist_tick_lab=5)

#################
# Download data #
#################
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
            'ATLAS:LUMI_TOT_INST',
            'CMS:BUNCH_LUMI_INST',
            'ATLAS:BUNCH_LUMI_INST'], t_start, t_stop))
print 'Downloaded luminosity'

bint = data['LHC.BCTFR.A6R4.B%d:BUNCH_INTENSITY'%beam][1]
t_stamps = data['LHC.BCTFR.A6R4.B%d:BUNCH_INTENSITY'%beam][0]

###################################
# Compute lumi on FBCT time stamp #
###################################
lumi_m2s = np.zeros_like(bint)
for ii in range(len(bint[0, :])):
    lumi_m2s[:, ii] += 1e34 * np.interp(t_stamps,
         data['CMS:BUNCH_LUMI_INST'][0],
         np.array(data['CMS:BUNCH_LUMI_INST'][1])[:, ii])
    lumi_m2s[:, ii] += 1e34 * 1e-3 * np.interp(t_stamps,
         data['ATLAS:BUNCH_LUMI_INST'][0],
         np.array(data['ATLAS:BUNCH_LUMI_INST'][1])[:, ii])



bint[bint<0.5e11]=1.

slots = np.arange(len(bint[0, :]))

bint_norm = bint*0

for ii in xrange(len(t_stamps)):
    bint_norm[ii, :] = (1.-bint[ii, :]/bint[0, :])*100.

t_mat = np.dot(np.ones((len(slots), 1)), np.atleast_2d(t_stamps)).T
loss_rate = -np.diff(bint, axis=0)/np.diff(t_mat, axis=0)
lifet_h = 1/(loss_rate/bint[:-1,:])/3600.
lifet_h[lifet_h<0]= 200

BO_loss_rate = 0.5*(lumi_m2s[:-1, :]+lumi_m2s[1:, :])*sigma_m2
loss_rate_woBO = loss_rate-BO_loss_rate
lifet_woBO_h = 1/(loss_rate_woBO/bint[:-1,:])/3600.
lifet_woBO_h[lifet_woBO_h<0]= 200
lifet_woBO_h[lifet_woBO_h>200]= 200


t_ref = t_stamps[0]
t_fbct_minutes = (t_stamps-t_ref)/60.

time_conv = th.TimeConverter(time_in='h', t_plot_tick_h=None, t_ref=t_stamps[0])
tc = time_conv.from_unix


fig = plt.figure(beam, figsize=(8*1.8,6*1.3))
fig.set_facecolor('w')
axlt = plt.subplot2grid(shape=(5, 5), loc=(0, 1), colspan=3, rowspan=4)

cc=axlt.pcolormesh(np.arange(3564), tc(t_stamps), lifet_woBO_h, 
    cmap=cm.jet_r, vmin=0, vmax=120)
axcb = plt.subplot2grid(shape=(5, 5), loc=(4, 1), colspan=3, rowspan=1)
plt.colorbar(cc, cax=axcb, label='Lifetime (BO corrected) [h]', orientation='horizontal')

axlt.set_xlabel('25ns slot')
axslot = axlt
axslot.set_xlim(0, 3500)

ax1 = plt.subplot2grid(shape=(5, 5), loc=(0, 0), colspan=1, rowspan=4, sharey=axlt)
ax1.step(fill_data['xing_angle'][1]*1e6/2, tc(t_stamps), lw=2.)
ax1.set_xlim(120, 170)
ax1.set_xlabel('Half crossing angle [urad]')
ax1.set_ylabel('Time [h]')
ax1.xaxis.label.set_color('b')
ax1.tick_params(axis='x', colors='b')
ax1.grid('on')

axlumi = ax1.twiny()
axlumi.plot(0.5*np.sum(lumi_data['CMS']['bunch_lumi']+lumi_data['ATLAS']['bunch_lumi'], 
            axis=1)/1e34*1e-4, tc(t_stamps), 'r', lw=2.)
axlumi.set_xlabel('Avg lumi [1e34 cm^-2.s-1]')
axlumi.xaxis.label.set_color('r')
axlumi.tick_params(axis='x', colors='r')
axlumi.xaxis.set_major_locator(MaxNLocator(5))

axqp = plt.subplot2grid(shape=(5, 5), loc=(0, 4), colspan=1, rowspan=4, sharey=axlt)
if enable_QP:
    for plane, styl in zip(['H', 'V'], ['--', '-']):
        parname = 'LHCBEAM%d/QP'%beam + plane
        thistrim = chromaTrims_end[parname]
        thistrim.data.append(thistrim.data[-1])
        thistrim.time.append(t_stop_SB)
        axqp.step(np.array(thistrim.data)+QP_offsets[parname], tc(thistrim.time), ls=styl, lw=2, color = 'k')
axqp.set_xlim(0, 16)
axqp.set_xlabel('Chromaticity')
axqp.grid('on')

axoct = axqp.twiny()
thisoct = data_timb['RPMBB.RR17.ROD.A12B%d:I_MEAS'%beam]
axoct.plot(np.abs(thisoct[1]), tc(thisoct[0]), 'r', lw=2)
axoct.set_xlim(200, 550)
axoct.xaxis.label.set_color('r')
axoct.tick_params(axis='x', colors='r')
axoct.xaxis.set_major_locator(MaxNLocator(5))
axoct.set_xlabel('I octupoles [A]')
axoct.set_ylim(bottom=0)

fig.subplots_adjust(right=.95, left=.05, bottom=.12, top=.81, hspace = 1)

fig.suptitle('Fill %d SB BurnOff corrected lifetime B%d\nSB started on %s\n'%(filln, beam, tref_string))



prrrrrrr


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

############################
# Look at specific instant #
############################

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
