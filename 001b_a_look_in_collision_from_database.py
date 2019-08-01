import sys, os
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator

from FillingPatterns import filling_pattern as fp

mode = 'loss_rate'

assert(mode in ['loss_rate', 'lifetime'])

plt.close('all')

BIN = os.path.expanduser("../LHC_fullRun2_analysis_scripts/")
sys.path.append(BIN)

group_definitions = [
    {'color':'mediumseagreen', 'slots_within_injection': [1,2,3,4,5], 'n_bunches_in_injection':144},
    {'color':'green', 'slots_within_injection': [18, 19, 20, 21, 22], 'n_bunches_in_injection':144},
    {'color':'indianred', 'slots_within_injection': [135, 136, 137, 138, 139], 'n_bunches_in_injection':144},
    {'color':'darkred', 'slots_within_injection': [153, 154, 155, 156, 157], 'n_bunches_in_injection':144},
    ]

# MD large telescope
filln = 7174
T_observ_h = 5
t_detail_h = .5
beam = 1

# # Physics fill
# filln = 7145
# T_observ_h = 5
# t_detail_h = .8
# beam = 1

# # Physics fill
# filln = 7236
# T_observ_h = 10
# t_detail_h = 3
# beam = 1

# Physics fill (constant angle)
filln = 7266
T_observ_h = 10
t_detail_h = 3
beam = 1

# # Only beam 1
# filln = 6966
# T_observ_h = 5
# t_detail_h = 1.7
# beam = 1



dt_minutes = 5 


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

tref_string=time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime(t_start))

data = {}

data.update(ldb.get([
            'LHC.BCTFR.A6R4.B1:BUNCH_INTENSITY',
            'LHC.BCTFR.A6R4.B2:BUNCH_INTENSITY'], t_start, t_stop))
print 'Downloaded Intensity'
data.update(ldb.get([
            'CMS:LUMI_TOT_INST',
            'ATLAS:LUMI_TOT_INST',
            'CMS:BUNCH_LUMI_INST',
            'ATLAS:BUNCH_LUMI_INST',
            'LHC.RUNCONFIG:IP1-XING-V-MURAD',
            'LHC.RUNCONFIG:IP5-XING-H-MURAD'], t_start, t_stop))
print 'Downloaded luminosity'

# Get angle first and  point
varnames = ['LHC.RUNCONFIG:IP1-XING-V-MURAD', 'LHC.RUNCONFIG:IP5-XING-H-MURAD']
for var in varnames:
    dpoint_start = ldb.get(var, t_start)
    dpoint_stop = ldb.get(var, t_stop)
    data[var] = list(data[var])
    data[var][0] = np.concatenate(([t_start], data[var][1], [t_stop]))
    data[var][1] = np.concatenate((dpoint_start[var][1], data[var][1], dpoint_stop[var][1]))

bint_raw = data['LHC.BCTFR.A6R4.B%d:BUNCH_INTENSITY'%beam][1]
t_stamps_raw = data['LHC.BCTFR.A6R4.B%d:BUNCH_INTENSITY'%beam][0]

t_stamps = np.arange(t_start, t_stop, dt_minutes*60)
bint = np.zeros((len(t_stamps), len(bint_raw[0, :])))
for ii in range(len(bint[0, :])):
    bint[:, ii] = np.interp(t_stamps, t_stamps_raw, bint_raw[:, ii])

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

mask_filled = bint[0, :] > 5e10
patt = fp.Filling_Pattern_Single_Beam(np.int_(mask_filled))

for gg in group_definitions:
    gg['mask'] = patt.belongs_to_group(slots_within_injection=gg['slots_within_injection'],
                                       n_bunches_in_injection=gg['n_bunches_in_injection'])

figpatt = plt.figure(1000)
axpatt = figpatt.add_subplot(1,1,1)
axpatt.plot(patt.pattern)
for gg in group_definitions:
    axpatt.plot(gg['mask'], '.')
    
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
axgr = plt.subplot2grid(shape=(5, 5), loc=(0, 4), colspan=1, rowspan=4, sharey=axlt)

if mode == 'lifetime':
    cc=axlt.pcolormesh(np.arange(3564), tc(t_stamps), lifet_woBO_h, 
        cmap=cm.jet_r, vmin=0, vmax=120)
    axcb = plt.subplot2grid(shape=(5, 5), loc=(4, 1), colspan=3, rowspan=1)
    plt.colorbar(cc, cax=axcb, label='Lifetime (BO corrected) [h]', orientation='horizontal')
elif mode == 'loss_rate':
    cc=axlt.pcolormesh(np.arange(3564), tc(t_stamps[:-1]), loss_rate_woBO*3600/1e9, 
        cmap=cm.jet, vmin=0, vmax=4)
    axcb = plt.subplot2grid(shape=(5, 5), loc=(4, 1), colspan=3, rowspan=1)
    plt.colorbar(cc, cax=axcb, label='Loss rate (BO corrected) [10^9 p/h]', orientation='horizontal')

for gg in group_definitions:
    axgr.plot(np.mean(loss_rate_woBO[:, gg['mask']], axis=1)*1e-9*3600, tc(t_stamps[:-1]),
            color=gg['color'], linewidth=2)
axgr.grid(True)
axgr.set_xlabel('Loss rate (BO corrected)\n[10^9 p/h]')
axgr.set_xlim(-0.5, 5)

axlt.set_xlabel('25ns slot')
axslot = axlt
axslot.set_xlim(0, 3500)

ax1 = plt.subplot2grid(shape=(5, 5), loc=(0, 0), colspan=1, rowspan=4, sharey=axlt)
ax1.step(data['LHC.RUNCONFIG:IP1-XING-V-MURAD'][1], 
        tc(data['LHC.RUNCONFIG:IP1-XING-V-MURAD'][0]), lw=2.)
ax1.set_xlim(120, 170)
ax1.set_xlabel('Half crossing angle [urad]')
ax1.set_ylabel('Time [h]')
ax1.xaxis.label.set_color('b')
ax1.tick_params(axis='x', colors='b')
ax1.grid('on')

axlumi = ax1.twiny()
axlumi.plot(0.5*np.sum(lumi_m2s, axis=1)/1e34*1e-4, tc(t_stamps), 'r', lw=2.)
axlumi.set_xlabel('Avg lumi [1e34 cm^-2.s-1]')
axlumi.xaxis.label.set_color('r')
axlumi.tick_params(axis='x', colors='r')
axlumi.xaxis.set_major_locator(MaxNLocator(5))


# axqp = plt.subplot2grid(shape=(5, 5), loc=(0, 4), colspan=1, rowspan=4, sharey=axlt)
# if enable_QP:
#     for plane, styl in zip(['H', 'V'], ['--', '-']):
#         parname = 'LHCBEAM%d/QP'%beam + plane
#         thistrim = chromaTrims_end[parname]
#         thistrim.data.append(thistrim.data[-1])
#         thistrim.time.append(t_stop_SB)
#         axqp.step(np.array(thistrim.data)+QP_offsets[parname], tc(thistrim.time), ls=styl, lw=2, color = 'k')
# axqp.set_xlim(0, 16)
# axqp.set_xlabel('Chromaticity')
# axqp.grid('on')

# axoct = axqp.twiny()
# thisoct = data_timb['RPMBB.RR17.ROD.A12B%d:I_MEAS'%beam]
# axoct.plot(np.abs(thisoct[1]), tc(thisoct[0]), 'r', lw=2)
# axoct.set_xlim(200, 550)
# axoct.xaxis.label.set_color('r')
# axoct.tick_params(axis='x', colors='r')
# axoct.xaxis.set_major_locator(MaxNLocator(5))
# axoct.set_xlabel('I octupoles [A]')
# axoct.set_ylim(bottom=0)

fig.subplots_adjust(right=.95, left=.05, bottom=.12, top=.81, hspace = 1)

fig.suptitle('Fill %d SB BurnOff corrected lifetime B%d\nFT started on %s\n'%(filln, beam, tref_string))

i_detail = np.argmin(np.abs(t_stamps - (t_start+t_detail_h*3600)))
t_h_detail = tc(t_stamps)[i_detail]
axd = None
axd2 = None
if t_h_detail<tc(t_stamps)[-1]:
    figd = plt.figure(100+beam, figsize=(8*1.4,6))
    figd.set_facecolor('w')
    axd = figd.add_subplot(111, sharex=axslot, sharey=axd)
    i_plot = np.argmin(np.abs(t_h_detail-tc(t_stamps)))
    BO_lr_det = BO_loss_rate[i_plot, :]
    lr_det = loss_rate[i_plot, :]
    bint_det = bint[i_plot, :]

    mask_filled_det = bint_det>2e10
    BO_lr_det[~mask_filled_det] = 0.
    lr_det[~mask_filled_det] = 0.

    lt_BO_avg = 1./(np.sum(BO_lr_det)/np.sum(bint[i_plot]))/3600.
    lt_other = 1./(np.sum(lr_det-BO_lr_det)/np.sum(bint[i_plot]))/3600.

    plt.fill_between(np.arange(3564),BO_lr_det, color='green', alpha=0.6, label='Burn off')
    plt.fill_between(np.arange(3564), BO_lr_det, lr_det, color='red', alpha=.6, label='Additional losses')
    ms.sciy()
    axd.set_ylabel('Loss rate [p/s]')
    axd.set_xlabel('25 ns slot')
    axd.legend(loc='lower right', prop={'size':14})
    axd.set_xlim(0,3500)
    axd.set_ylim(bottom=0)
    axd.grid('on')
    figd.subplots_adjust(right=.94, left=.1, bottom=.12, top=.86)
    figd.suptitle('Fill %d SB Loss Rate at %.1fh for B%d\nSB started on %s\n'%(filln, t_h_detail, beam, tref_string))
    

    figd2 = plt.figure(200+beam, figsize=(8*1.4,6))
    figd2.set_facecolor('w')
    axd2 = figd2.add_subplot(111, sharex=axslot, sharey=axd2)

    axd2.plot(np.arange(3564), BO_lr_det/bint_det*3600*100, 'b', lw=2, alpha=.4, label='Burn-off (tau=%.1fh)'%lt_BO_avg)
    axd2.plot(np.arange(3564), (lr_det-BO_lr_det)/bint_det*3600*100, 'r', lw=2, label='Other losses (tau=%.1fh)'%lt_other)
    ms.sciy()
    axd2.set_ylabel('Loss rate [%/h]')
    axd2.set_xlabel('25 ns slot')
    axd2.legend(loc='upper right', prop={'size':14})
    axd2.set_xlim(0,3500)
    axd2.set_ylim(bottom=0, top=5.)
    axd2.grid('on')
    figd2.subplots_adjust(right=.94, left=.1, bottom=.12, top=.86)
    figd2.suptitle('Fill %d Loss Rates at %.1fh for B%d\nFT started on %s\n'%(filln, t_h_detail, beam, tref_string))


plt.show()


