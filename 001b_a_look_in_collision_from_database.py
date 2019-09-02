import sys, os
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator

from FillingPatterns import filling_pattern as fp

mode = 'loss_rate'

assert(mode in ['loss_rate', 'lifetime', 'integrated'])

plt.close('all')

BIN = os.path.expanduser("../LHC_fullRun2_analysis_scripts/")
sys.path.append(BIN)

# For 25 ns 3x48b
group_definitions = [
    {'color':'mediumseagreen', 'slots_within_injection': [1,2,3,4,5], 'n_bunches_in_injection':144},
    {'color':'green', 'slots_within_injection': [18, 19, 20, 21, 22], 'n_bunches_in_injection':144},
    {'color':'indianred', 'slots_within_injection': [135, 136, 137, 138, 139], 'n_bunches_in_injection':144},
    {'color':'darkred', 'slots_within_injection': [153, 154, 155, 156, 157], 'n_bunches_in_injection':144},
    ]

# #For 8b+4e 2x96b
# group_definitions = [
#     {'color':'mediumseagreen', 'slots_within_injection': [1,2,3,4,5], 'n_bunches_in_injection':96},
#     {'color':'green', 'slots_within_injection': [36, 37, 38, 39, 40], 'n_bunches_in_injection':96},
#     {'color':'indianred', 'slots_within_injection': [112, 113, 114, 115, 116], 'n_bunches_in_injection':96},
#     {'color':'darkred', 'slots_within_injection': [136, 137, 138, 139, 140], 'n_bunches_in_injection':96},
#     ]

# #For 8b+4e 4x32b
# group_definitions = [
#     {'color':'mediumseagreen', 'slots_within_injection': [1,2,3,4,5], 'n_bunches_in_injection':128},
#     {'color':'green', 'slots_within_injection': [36, 37, 38, 39, 40], 'n_bunches_in_injection':128},
#     {'color':'indianred', 'slots_within_injection': [114, 115, 116, 117, 118], 'n_bunches_in_injection':128},
#     {'color':'darkred', 'slots_within_injection': [138, 139, 140, 141, 142], 'n_bunches_in_injection':128},
#     ]


# outp_folder = '/eos/user/g/giadarol/temp/20190801_losses_wp2/'
outp_folder = None

# # Test injection 
# filln = 6610
# T_download_h = 6
# slotrange = np.array([0, 520]) + 520
# T_h_range_colorplot = [0, T_download_h]
# beam = 1
# t_detail_h_list = [0.] # at the beginning of the fill
# bmode_start = 'INJPHYS'
# delay_start_h = .4

# # MD large telescope
# filln = 7174
# T_download_h = 5
# slotrange = np.array([0, 520]) + 520
# T_h_range_colorplot = [.2, 5.]
# beam = 1
# t_detail_h_list = [.3, 2.] # at the beginning of the fill
# bmode_start = 'FLATTOP'
# delay_start_h = 0.

# # Physics fill (2017)
# filln = 6060
# T_download_h = 15
# slotrange = np.array([0, 520])+1270
# T_h_range_colorplot = [.3, 14.1]
# t_detail_h_list = [0.83, 1.0, 2., 3.5, 3.58, 5.33, 7.17, 8.33, 10, 13]
# beam = 1 
# bmode_start = 'FLATTOP'
# delay_start_h = 0.

# # Physics fill (2017, 8b+4e)
# filln = 6315# 6305
# T_download_h = 14
# slotrange = np.array([0, 520])+1270
# T_h_range_colorplot = [.3, 14.1]
# t_detail_h_list = [2,5,8]
# beam = 1 
# bmode_start = 'FLATTOP'
# delay_start_h = 0.

# # Physics fill (2017, 50 ns)
# filln = 5980 
# T_download_h = 7
# slotrange = np.array([0, 520])+1270
# T_h_range_colorplot = [.3, 14.1]
# t_detail_h_list = [2,5,8]
# beam = 1 
# bmode_start = 'FLATTOP'
# delay_start_h = 0.

# # High-intensity 8b+4e 
# filln = 7366
# T_download_h = 5
# slotrange = np.array([0, 520])+1270
# T_h_range_colorplot = [.3, 2.5]
# t_detail_h_list = [.9]
# beam = 1 
# bmode_start = 'FLATTOP'
# delay_start_h = 0.

# Physics fill
filln = 7236
T_download_h = 10
slotrange = np.array([0, 520])+1270
T_h_range_colorplot = [.3, 9.5]
t_detail_h_list = [0.58, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9.2]
t_detail_h_list = [2]
beam = 1 
bmode_start = 'FLATTOP'
delay_start_h = 0.

# # Physics fill (constant angle)
# filln = 7266
# T_download_h = 10
# slotrange = np.array([0, 520])+1270
# T_h_range_colorplot = [.3, 9.2]
# beam = 1
# t_detail_h_list = [0.5, 2, 3, 4, 5, 6, 7]
# bmode_start = 'FLATTOP'
# delay_start_h = 0.

# # Physics fill (very long)
# filln = 7056
# T_download_h = 30
# slotrange = np.array([0, 520])+1270
# T_h_range_colorplot = [.3, 25.9]
# t_detail_h_list = [0.33, 1, 2,3,4,5,6,7,7.5,8,9,9.5,10,11,12,12.5,13,14,15,16,17,18,19,20,21,22,23,24,25]
# beam = 1 
# bmode_start = 'FLATTOP'
# delay_start_h = 0.

# # Only beam 1
# filln = 6966
# T_download_h = 5
# slotrange = np.array([0, 520])+1270
# T_h_range_colorplot = [.3, 2.77]
# t_detail_h_list = [1.7, 2.55, 2.7]
# # t_detail_h = 1.7 # 30 cm, 130 urad, nominal tunes
# # t_detail_h = 2.55 # 25 cm, 130 urad, nominal tunes
# # t_detail_h = 2.7 # 25 cm, 130 urad, modified tunes
# beam = 1
# bmode_start = 'FLATTOP'
# delay_start_h = 0.
  
# # Only beam 2 
# filln = 6967
# T_download_h = 5
# t_detail_h = 1.2 # 30 cm, 130 urad, nominal tunes
# t_detail_h = 1.3 # 25 cm ??
# slotrange = np.array([0, 520])+1270
# beam = 2
# bmode_start = 'FLATTOP'
# delay_start_h = 0.

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
    if bm['mode'] == bmode_start:
        t_start = bm['startTime'] + delay_start_h * 3600.

t_stop = t_start + T_download_h*3600

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
    data[var][0] = np.concatenate(([t_start], data[var][0], [t_stop]))
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

bint[bint<0.2e11]=1.

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

acc_BO_loss_rate = np.cumsum(BO_loss_rate*np.diff(t_mat, axis=0), axis=0)
acc_loss_rate_woBO = np.cumsum(loss_rate_woBO*np.diff(t_mat, axis=0), axis=0)

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
elif mode == 'integrated':
    cc=axlt.pcolormesh(np.arange(3564), tc(t_stamps[:-1]), acc_loss_rate_woBO/1e10, 
        cmap=cm.jet, vmin=0, vmax=4)
    axcb = plt.subplot2grid(shape=(5, 5), loc=(4, 1), colspan=3, rowspan=1)
    plt.colorbar(cc, cax=axcb, label='Accumulated losses (BO corrected) [10^10 p]', orientation='horizontal')

for gg in group_definitions:
    if mode == 'loss_rate' or mode == 'lifetime':
        axgr.plot(np.mean(loss_rate_woBO[:, gg['mask']], axis=1)*1e-9*3600, tc(t_stamps[:-1]),
            color=gg['color'], linewidth=2)
    elif mode == 'integrated':
        axgr.plot(np.mean(acc_loss_rate_woBO[:, gg['mask']], axis=1)*1e-10, tc(t_stamps[:-1]),
            color=gg['color'], linewidth=2)
 

axgr.grid(True)
if mode == 'loss_rate' or mode == 'lifetime':
    axgr.set_xlabel('Loss rate (BO corrected)\n[10^9 p/h]')
    axgr.set_xlim(-0.5, 5)
elif mode == 'integrated':
    axgr.set_xlabel('Accumulated losses (BO corrected)\n[10^10 p]')
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

fig.suptitle('Fill %d SB BurnOff corrected lifetime B%d\n%s started on %s\n'%(filln, beam, bmode_start, tref_string))

for t_detail_h in t_detail_h_list:
    i_detail = np.argmin(np.abs(t_stamps - (t_start+t_detail_h*3600)))
    t_det_h_found = tc(t_stamps)[i_detail]
    axd = None
    axd2 = None
    if t_det_h_found<tc(t_stamps)[-1]:
    #    figd = plt.figure(100+beam, figsize=(8*1.4,6))
    #    figd.set_facecolor('w')
    #    axd = figd.add_subplot(111, sharex=axslot, sharey=axd)
    
        i_plot = np.argmin(np.abs(t_det_h_found-tc(t_stamps)))
        BO_lr_det = BO_loss_rate[i_plot, :]
        lr_det = loss_rate[i_plot, :]
        bint_det = bint[i_plot, :]
    
        mask_filled_det = bint_det>2e10
        BO_lr_det[~mask_filled_det] = 0.
        lr_det[~mask_filled_det] = 0.

        acc_loss_rate_woBO_det = acc_loss_rate_woBO[i_plot, :]
        acc_loss_rate_woBO_det[~mask_filled_det] = 0.
        acc_BO_loss_det = acc_BO_loss_rate[i_plot, :]
        acc_BO_loss_det[~mask_filled_det] = 0.
    
        lt_BO_avg = 1./(np.sum(BO_lr_det)/np.sum(bint[i_plot]))/3600.
        lt_other = 1./(np.sum(lr_det-BO_lr_det)/np.sum(bint[i_plot]))/3600.
    
        # plt.fill_between(np.arange(3564),BO_lr_det, color='green', alpha=0.6, label='Burn off')
        # plt.fill_between(np.arange(3564), BO_lr_det, lr_det, color='red', alpha=.6, label='Additional losses')
        # ms.sciy()
        # axd.set_ylabel('Loss rate [p/s]')
        # axd.set_xlabel('25 ns slot')
        # axd.legend(loc='lower right', prop={'size':14})
        # axd.set_xlim(0,3500)
        # axd.set_ylim(bottom=0)
        # axd.grid('on')
        # figd.subplots_adjust(right=.94, left=.1, bottom=.12, top=.86)
        # figd.suptitle('Fill %d SB Loss Rate at %.1fh for B%d\nSB started on %s\n'%(filln, t_det_h_found, beam, tref_string))
        
    
        figd2 = plt.figure(200+beam, figsize=(8*1.4,6))
        figd2.clear()
        figd2.set_facecolor('w')
        axd2 = figd2.add_subplot(111, sharex=axslot, sharey=axd2)
         
        if mode == 'loss_rate' or mode == 'lifetime':
            axd2.plot(np.arange(3564), BO_lr_det*3600*1e-9, 'b', lw=2, alpha=.4, label='Burn-off')
            axd2.plot(np.arange(3564), (lr_det-BO_lr_det)*3600*1e-9, 'red', lw=2, label='Other losses')
            ms.sciy()
            axd2.set_ylabel(r'Loss rate [10$^9$ p/h]')
        elif mode == 'integrated':
            axd2.plot(np.arange(3564), acc_BO_loss_det*1e-10, 'b', lw=2, alpha=.4, label='Burn-off')
            axd2.plot(np.arange(3564), acc_loss_rate_woBO_det*1e-10, 'red', lw=2, label='Other losses')
            ms.sciy()
            axd2.set_ylabel(r'Accumulated losses [10$^{10}$ p]')
        
        
        axd2.set_xlabel('25 ns slot')
        axd2.legend(loc='upper right', prop={'size':14})
        axd2.set_xlim(0,3500)
        axd2.set_ylim(bottom=-.2, top=5.)
        axd2.grid('on')
        figd2.subplots_adjust(right=.94, left=.1, bottom=.12, top=.86)
        figd2.suptitle('Fill %d Loss Rates at %.1fh for B%d\nFT started on %s\n'%(filln, t_det_h_found, beam, tref_string))
    
    if slotrange is not None:
        axslot.set_xlim(slotrange)
    
    if T_h_range_colorplot is not None:
        axgr.set_ylim(T_h_range_colorplot)
    
    
    def savefigures(folder, beam, filln, t_det_h_found, tag=''):
        fig.savefig(folder+'/fill%d_%s_beam%s_evol.png'%(filln, tag, beam), dpi=200)
        figd2.savefig(folder+'/fill%d_%s_beam%s_detals_at_%.2fh.png'%(filln, tag, beam, t_det_h_found), dpi=200)
    
    
    if outp_folder is not None:
        savefigures(folder=outp_folder, beam=beam, filln=filln, t_det_h_found=t_det_h_found)

plt.show()
