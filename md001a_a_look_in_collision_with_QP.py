import sys, os
BIN = os.path.expanduser("../LHC_fullRun2_analysis_scripts/")
sys.path.append(BIN)

import LHCMeasurementTools.mystyle as ms
import LHCMeasurementTools.TimestampHelpers as th

import pickle
import gzip
import time
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator

import pjlsa
import pytimber

import numpy as np


t_start_string = '2018_07_24 00:01:00'
t_stop_string = '2018_07_24 02:00:00'

t_h_detail = 0.5

physics_BP_end = "PHYSICS-6.5TeV-30cm-120s-2018_V1@120_[END]"
physics_BP = "PHYSICS-6.5TeV-30cm-120s-2018_V1"



def with_empty_slots(lifet_woBO_h, slots):
    
    mat = lifet_woBO_h.copy()
    outp = np.zeros((mat.shape[0], 3564), dtype=float)
    outp[:, slots]=mat

    return outp



plt.close('all')
ms.mystyle_arial(fontsz=14, dist_tick_lab=5)

t_start_SB = th.localtime2unixstamp(t_start_string)
t_stop_SB = th.localtime2unixstamp(t_stop_string)

time_conv = th.TimeConverter(time_in='h', t_plot_tick_h=None, t_ref=np.float_(t_start_SB))
tc = time_conv.from_unix


lsa = pjlsa.LSAClient()

QP_offsets = {
'LHCBEAM1/QPH': 15.-7.95981,
'LHCBEAM1/QPV': 15.-20.1143,
'LHCBEAM2/QPH': 15.-8.51945,
'LHCBEAM2/QPV': 15.-18.32254,
}

chromaTrims_end = lsa.getTrims(parameter=['LHCBEAM1/QPH','LHCBEAM1/QPV','LHCBEAM2/QPH','LHCBEAM2/QPV'], 
                             beamprocess=physics_BP_end,  
                             start=t_start_SB-30*60, end=t_stop_SB)

ldb = pytimber.LoggingDB(source='ldb')
data_timb = ldb.getScaled(['RPMBB.RR17.ROD.A12B1:I_MEAS', 'RPMBB.RR17.ROD.A12B2:I_MEAS'],
                    t_start_SB, t_stop_SB, scaleAlgorithm='AVG', scaleInterval='MINUTE',scaleSize='1')
                    
print 'Start downloading intensity...'
data_timb.update(ldb.get([
            'LHC.BCTFR.A6R4.B1:BUNCH_INTENSITY',
            'LHC.BCTFR.A6R4.B2:BUNCH_INTENSITY'], t_start_SB, t_stop_SB))
print 'Downloaded Intensity'

print 'Start downloading xang'
data_timb.update(ldb.getScaled(['LHC.RUNCONFIG:IP1-XING-V-MURAD'],
                    t_start_SB, t_stop_SB, scaleAlgorithm='REPEAT', scaleInterval='MINUTE',scaleSize='1'))
print 'Downloaded xang'


print 'Start downloading betastar'
data_timb.update(ldb.getScaled(['HX:BETASTAR_IP1'],
                    t_start_SB, t_stop_SB, scaleAlgorithm='REPEAT', scaleInterval='MINUTE',scaleSize='1'))
print 'Downloaded betastar'


# plt.figure(2)
# plt.pcolormesh(fill_data['b_inten_interp_coll'][beam])
# plt.colorbar(label='Bunch intensity')

figlist = []

axslot = None
axd = None
axd2 = None
for beam in [1,2]:
    
    fbct = data_timb['LHC.BCTFR.A6R4.B%d:BUNCH_INTENSITY'%beam][1]
    t_stamps = data_timb['LHC.BCTFR.A6R4.B%d:BUNCH_INTENSITY'%beam][0]

    slots = np.where(fbct[0]>5e10)[0]
    bint = np.take(np.array(fbct), slots, axis=1)
    

    t_mat = np.dot(np.ones((len(slots), 1)), np.atleast_2d(t_stamps)).T
    loss_rate = -np.diff(bint, axis=0)/np.diff(t_mat, axis=0)
    lifet_h = 1/(loss_rate/bint[:-1,:])/3600.
    lifet_h[lifet_h<0]= 200


    loss_rate_woBO = loss_rate

    lifet_woBO_h = 1/(loss_rate_woBO/bint[:-1,:])/3600.
    lifet_woBO_h[lifet_woBO_h<0]= 200
    lifet_woBO_h[lifet_woBO_h>200]= 200

    # plt.figure(3)
    # plt.pcolormesh(lifet_h, cmap=cm.jet_r)
    # plt.colorbar()

    fig = plt.figure(beam, figsize=(8*1.8,6*1.3))
    fig.set_facecolor('w')
    axlt = plt.subplot2grid(shape=(5, 5), loc=(0, 1), colspan=3, rowspan=4, sharex=axslot)
    lifet_woBO_h_allslots = with_empty_slots(lifet_woBO_h, slots)

    cc=axlt.pcolormesh(np.arange(3564), tc(t_stamps), lifet_woBO_h_allslots, 
        cmap=cm.jet_r, vmin=0, vmax=200)
    axcb = plt.subplot2grid(shape=(5, 5), loc=(4, 1), colspan=3, rowspan=1)
    plt.colorbar(cc, cax=axcb, label='Lifetime [h]', orientation='horizontal')
    
    axlt.set_xlabel('25ns slot')
    axslot = axlt
    axslot. set_xlim(0, 3500)

    ax1 = plt.subplot2grid(shape=(5, 5), loc=(0, 0), colspan=1, rowspan=4, sharey=axlt)
    
    var = 'LHC.RUNCONFIG:IP1-XING-V-MURAD';ax1.step(data_timb[var][1], tc(data_timb[var][0]), lw=2.)
    ax1.set_xlim(120, 170)
    ax1.set_xlabel('Half crossing angle [urad]')
    ax1.set_ylabel('Time [h]')
    ax1.xaxis.label.set_color('b')
    ax1.tick_params(axis='x', colors='b')
    ax1.grid('on')

    axlumi = ax1.twiny()
    var = 'HX:BETASTAR_IP1';axlumi.step(data_timb[var][1], tc(data_timb[var][0]), lw=2., color='r')
    axlumi.set_xlabel('Beta* [m]')
    axlumi.xaxis.label.set_color('r')
    axlumi.tick_params(axis='x', colors='r')
    axlumi.xaxis.set_major_locator(MaxNLocator(5))

    axqp = plt.subplot2grid(shape=(5, 5), loc=(0, 4), colspan=1, rowspan=4, sharey=axlt)
    for plane, styl in zip(['H', 'V'], ['--', '-']):
        if len(chromaTrims_end)>0:
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

    tref_string=time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime(t_stamps[0]))
    fig.suptitle('Lifetime B%d\nTest started on %s\n'%(beam, tref_string))

    figlist.append(fig)


    if t_h_detail<tc(t_stamps)[-1]:

        i_plot = np.argmin(np.abs(t_h_detail-tc(t_stamps)))
        lr_det = with_empty_slots(loss_rate, slots)[i_plot, :]

        lt_other = 1./(np.sum(lr_det)/np.sum(bint[i_plot]))/3600.


        figd2 = plt.figure(200+beam, figsize=(8*1.4,6))
        figd2.set_facecolor('w')
        axd2 = figd2.add_subplot(111, sharex=axslot, sharey=axd2)

        axd2.plot(np.arange(3564), lr_det, 'r', lw=2, label='Other losses (tau=%.1fh)'%lt_other)
        ms.sciy()
        axd2.set_ylabel('Loss rate [p/s]')
        axd2.set_xlabel('25 ns slot')
        axd2.legend(loc='upper right', prop={'size':14})
        axd2.set_xlim(0,3500)
        axd2.set_ylim(bottom=0)
        axd2.grid('on')
        figd2.subplots_adjust(right=.94, left=.1, bottom=.12, top=.86)
        figd2.suptitle('Loss Rates at %.1fh for B%d\nTest started on %s\n'%(t_h_detail, beam, tref_string))
        figlist.append(figd2)

    # plt.figure(40)
    # plt.pcolormesh(BO_loss_rate/loss_rate, cmap=cm.jet_r)
    # plt.colorbar()

    # plt.figure(5)
    # plt.pcolormesh(tot_lumi, cmap=cm.jet)
    # plt.colorbar()

#~ ff = '/eos/user/g/giadarol/temp/20180513_lossesbbb'
#~ [fg.savefig(ff+'/'+fg._suptitle.get_text().split('\n')[0].replace(' ','_')+'.png', dpi=200) for fg in figlist]
#~ axslot.set_xlim(770, 1270);
#~ [fg.savefig(ff+'/'+fg._suptitle.get_text().split('\n')[0].replace(' ','_')+'_zoom.png', dpi=200) for fg in figlist]


plt.show()

# To save
# for fg in figlist: fg.savefig(ff+'/'+fg._suptitle.get_text().split('\n')[0].replace(' ','_')+'.png', dpi=200)

