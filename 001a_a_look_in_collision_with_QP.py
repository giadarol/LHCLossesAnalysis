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

t_h_detail = 1

filln = 6641
filln = 6640
filln = 6648

filln = 6650
filln = 6659
filln = 6662
# filln = 6663


t_h_detail = 1

folder = '/afs/cern.ch/user/l/lumimod/lumimod/2018/procdata/fill_%d'%filln

with gzip.open(folder+'/fill_%d.pkl.gz'%filln, 'r') as fid:
    fill_data = pickle.load(fid)
print('Loaded fill data.')    

with gzip.open(folder+'/fill_%d_lumi_meas.pkl.gz'%filln, 'r') as fid:
    lumi_data = pickle.load(fid)
print('Loaded lumi data.') 

N_traces = len(fill_data['time_range']) 

n_dec = 10 

beam =1

sigma_m2 = 80e-3*1e-28

def with_empty_slots(lifet_woBO_h, slots):
    
    mat = lifet_woBO_h.copy()
    outp = np.zeros((mat.shape[0], 3564), dtype=float)
    outp[:, slots]=mat

    return outp



plt.close('all')
ms.mystyle_arial(fontsz=14, dist_tick_lab=5)
# Instantaneous lumi is in m^-2 s^-1

t_stamps = fill_data['time_range']
time_conv = th.TimeConverter(time_in='h', t_plot_tick_h=None, t_ref=t_stamps[0])
tc = time_conv.from_unix

t_start_SB = t_stamps[0]
t_stop_SB = t_stamps[-1]

physics_BP_end = "PHYSICS-6.5TeV-30cm-120s-2018_V1@120_[END]"
physics_BP = "PHYSICS-6.5TeV-30cm-120s-2018_V1"
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


# plt.figure(2)
# plt.pcolormesh(fill_data['b_inten_interp_coll'][beam])
# plt.colorbar(label='Bunch intensity')

figlist = []

axslot = None
axd = None

for beam in [1,2]:

    slots = fill_data['slots_filled_coll'][beam]
    bint = fill_data['b_inten_interp_coll'][beam]

    t_mat = np.dot(np.ones((len(slots), 1)), np.atleast_2d(t_stamps)).T
    loss_rate = -np.diff(bint, axis=0)/np.diff(t_mat, axis=0)
    lifet_h = 1/(loss_rate/bint[:-1,:])/3600.
    lifet_h[lifet_h<0]= 200


    tot_lumi = (lumi_data['ATLAS']['bunch_lumi']+lumi_data['CMS']['bunch_lumi'])[:-1, :]
    BO_loss_rate = tot_lumi*sigma_m2

    loss_rate_woBO = loss_rate-BO_loss_rate

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
        cmap=cm.jet_r, vmin=0, vmax=120)
    axcb = plt.subplot2grid(shape=(5, 5), loc=(4, 1), colspan=3, rowspan=1)
    plt.colorbar(cc, cax=axcb, label='Lifetime (BO corrected) [h]', orientation='horizontal')
    
    axlt.set_xlabel('25ns slot')
    axslot = axlt


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

    tref_string=time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime(t_stamps[0]))
    fig.suptitle('Fill %d SB BurnOff corrected lifetime B%d\nSB started on %s\n'%(filln, beam, tref_string))

    figlist.append(fig)

    figd = plt.figure(100+beam, figsize=(8*1.4,6))
    figd.set_facecolor('w')
    axd = figd.add_subplot(111, sharex=axslot, sharey=axd)
    i_plot = np.argmin(np.abs(t_h_detail-tc(t_stamps)))
    BO_lr_det = with_empty_slots(BO_loss_rate, slots)[i_plot, :]
    lr_det = with_empty_slots(loss_rate, slots)[i_plot, :]
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
    figd.suptitle('Fill %d SB Loss Rates at %.1fh for B%d\nSB started on %s\n'%(filln, t_h_detail, beam, tref_string))
    
    figlist.append(figd)

    # plt.figure(40)
    # plt.pcolormesh(BO_loss_rate/loss_rate, cmap=cm.jet_r)
    # plt.colorbar()

    # plt.figure(5)
    # plt.pcolormesh(tot_lumi, cmap=cm.jet)
    # plt.colorbar()



plt.show()

# To save
# for fg in figlist: fg.savefig(ff+'/'+fg._suptitle.get_text().split('\n')[0].replace(' ','_')+'.png', dpi=200)

