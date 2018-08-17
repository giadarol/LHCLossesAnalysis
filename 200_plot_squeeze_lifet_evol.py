import sys, os
BIN = os.path.expanduser("../LHC_fullRun2_analysis_scripts/")
sys.path.append(BIN)

import LHCMeasurementTools.mystyle as ms
import numpy as np

beam = 1

beta_obs_cm = 33.

flag_avgs = False
first_fill_after_TS = 6860

filln_list_str = '''6639
6640
6641
6642
6643
6645
6646
6648
6650
6654
6659
6662
6663
6666
6672
6674
6675
6677
6681
6683
6696
6700
6702
6709
6710
6711
6712
6714
6719
6724
6729
6731
6733
6737
6738
6740
6741
6744
6749
6751
6752
6755
6757
6759
6761
6762
6763
6768
6770
6772
6773
6774
6776
6778
6860
6874
6904
6909
6911
6912
6919
6921
6923
6924
6925
6931
6940
6942
6944
6946
6953
6956
6957
6960
6961
7006
7008
7013
7017
7020
7024
7031
7033
7035
'''



filln_list = map(int, filln_list_str.split())

oct_curr = {1:[], 2:[]}
lifet_h = {1:[], 2:[]}

for filln in filln_list:

        with open('avglifet/avg_lifet_fill%d_beam%d_betast%.1fcm.txt'%(filln, beam, beta_obs_cm), 'r') as fid:
            line = fid.read()
            oct_curr[beam].append(float(line.split()[0]))
            lifet_h[beam].append(float(line.split()[1]))


import matplotlib.pyplot as plt
plt.close('all')
ms.mystyle_arial(fontsz=14, dist_tick_lab=5)

fig = plt.figure(1, figsize=(1.7*8,6))
fig.set_facecolor('w')
sp1=plt.subplot(2,1,1)

plt.plot(np.arange(len(oct_curr[beam])), np.array(oct_curr[beam]), '.', markersize=10)
plt.xticks(np.arange(len(filln_list)), filln_list, rotation='vertical')


sp2=plt.subplot(2,1,2, sharex=sp1)

plt.plot(np.arange(len(oct_curr[beam])), np.array(lifet_h[beam]), '.', markersize=10)

plt.xticks(np.arange(len(filln_list)), filln_list, rotation='vertical')

fig.subplots_adjust(hspace=.33)
fig.suptitle('Beam %s, at beta_st %.2f cm\nLifetime measured over 3 min.'%(beam, beta_obs_cm))

for sp in [sp1, sp2]:
	sp.grid('on')

sp1.set_ylabel('Octupole current [A]')
sp2.set_ylabel('Average bunch lifetime [h]')


figl = plt.figure(2, figsize=(1.7*8,6))
figl.set_facecolor('w')
spl1=plt.subplot(2,1,1)

plt.plot(np.arange(len(oct_curr[beam])), np.array(oct_curr[beam]), '.', markersize=10)
plt.xticks(np.arange(len(filln_list)), filln_list, rotation='vertical')


spl2=plt.subplot(2,1,2, sharex=spl1)

plt.plot(np.arange(len(oct_curr[beam])), np.array(1./np.array(lifet_h[beam])*100), '.', markersize=10)

plt.xticks(np.arange(len(filln_list)), filln_list, rotation='vertical')

figl.subplots_adjust(hspace=.33)
figl.suptitle('Beam %s, at beta_st %.2f cm\nLifetime measured over 3 min.'%(beam, beta_obs_cm))

for sp in [spl1, spl2]:
	sp.grid('on')

spl1.set_ylabel('Octupole current [A]')
spl2.set_ylabel('Average loss rate [%/h]')

if flag_avgs:
	mask_before = np.array(filln_list)<first_fill_after_TS
	sp2.plot([0, np.sum(mask_before)], np.array([1.,1.])*np.mean(np.array(lifet_h[beam])[mask_before]))
	sp2.plot([np.sum(mask_before), len(mask_before)], np.array([1.,1.])*np.mean(np.array(lifet_h[beam])[~mask_before]))

for ff in [fig, figl]:
	ff.subplots_adjust(left=.06, right=.99)

plt.show()



