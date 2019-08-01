import sys, os
BIN = os.path.expanduser("../LHC_fullRun2_analysis_scripts/")
sys.path.append(BIN)

import LHCMeasurementTools.mystyle as ms
import numpy as np

beam = 2

xang_obs_urad = 140.


filln_list_str = '''
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
 6909
 6911
 6912
 6919
 6921
 6923
 6924
 6925
 6927
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
 7018
 7020
 7024
 7026
 7031
 7033
 7035
'''



filln_list = map(int, filln_list_str.split())

lifet_h = {1:[], 2:[]}

for filln in filln_list:
    for err_xang in [0, 1., -1.]:
        found = False
        try:
            with open('avglifet_SB/avg_lifet_fill%d_beam%d_%.1furad.txt'%(filln, beam, xang_obs_urad+err_xang), 'r') as fid:
                line = fid.read()
                lifet_h[beam].append(float(line.split()[0]))
                found = True
                break
        except IOError:
            pass
    if not found:
        print('skipped %d'%filln)
        lifet_h[beam].append(np.nan)


import matplotlib.pyplot as plt
plt.close('all')
ms.mystyle_arial(fontsz=14, dist_tick_lab=5)

fig = plt.figure(1, figsize=(1.7*8,6))
fig.set_facecolor('w')
sp1=plt.subplot(2,1,1)



sp2=plt.subplot(2,1,2, sharex=sp1)

plt.plot(np.arange(len(lifet_h[beam])), np.array(lifet_h[beam]), '.', markersize=10)

plt.xticks(np.arange(len(filln_list)), filln_list, rotation='vertical')

fig.subplots_adjust(hspace=.33)
fig.suptitle('Beam %s, at %.1f urad'%(beam, xang_obs_urad))

for sp in [sp1, sp2]:
	sp.grid('on')

sp1.set_ylabel('Octupole current [A]')
sp2.set_ylabel('Average bunch lifetime [h]')


figl = plt.figure(2, figsize=(1.7*8,6))
figl.set_facecolor('w')
spl1=plt.subplot(2,1,1)



spl2=plt.subplot(2,1,2, sharex=spl1)

plt.plot(np.arange(len(lifet_h[beam])), np.array(1./np.array(lifet_h[beam])*100), '.', markersize=10)

plt.xticks(np.arange(len(filln_list)), filln_list, rotation='vertical')

figl.subplots_adjust(hspace=.33)
figl.suptitle('Beam %s, at beta_st %.2f cm\nLifetime measured over 3 min.'%(beam, xang_obs_urad))

for sp in [spl1, spl2]:
	sp.grid('on')

spl1.set_ylabel('Octupole current [A]')
spl2.set_ylabel('Average loss rate [%/h]')


for ff in [fig, figl]:
	ff.subplots_adjust(left=.06, right=.99)

plt.show()



