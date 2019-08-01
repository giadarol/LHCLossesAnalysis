

fname_lpc_table = 'fill_list_from_lpc.txt'

with open(fname_lpc_table, 'r') as fid:
    lines = fid.readlines()
    

filln_list = []

for line in lines:
    fn = int(line.split()[0])
    fs = line.split()[-1]
    
    if '2556b_2544' in fs:
        filln_list.append(fn)
    
