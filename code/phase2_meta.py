# Header fields of the Phase 2 recordings (start time, transmitter and receiver positions, sample rate) used by
# phase2_analyze.py and catalog_association.py, written to results/phase2_meta.json. Historical configuration:
# the transmitter position is the one written in each file, as in the original search (README, Configurations).
import os as _os; _os.environ['BRAMS_TX_HEADER']='1'
import os, glob, json, bwav
D=os.environ.get('P2DIR','../data/raw_wav_phase2/')
sel=json.load(open('../data/planning/phase2_files.json'))
path={os.path.basename(p)[11:31]:p for d in D.split(':') if d for p in glob.glob(os.path.join(d,'*.wav'))}
K=('start_us','beacon_lat','beacon_long','beacon_alt','ant_lat','ant_long','ant_alt','samplerate')
out={}
for f in sorted(x for x in sel if x in path):
    md=bwav.read(path[f])[0]; out[f]={k:(int(md[k]) if k=='start_us' else float(md[k])) for k in K}
json.dump(out,open('../results/phase2_meta.json','w'))
print(len(out),'recordings')
