import os, struct, numpy as np
# Dourbes transmitter: one position for every receiver. The BRA1 headers of the version-5 files (Humain, Uccle, Maubeuge)
# give 50.09756 N, 4.58852 E, 236 m; the version-4 headers differ (Oudsbergen by 37 m in altitude, Ottignies by 276 m,
# with coordinates rounded to 1e-4 deg). The header values are kept as beacon_*_header; BRAMS_TX_HEADER=1 uses them.
DOURBES=tuple(float(v) for v in os.environ['BRAMS_TX'].split(',')) if os.environ.get('BRAMS_TX') else (50.09756,4.58852,236.0)
FMT="<H 2d 2Q 5d 2H 5d 6s 6s 6s 227s H 3B H"
def read(path):
    b=open(path,'rb').read(); i=12; ch={}
    while i<len(b)-8:
        cid=b[i:i+4].decode('latin1'); n=struct.unpack('<I',b[i+4:i+8])[0]; ch[cid]=b[i+8:i+8+n]; i+=8+n+(n&1)
    h=struct.unpack(FMT,ch['BRA1'][:struct.calcsize(FMT)])
    keys="version samplerate lo_freq start_us pps_count beacon_lat beacon_long beacon_alt beacon_freq beacon_power beacon_polar ant_id ant_lat ant_long ant_alt ant_az ant_el beacon_code observer_code station_code".split()
    md=dict(zip(keys,h[:len(keys)]))
    for k,v in zip(('beacon_lat','beacon_long','beacon_alt'),DOURBES):
        md[k+'_header']=md[k]
        if os.environ.get('BRAMS_TX_HEADER')!='1': md[k]=v
    x=np.frombuffer(ch['data'],np.int16).astype(float)
    pps=np.frombuffer(ch.get('BRA2',b''),dtype=np.uint64).reshape(-1,2) if 'BRA2' in ch else None
    return md,x,pps
