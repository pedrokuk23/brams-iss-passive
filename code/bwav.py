import struct, numpy as np
FMT="<H 2d 2Q 5d 2H 5d 6s 6s 6s 227s H 3B H"
def read(path):
    b=open(path,'rb').read(); i=12; ch={}
    while i<len(b)-8:
        cid=b[i:i+4].decode('latin1'); n=struct.unpack('<I',b[i+4:i+8])[0]; ch[cid]=b[i+8:i+8+n]; i+=8+n+(n&1)
    h=struct.unpack(FMT,ch['BRA1'][:struct.calcsize(FMT)])
    keys="version samplerate lo_freq start_us pps_count beacon_lat beacon_long beacon_alt beacon_freq beacon_power beacon_polar ant_id ant_lat ant_long ant_alt ant_az ant_el beacon_code observer_code station_code".split()
    md=dict(zip(keys,h[:len(keys)]))
    x=np.frombuffer(ch['data'],np.int16).astype(float)
    pps=np.frombuffer(ch.get('BRA2',b''),dtype=np.uint64).reshape(-1,2) if 'BRA2' in ch else None
    return md,x,pps
