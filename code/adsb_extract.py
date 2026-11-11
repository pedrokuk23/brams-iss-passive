"""Extrai do arquivo diário do adsb.lol (globe_history) as trajetórias dentro de uma caixa lat/lon e de uma janela de tempo.

Uso (o arquivo do dia vem dividido em partes .tar.aa, .tar.ab ...):
    cat v2026.11.03-planes-readsb-prod-0.tar.a* | python3 adsb_extract.py 2026-11-03T09:55 2026-11-03T10:35 ../data/adsb/adsb_be_nov03.jsonl
Fonte: https://github.com/adsblol/globe_history_2026 (licença ODbL 1.0)
Saída: uma linha JSON por aeronave {icao, r, t, pts:[[t_unix, lat, lon, alt_baro, gs, track, alt_geom], ...]}
"""
import sys, tarfile, gzip, json, datetime as dt
lat0,lat1,lon0,lon1=48.8,52.0,2.0,7.2          # Bélgica e arredores
t0=dt.datetime.fromisoformat(sys.argv[1]).replace(tzinfo=dt.timezone.utc).timestamp()
t1=dt.datetime.fromisoformat(sys.argv[2]).replace(tzinfo=dt.timezone.utc).timestamp()
out=open(sys.argv[3],'w'); n=0; kept=0
tf=tarfile.open(fileobj=sys.stdin.buffer,mode='r|')
for m in tf:
    if not m.isfile() or 'trace_full' not in m.name: continue
    n+=1
    raw=tf.extractfile(m).read()
    try:
        d=json.loads(gzip.decompress(raw) if raw[:2]==b'\x1f\x8b' else raw)
    except Exception: continue
    base=d.get('timestamp',0); pts=[]
    for p in d.get('trace',[]):
        ts=base+p[0]
        if t0<=ts<=t1 and p[1] is not None and lat0<=p[1]<=lat1 and lon0<=p[2]<=lon1:
            pts.append([ts,p[1],p[2],p[3],p[4],p[5],p[10] if len(p)>10 else None])
    if pts:
        kept+=1
        out.write(json.dumps({'icao':d.get('icao'),'r':d.get('r'),'t':d.get('t'),'pts':pts})+'\n')
    if n%50000==0: print(n,kept,flush=True,file=sys.stderr)
print('total',n,'kept',kept,file=sys.stderr)
