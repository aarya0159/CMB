import json, os, sys, time, numpy as np, core
W,H,NS = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
BUDGET = float(sys.argv[4]) if len(sys.argv)>4 else 130.0
path = "cell_%dx%d_%d.json" % (W,H,NS)
vals = json.load(open(path))["dT"] if os.path.exists(path) else []
t0 = time.time()
while len(vals) < NS and time.time()-t0 < BUDGET:
    vals.append(core.field_std(W,H,384,len(vals)+1,5.0))
json.dump({"W":W,"H":H,"n_target":NS,"dT":vals}, open(path,"w"))
v=np.array(vals)
print("%4dx%-4d %d/%d done  mean=%.6f cov=%.6f  [%.0fs]"%(W,H,len(vals),NS,v.mean(),v.std(ddof=1)/v.mean(),time.time()-t0))
