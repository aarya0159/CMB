import json, math, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, FixedLocator, FixedFormatter

plt.rcParams.update({"font.size":8,"font.family":"serif","axes.linewidth":0.7,
                     "xtick.direction":"in","ytick.direction":"in",
                     "xtick.top":True,"ytick.right":True,"savefig.dpi":400})
BLUE,RED,GREY = "#1f4e79","#c0392b","#7f7f7f"
wa=json.load(open("wave_analysis.json")); wd=json.load(open("width_analysis.json")); ar=json.load(open("area_analysis.json"))

# authored at the width it is PLACED (NHSJS single column ~ 8.5 cm per panel)
fig,axes = plt.subplots(1,3, figsize=(7.2,2.35))

# (a) CoV vs N
ax=axes[0]
N=np.array([r["N"] for r in wa["rows"]],float)
y=np.array([r["cov"] for r in wa["rows"]]); e=np.array([r["boot"] for r in wa["rows"]])
xs=np.logspace(math.log10(4),math.log10(600),200)
ax.plot(xs,0.806*xs**-0.5,"--",color=RED,lw=1.1,label=r"predicted $0.806\,N^{-1/2}$")
ax.plot(xs,np.sqrt(wa["A"]/xs+wa["B"]),"-",color=BLUE,lw=1.2,label=r"model $\sqrt{A/N+B}$")
ax.axhline(wa["floor"],ls=":",color=GREY,lw=1.0)
ax.errorbar(N,y,yerr=e,fmt="o",ms=3.4,color=BLUE,ecolor=BLUE,elinewidth=0.9,capsize=1.6,
            mfc="white",mew=0.9,label="measured",zorder=5)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("number of waves  $N$"); ax.set_ylabel(r"seed-to-seed scatter  $\sigma_{\Delta T}/\langle\Delta T\rangle$")
ax.text(0.97,0.14,r"finite-map floor $\sqrt{B}$",color=GREY,fontsize=6.3,
        transform=ax.transAxes,ha="right")
ax.legend(frameon=False,fontsize=6.3,loc="upper right")
ax.set_title("(a)  scatter against wave count",fontsize=8,pad=4)

# (b) floor vs AREA, all seven maps
ax=axes[1]
rows=ar["rows"]; A=np.array([r["area"] for r in rows],float)
F=np.array([r["floor"] for r in rows]); E=np.array([r["err"] for r in rows])
xs=np.logspace(math.log10(1.8e4),math.log10(4.6e5),100)
ax.plot(xs,ar["C"]/np.sqrt(xs),"--",color=RED,lw=1.1,
        label=r"$%.1f\,/\sqrt{\mathrm{area}}$"%ar["C"])
ax.errorbar(A,F,yerr=E,fmt="o",ms=3.4,color=BLUE,ecolor=BLUE,elinewidth=0.9,capsize=1.6,
            mfc="white",mew=0.9,zorder=5)
for r in rows:
    if r["area"]==184320:
        ax.annotate("%d$\\times$%d"%(r["W"],r["H"]),(r["area"],r["floor"]),
                    textcoords="offset points",xytext=(-7,6) if r["W"]==384 else (-7,-12),fontsize=5.8,color=GREY,ha="right")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"map area  $W\!\times\!H$  (pixels$^2$)"); ax.set_ylabel(r"floor  $\sqrt{B}$")
ax.legend(frameon=False,fontsize=6.3,loc="upper right")
ax.set_title("(b)  the floor is set by area alone",fontsize=8,pad=4)
ax.text(0.03,0.06,"slope $%.3f\\pm%.3f$\n(predicted $-0.500$)"%(ar["area_slope"],ar["area_slope_err"]),
        transform=ax.transAxes,fontsize=6.3)

# (c) the three geometries
ax=axes[2]
cols={"both":BLUE,"width":"#2e7d32","height":"#8e44ad"}
labs={"both":r"both dims: $-1.00$","width":r"width only: $-0.50$","height":r"height only: $-0.50$"}
for g in ("both","width","height"):
    d=wd[g]; x=np.array(d["x"],float); f=np.array(d["floor"]); e=np.array(d["floor_err"])
    ax.errorbar(x,f,yerr=e,fmt="o-",ms=3.2,lw=0.9,color=cols[g],ecolor=cols[g],
                elinewidth=0.8,capsize=1.5,mfc="white",mew=0.9,
                label="%s\nfitted $%.3f\\pm%.3f$"%(labs[g],d["slope"],d["slope_err"]))
ax.set_xscale("log"); ax.set_yscale("log")
ax.xaxis.set_minor_formatter(NullFormatter())
ax.xaxis.set_major_locator(FixedLocator([120,192,384,768]))
ax.xaxis.set_major_formatter(FixedFormatter(["120","192","384","768"]))
ax.xaxis.set_minor_locator(FixedLocator([240,480]))
ax.set_xlabel("scaled dimension (pixels)"); ax.set_ylabel(r"floor  $\sqrt{B}$")
ax.legend(frameon=False,fontsize=5.6,loc="lower left")
ax.set_title("(c)  each geometry on its own law",fontsize=8,pad=4)

for a in axes: a.tick_params(labelsize=6.8)
fig.tight_layout(pad=0.5,w_pad=1.4)
fig.savefig("figure1_revised.png",bbox_inches="tight")
fig.savefig("figure1_revised.pdf",bbox_inches="tight")
print("wrote figure1_revised.png / .pdf")
