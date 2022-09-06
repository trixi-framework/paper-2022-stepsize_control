
using LaTeXStrings
import PyCall
import PyPlot; plt = PyPlot

# line cyclers adapted to colorblind people
let
  cycler = PyCall.pyimport("cycler").cycler

  global const line_cycler = (
    cycler(color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]) +
    cycler(linestyle=["-", "--", "-.", ":", "-", "--", "-."]))

  global const marker_cycler = (
    cycler(color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]) +
    cycler(linestyle=["none", "none", "none", "none", "none", "none", "none"]) +
    cycler(marker=["4", "2", "3", "1", "+", "x", "."]))

  # matplotlib's standard cycler
  global const standard_cycler = cycler("color",
    ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
    "#7f7f7f", "#bcbd22", "#17becf"])

  plt.rc("axes", prop_cycle=line_cycler)

  # plt.rc("text", usetex=true)
  # plt.rc("text.latex", preamble="\\usepackage{newpxtext}\\usepackage{newpxmath}\\usepackage{commath}\\usepackage{mathtools}")
  plt.rc("font", family="serif", size=18.)
  plt.rc("savefig", dpi=200)
  plt.rc("legend", loc="best", fontsize="medium", fancybox=true, framealpha=0.5)
  plt.rc("lines", linewidth=2.5, markersize=10, markeredgewidth=2.5)

#   SMALL_SIZE = 8
#   MEDIUM_SIZE = 10
#   BIGGER_SIZE = 12
#   plt.rc("font", size=SMALL_SIZE)         # controls default text sizes
#   plt.rc("axes", titlesize=SMALL_SIZE)    # fontsize of the axes title
#   plt.rc("axes", labelsize=MEDIUM_SIZE)   # fontsize of the x and y labels
#   plt.rc("xtick", labelsize=SMALL_SIZE)   # fontsize of the tick labels
#   plt.rc("ytick", labelsize=SMALL_SIZE)   # fontsize of the tick labels
#   plt.rc("legend", fontsize=SMALL_SIZE)   # legend fontsize
#   plt.rc("figure", titlesize=BIGGER_SIZE) # fontsize of the figure title
end
