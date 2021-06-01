import os 
import sys 
from Bio import SeqIO
from Bio import pairwise2
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd 
import logomaker
matplotlib.rcParams['font.family']       = 'sans-serif'
matplotlib.rcParams['font.sans-serif']   = ["Arial","DejaVu Sans","Lucida Grande","Verdana"]
matplotlib.rcParams['figure.figsize']    = [3,3]
matplotlib.rcParams['font.size']         = 10
matplotlib.rcParams["axes.labelcolor"]   = "#000000"
matplotlib.rcParams["axes.linewidth"]    = 1.0 
matplotlib.rcParams["xtick.major.width"] = 1.0
matplotlib.rcParams["ytick.major.width"] = 1.0

_atgc_dict = {0:"A", 1:"T", 2:"G", 3:"C"}

def abi_to_dict(filename):
    record   = SeqIO.read(filename,'abi')
    abi_data = {"conf":[],
                "channel":{"A":[],
                           "T":[],
                           "G":[],
                           "C":[],
                          },
                "_channel":{"A":[],
                            "T":[],
                            "G":[],
                            "C":[],
                          }
                }
    for i, (pos, conf) in enumerate(zip(record.annotations['abif_raw']["PLOC1"], record.annotations['abif_raw']["PCON1"])):
        if pos > 4 and pos < len(record.annotations['abif_raw']["DATA9"])-5: 
            abi_data["conf"].append(conf)
            abi_data["channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos])
            abi_data["channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos])
            abi_data["channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos])
            abi_data["channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos])

            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos-5])
            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos-3])
            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos])
            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos+3])
            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos+5])
            
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos-5])
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos-3])
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos])
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos+3])
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos+5])
            
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos-5])
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos-3])
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos])
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos+3])
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos+5])

            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos-5])
            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos-3])
            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos])
            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos+3])
            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos+5])

    return abi_data 

def generate_consensusseq(abidata):
    consensus_seq = "" 
    
    for values in zip(abidata["channel"]["A"], abidata["channel"]["T"], abidata["channel"]["G"], abidata["channel"]["C"]):
        consensus_seq += _atgc_dict[values.index(max(values))]
     
    return (consensus_seq, consensus_seq.translate(str.maketrans("ATGC","TACG"))[::-1]) 

def generate_pwm(abidata):
    pwm = {"A":[], "T":[], "G":[], "C":[]} 
    for values in zip(abidata["channel"]["A"], abidata["channel"]["T"], abidata["channel"]["G"], abidata["channel"]["C"]):
        v = 100000 / (sum(values)+1) 
        new_values = (v*values[0], v*values[1], v*values[2], v*values[3])
        new_values = list(map(int, new_values))
        
        while sum(new_values) < 100000:
            for i in range(len(new_values)):
                new_values[i] += 1
                if sum(new_values) == 100000:
                    break 
        
        pwm["A"].append(new_values[0])
        pwm["T"].append(new_values[1])
        pwm["G"].append(new_values[2])
        pwm["C"].append(new_values[3])
    
    pwm=pd.DataFrame(pwm)
    return pwm 

def _colorbar(ax, ref, matches=None, char=True, fontsize=10):
    bars = ax.bar(list(range(len(ref))), [0.9] * (len(ref)), width=1.0, edgecolor="#BBBBBB", linewidth=0.5, align="edge",bottom=0.05)
    ax.set_xlim(0,len(ref))
    ax.set_ylim(0,1.00)
    p = 0
    if matches is None:
        for bar, c in zip(bars,ref):
            bar.set_facecolor("w")
            if char == True:
                ax.text(p+0.5,0.45,c,va="center",ha="center",fontsize=fontsize,zorder=100)     
            p += 1
    
    else:
        for m, bar, c in zip(matches, bars,ref):
            if m == -1:
                bar.set_alpha(0.5)
                bar.set_facecolor("r")
                bar.set_edgecolor("#BBBBBB")
                bar.set_linewidth(0.5)
            else:
                bar.set_facecolor("w")
                bar.set_edgecolor("#BBBBBB")
                bar.set_linewidth(0.5)
            
            if char == True:
                ax.text(p+0.5,0.45,c,va="center",ha="center",fontsize=fontsize,zorder=100)     
            p += 1
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    #ax.patch.set_alpha(0.0)
    return ax


def visualize(abidata, template=None, strand=1, fig=None, region="all"):  
    avalues = abidata["_channel"]["A"]
    tvalues = abidata["_channel"]["T"]
    gvalues = abidata["_channel"]["G"]
    cvalues = abidata["_channel"]["C"]
    consensus_seq_set = generate_consensusseq(abidata)
    
    if strand == 1: 
        subject = consensus_seq_set[0]
    
    if strand == -1:
        subject = consensus_seq_set[1] 
        avalues, tvalues, gvalues, cvalues = tuple(map(lambda x:list(reversed(x)), tvalues, avalues, cvalues, gvalues)) 
   
    if template is not None:
        alignments  = pairwise2.align.globalms(template, subject, 2, 0, -10, -1, penalize_end_gaps=False)
        atemplate   = alignments[0][0]
        asubject    = alignments[0][1]

        new_avalues = []
        new_tvalues = [] 
        new_gvalues = []
        new_cvalues = [] 
        
        ts = 0 
        ss = 0 
        for t, s in zip(atemplate, asubject):
            if t != "-" and s == "-": 
                ts += 1
            
            elif t == "-" and s != "-": 
                ss += 1  
            
            elif t == "-" and s == "-":
                ts += 1
                ss += 1

            else:
                break 
        
        te = 0 
        se = 0 
        for t, s in zip(atemplate[::-1], asubject[::-1]):
            if t != "-" and s == "-": 
                te += 1
            elif t == "-" and s != "-": 
                se += 1  
            elif t == "-" and s == "-":
                te += 1
                se += 1
            else:
                break 
       
        pos     = 0 
        matches = [] 
        for t, s in zip(atemplate, asubject):
            if s == "-":
                new_avalues.extend([0,0,0,0,0]) 
                new_tvalues.extend([0,0,0,0,0]) 
                new_gvalues.extend([0,0,0,0,0]) 
                new_cvalues.extend([0,0,0,0,0]) 
            else:
                new_avalues.append(avalues[5*pos]) 
                new_avalues.append(avalues[5*pos+1]) 
                new_avalues.append(avalues[5*pos+2])
                new_avalues.append(avalues[5*pos+3]) 
                new_avalues.append(avalues[5*pos+4])

                new_tvalues.append(tvalues[5*pos])
                new_tvalues.append(tvalues[5*pos+1]) 
                new_tvalues.append(tvalues[5*pos+2])
                new_tvalues.append(tvalues[5*pos+3]) 
                new_tvalues.append(tvalues[5*pos+4])
                
                new_gvalues.append(gvalues[5*pos])
                new_gvalues.append(gvalues[5*pos+1]) 
                new_gvalues.append(gvalues[5*pos+2])
                new_gvalues.append(gvalues[5*pos+3]) 
                new_gvalues.append(gvalues[5*pos+4])

                new_cvalues.append(cvalues[5*pos])
                new_cvalues.append(cvalues[5*pos+1])
                new_cvalues.append(cvalues[5*pos+2])
                new_cvalues.append(cvalues[5*pos+3]) 
                new_cvalues.append(cvalues[5*pos+4])

                pos += 1
            
            if t == s: 
                matches.append(1) 
            else:
                matches.append(-1) 
        
        if region == "all":
            asubject  = asubject[ts:len(asubject)  - te]
            atemplate = atemplate[ts:len(atemplate) - te]
            matches   = matches[ts:len(asubject)   - te] 
            ss = 0
            se = None
            ts = ts-atemplate[:ts].count("-")
            te = -1 * (te-atemplate[-1*te:].count("-"))
            if te == 0:
                te = None
        else:
            asubject  = asubject[ss:len(asubject)  - se]
            atemplate = atemplate[ss:len(atemplate) - se]
            matches   = matches[ss:len(asubject)   - se] 
            ss = ss-asubject[:ss].count("-")
            se = -1 * (se-asubject[-1*se:].count("-"))
            if se == 0:
                se = None
            ts = 0
            te = None

        avalues = new_avalues 
        tvalues = new_tvalues
        gvalues = new_gvalues
        cvalues = new_cvalues
    else:
        asubject  = subject
        atemplate = template 

    if fig is None:
        fig = plt.figure(figsize=(3,3)) 
    
    ax  = fig.add_axes([0, 0.5, 1.0 * len(asubject)/20,  0.3])
    axs = fig.add_axes([0, 0.42, 1.0 * len(asubject)/20, 0.08])
    
    ax.bar(list(range(len(asubject))), abidata["conf"][ss:se], width=1.0, edgecolor="#BBBBBB", linewidth=0.5, facecolor="#F7F7FF", align="edge", zorder=0)
    ax.set_xlim(0,len(asubject))
    
    ax2 = ax.twinx()
    if se is None:
        _se = None
    else:
        _se = 5*se
    positions = list(map(lambda x: (x+0.5)/5, list(range(len(tvalues[5*ss:_se])))))
    ax2.plot(positions, tvalues[5*ss:_se], color="#FC58FE", lw=1, zorder=1)  
    ax2.plot(positions, avalues[5*ss:_se], color="#33CC33", lw=1, zorder=1) 
    ax2.plot(positions, gvalues[5*ss:_se], color="#303030", lw=1, zorder=1)
    ax2.plot(positions, cvalues[5*ss:_se], color="#395CC5", lw=1, zorder=1) 
    ax2.set_ylim(min(list(tvalues[5*ss:_se]) + list(avalues[5*ss:_se]) + list(gvalues[5*ss:_se]) + list(cvalues[5*ss:_se])), 1.01*max(list(tvalues[5*ss:_se]) + list(avalues[5*ss:_se]) + list(gvalues[5*ss:_se]) + list(cvalues[5*ss:_se])))
    #ax2.set_xlim(0,3*len(asubject))

    ax.set_ylabel("Quality")
    ax2.set_ylabel("Peak height") 
    
    ax.spines["top"].set_visible(False) 
    ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False) 
    ax2.spines["top"].set_visible(False) 
    ax2.spines["bottom"].set_visible(False)
    ax2.spines["left"].set_visible(False) 

    axs = _colorbar(axs, asubject, matches=matches, char=True, fontsize=10)
    length = max([len(asubject), len(atemplate)]) 
    if length < 100:
        tick_space = 10
    elif length < 200:
        tick_space = 20
    elif length < 500:
        tick_space = 50
    elif length < 1000:
        tick_space = 50
    else:
        tick_space = 50
    
    num = 0 
    if ss == 0:
        ticks      = [0.5] 
        ticklabels = ["1"] 
    else:
        ticks      = [] 
        ticklabels = []
    
    for letter in asubject:
        if (ss+num+1) % tick_space == 0: 
            ticks.append(num+0.5) 
            ticklabels.append(str(ss+num+1))
        if letter == "-":
            pass 
        else:
            num += 1 

    axs.set_xticks(ticks) 
    axs.set_xticklabels(ticklabels) 
    
    ax.set_xticks([]) 
    ax.set_xticklabels([]) 
    axs.set_ylabel("Consensus", rotation=0, va="top", ha="right") 

    if atemplate is not None:
        axt = fig.add_axes([0, 0.8, 1.0 * len(atemplate)/20, 0.08])
        axt = _colorbar(axt, atemplate, matches=matches, char=True, fontsize=10)
        num = 0 
        if ts == 0:
            ticks      = [0.5] 
            ticklabels = ["1"] 
        else:
            ticks      = [] 
            ticklabels = []
        
        for letter in atemplate:
            if (ts+num+1) % tick_space == 0: 
                ticks.append(ts+num+0.5) 
                ticklabels.append(str(ts+num+1))
            if letter == "-":
                pass 
            else:
                num += 1 
        axt.xaxis.tick_top() 
        axt.set_xticks(ticks) 
        axt.set_xticklabels(ticklabels) 
        
        axt.set_ylabel("Template", rotation=0, va="center", ha="right") 
    return fig 

if __name__ == "__main__":
    import regex as re
    
    abidata = abi_to_dict(sys.argv[1])  
    fseq, rseq = generate_consensusseq(abidata)  
    fig = visualize(abidata, template="AGCCGGCTGGCTGCAGGCGT", region="aligned") 
    fig.savefig("test.pdf", bbox_inches="tight") 
    
    abidata = abi_to_dict(sys.argv[1])  
    fseq, rseq = generate_consensusseq(abidata)  
    fig = visualize(abidata, template="AGCCGGCTGGCTGCAGGCGT", region="all") 
    fig.savefig("test2.pdf", bbox_inches="tight") 

    pwm = generate_pwm(abidata) 
    fseq, rseq = generate_consensusseq(abidata) 
    match = re.search("(AGCCGGCTGGCTGCAGGCGT){e<=1}", fseq.upper())
    s,e   = match.span()
    fig = plt.figure(figsize=(0.25,1))
    ax  = fig.add_axes([0.1, 0.1, e-s, 0.75])
    pwm = logomaker.transform_matrix(pwm.iloc[s:e, :], from_type="counts", to_type="probability")
    #pwm = logomaker.transform_matrix(pwm.iloc[s:e, :], from_type="counts", to_type="information")
    logo = logomaker.Logo(pwm,
        font_name='Helvetica',
        color_scheme='classic',
        vpad=.0,
        width=.8,
        ax=ax)
    ax.set_xticks([]) 
    fig.savefig("test_logo.pdf", bbox_inches="tight") 
