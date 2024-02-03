def legendFunction(runsdict):
    conversion={"d0":"cavityd0","d-200":"cavityd-200","slope200":"cavityd200","slope375":"cavityd375"}
    for k in runsdict.keys():
        for l in range(len(runsdict[k]["specialstring"])):
            if not runsdict[k]["specialstring"][0]:
                ss = k.split("-")[0]
            else:
                ss = runsdict[k]["specialstring"][l]
            if ss in conversion.keys():
                ss = conversion[ss]
            plt.scatter(1,1,marker=runsdict[k]["marker"][l],color=runsdict[k]["color"][l],label=ss)
    plt.legend()
    plt.show()
            
def folderMapTimeSeries(runsdict,save=True):
    fig,axises = plt.subplots(2,3,figsize=(8,7))
    for k in runsdict.keys():
        for f in glob.glob(str("/home/garrett/Projects/MITgcm_ISC/experiments/"+k+"/*"), recursive = True):
            for l in range(len(runsdict[k]["specialstring"])):
                key=runsdict[k]["specialstring"][l]
                if key and key in f:
                    try:
                        timeSeriesDashboard(f+"/results",key+f[-6:-10],fig,axises,color=runsdict[k]["color"][l])
                    except:
                        print("yeesh")
                else:
                    try:
                        timeSeriesDashboard(f+"/results","",fig,axises,color=runsdict[k]["color"][l])
                    except:
                        print("yeesh")
    axises[0][0].legend()

    plt.show()

def fastExplosionCheck(runsdict,save=True):
    fig,axises = plt.subplots(2,3,figsize=(8,7))
    count = [0,0]
    for k in tqdm(runsdict.keys()):
        for f in glob.glob(str("/home/garrett/Projects/MITgcm_ISC/experiments/"+k+"/*"), recursive = True):
            for l in range(len(runsdict[k]["specialstring"])):
                key=runsdict[k]["specialstring"][l]
                if key and key in f:
                    
                    try:
                        with open(f+"/results/STDOUT.0000") as myfile:
                            if 'NaN' in myfile.read():
                                print(f, "oh no a NaN!")
                            else:
                                count[0]=count[0]+1
                            count[1]=count[1]+1
                    except:
                        print("no STDOUT")
    print(count)

def generateRunsTable(runsdict):
    table = []
    prettynames = {'shelf_depth':"Nominal depth of shelf", \
            'rng_seed':"Random bathymetry seed used", 'random_amplitude':"Amplitude of random bathymetry",\
            'cavity_depth':"Depth of cavity relative to depth of shelf", 'cavity_width': "Width of cavity",\
            'yicefront':"Northward extent of ice shelf"}#, 'tcline_atshelf_depth': "Depth of temperature maximum"}
    for k in runsdict.keys():
        for f in glob.glob(str("/home/garrett/Projects/MITgcm_ISC/experiments/"+k+"/*"), recursive = True):
            for l in range(len(runsdict[k]["specialstring"])):
                key=runsdict[k]["specialstring"][l]
                fname = f+"/" #+"results/"
                variables = grabMatVars(fname,("experiment_parameters"))
                fields = variables["experiment_parameters"][0][0].__dict__
                conversion={"d0":"cavityd0","d-200":"cavityd-200","slope200":"cavityd200","slope375":"cavityd375"}
                d={}
                if not runsdict[k]["specialstring"][0]:
                    ss = k.split("-")[0]
                else:
                    ss = runsdict[k]["specialstring"][l]
                if ss in conversion.keys():
                    ss = conversion[ss]
                d["Experiment Name"] = ss
                for j in fields.keys():
                    if j in prettynames.keys():
                        d[prettynames[j]] = fields[j][0][0]
                if "Northward extent of ice shelf" not in d:
                    d["Northward extent of ice shelf"]=150
                if d not in table:
                    table.append(d)
    for k in table:print(k.keys())
    print(tabulate(table,headers="keys",tablefmt="latex",maxcolwidths=[9]*len(table[0].keys())))
    return table


