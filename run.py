import os
import argparse
import multiprocessing
import pandas as pd
import io
import progressbar
import time

flatten = lambda l: [item for sublist in l for item in sublist]
par = argparse.ArgumentParser()
par.add_argument("num",type=int)
par.add_argument("terrSize",type=int,nargs=2,metavar=("x","y"))
par.add_argument("gridSize",type=int,nargs=2,metavar=("x","y"))
par.add_argument("gw",type=int)
par.add_argument("relay",type=int)
par.add_argument("simTime",type=int)
par.add_argument("-s","--scale",type=int,help="Percentual value of scale, from 1% to 100%, default 15%")
par.add_argument("--visual",action="store_true")

args = par.parse_args()
xTerr,yTerr=args.terrSize
xGrid,yGrid=args.gridSize

def execute(i):
    print(f"Started simulation #{i}")
    if args.scale is None:
        if args.visual:
            #expRes.append(subprocess.run(["python" , f"simulator.py {xTerr} {yTerr} {xGrid} {yGrid} {args.gw} {args.relay} {args.simTime} --visual"], stdout=subprocess.PIPE).stdout.decode('utf-8'))
            return os.popen(f"python simulator.py {xTerr} {yTerr} {xGrid} {yGrid} {args.gw} {args.relay} {args.simTime} --visual").read()
        else:
            #expRes.append(subprocess.run(["python" , f"simulator.py {xTerr} {yTerr} {xGrid} {yGrid} {args.gw} {args.relay} {args.simTime}"], stdout=subprocess.PIPE).stdout.decode('utf-8'))
            return os.popen(f"python simulator.py {xTerr} {yTerr} {xGrid} {yGrid} {args.gw} {args.relay} {args.simTime}").read()
    else:
        if args.visual:
            #expRes.append(subprocess.run(["python" , f"simulator.py {xTerr} {yTerr} {xGrid} {yGrid} {args.gw} {args.relay} {args.simTime} -s {args.scale} --visual"], stdout=subprocess.PIPE).stdout.decode('utf-8'))
            return os.popen(f"python simulator.py {xTerr} {yTerr} {xGrid} {yGrid} {args.gw} {args.relay} {args.simTime} -s {args.scale} --visual").read()
        else:
            #expRes.append(subprocess.run(["python",f"simulator.py {xTerr} {yTerr} {xGrid} {yGrid} {args.gw} {args.relay} {args.simTime} -s {args.scale}"], stdout=subprocess.PIPE).stdout.decode('utf-8'))
            return os.popen(f"python simulator.py {xTerr} {yTerr} {xGrid} {yGrid} {args.gw} {args.relay} {args.simTime} -s {args.scale}").read()

res=[]
process_pool = multiprocessing.Pool()

print(f"Starting {args.num} simulations")
bar = progressbar.ProgressBar(max_value=args.num)
res=process_pool.map_async(execute,list(range(args.num)))
process_pool.close()
process_pool.join()
expRes=res.get()

expRes.insert(0,"TS:T:G:Rel:S:R:MNoC:MNoS:NC")
df = pd.read_csv(io.StringIO('\n'.join(expRes)), sep=":")
if not os.path.exists('tests'):
    os.makedirs('tests')
filename=r"tests/xs"+f"_{xTerr}_ys_{yTerr}_xg_{xGrid}_yg_{yGrid}_gw_{args.gw}_rel_{args.relay}_t_{args.simTime}.csv"
with open(filename, 'a') as f:
    df.to_csv(f, mode='a', header=f.tell()==0,index=False)
#df.to_csv(r"test/xs"+f"_{xTerr}_ys_{yTerr}_xg_{xGrid}_yg_{yGrid}_gw_{args.gw}_rel_{args.relay}_t_{args.simTime}.csv", index = False)
#print(df)
#df = pd.DataFrame(expRes[0], sep=":")
print("Done.")