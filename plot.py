import matplotlib.pyplot as plt
import pandas as pd
import os
import csv
#plotta solo i risultachi che hanno #nc*100/#tot<30
fn=[]
pdr=[]
nc=[]
sm=[]
tot=[]
gw=[]
rel=[]
noc=[]
toShow=[]
directory = os.path.join("./","testc")
for root,dirs,files in os.walk(directory):
    for file in files:
       if file.endswith(".csv"):
            f=open(directory+"/"+file, 'r')
            df_in=pd.read_csv(f)
            df=pd.DataFrame(df_in[df_in['NC']/(df_in['T'][0]-df_in['G'][0])<0.3]).reset_index(drop=True)

            fn.append(f.name)
            #df=pd.read_csv(f)
            #print(df)
            tot.append(df['T'][0])
            gw.append(df['G'][0])
            rel.append(df['Rel'][0])
            noc.append(df['MNoC'].mean()/(df['Rel'][0]+df['G'][0]))
            df['Pdr']=df['R']/df['S']*100
            pdr.append(df['Pdr'].mean())
            sm.append(df['MNoS'].mean()*100/(df['Rel'][0]+df['G'][0]))
            nc.append(df['NC'].mean()*100/(df['T'][0]-df['G'][0]))
            
            f.close()

for i in range(len(fn)):
    toShow.append((tot[i],gw[i],rel[i],pdr[i],noc[i],nc[i],sm[i]))
toShow.sort(key=lambda tup:tup[0])
xNames=[f"{i[0]},{i[1]},{i[2]}" for i in toShow]
#for i in toShow:
#    xNames.append(f"{i[0]},{i[1]},{i[2]}")
x=xNames
#print(toShow)

y1=[i[3] for i in toShow]
y2=[i[4] for i in toShow]
y3=[i[6] for i in toShow]
y4=[i[5] for i in toShow]
for c in toShow:
    print(f"Tot:{c[0]},Gws:{c[1]},Rep:{c[2]},PdrM:{c[3]},CollM:{c[4]},NCM:{c[5]},SNM:{c[6]}")
fig, axs = plt.subplots(4,sharex=True)
axs[0].set_title("Packet delivery rateo")
axs[0].plot(x,y1)
axs[1].set_title("Mean collisions # per node")
axs[1].plot(x,y2)
axs[2].set_title("Mean sons # per node")
axs[2].plot(x,y3)
axs[3].set_title("% Not Connected")
axs[3].plot(x,y4)
fig.tight_layout(pad=2.0)
fig.canvas.set_window_title('Performance')
#plt.legend(loc='upper left')
plt.show()