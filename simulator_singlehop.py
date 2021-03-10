import argparse
import numpy as np
import math
from scipy.stats import norm
import random
from collections import deque
import wsnsimpy.wsnsimpy_tk as wsp
import pandas as pd
import matplotlib.pyplot as plt
import time


debugNodeId=8
# this is an array with measured values for sensitivity
# see paper, Table 3
#sf7 = np.array([7,-126.5,-124.25,-120.75])
#sf8 = np.array([8,-127.25,-126.75,-124.0])
#sf9 = np.array([9,-131.25,-128.25,-127.5])
#sf10 = np.array([10,-132.75,-130.25,-128.75])
#sf11 = np.array([11,-134.5,-132.75,-130])
#sf12 = np.array([12,-133.25,-132.25,-132.25])
sf7 = np.array([7,-123,-120,-117.0])
sf8 = np.array([8,-126,-123,-120.0])
sf9 = np.array([9,-129,-126,-123.0])
sf10 = np.array([10,-132,-129,-126.0])
sf11 = np.array([11,-134.53,-131.52,-128.51])
sf12 = np.array([12,-137,-134,-131.0])
#sf array
sensi = np.array([sf7,sf8,sf9,sf10,sf11,sf12])
#ISO threshold values
IS7 = np.array([1,-8,-9,-9,-9,-9])
IS8 = np.array([-11,1,-11,-12,-13,-13])
IS9 = np.array([-15,-13,1,-13,-14,-15])
IS10 = np.array([-19,-18,-17,1,-17,-18])
IS11 = np.array([-22,-22,-21,-20,1,-20])
IS12 = np.array([-25,-25,-25,-24,-23,1])
IsoThresholds = np.array([IS7,IS8,IS9,IS10,IS11,IS12])

Ptx = 9.5
#global values for airtime calculation 
gamma = 2.08
d0 = 40.0
Lpld0 = 127.41
#Bandwiths
BLOW=125
BMED=250
BHIGH=500
#Crc coding rate
CodingRate = 1

#Variance
var=0
# packet size per SFs
PcktLength_SF = [20,20,20,20,20,20]
LorawanHeader = 7
GL=0 #Still not clear
BROADCAST_ADDR = 0xFFFF #Dummy

#Returns the airtime of a packet
def airtime(sf,cr,pl,bw):
    H = 0        # implicit header disabled (H=0) or not (H=1)
    DE = 0       # low data rate optimization enabled (=1) or not (=0)
    Npream = 8   # number of preamble symbol (12.25  from Utz paper)

    if bw == 125 and sf in [11, 12]:
        # low data rate optimization mandated for BW125 with SF11 and SF12
        DE = 1
    if sf == 6:
        # can only have implicit header with SF6
        H = 1

    Tsym = (2.0**sf)/bw  # msec
    Tpream = (Npream + 4.25)*Tsym
    payloadSymbNB = 8 + max(math.ceil((8.0*pl-4.0*sf+28+16-20*H)/(4.0*(sf-2*DE)))*(cr+4),0)
    Tpayload = payloadSymbNB * Tsym
    return ((Tpream + Tpayload)/1000.0)  # to secs

#################################################################
#Computes the interference packet error

def ber_reynders(eb_no, sf):
    """Given the energy per bit to noise ratio (in db), compute the bit error for the SF"""
    return norm.sf(math.log(sf, 12)/math.sqrt(2)*eb_no)

def ber_reynders_snr(snr, sf, bw, cr):
    """Compute the bit error given the SNR (db) and SF"""
    Temp = [4.0/5,4.0/6,4.0/7,4.0/8]
    CR = Temp[cr-1]
    BW = bw*1000.0
    eb_no =  snr - 10*math.log10(BW/2**sf) - 10*math.log10(sf) - 10*math.log10(CR) + 10*math.log10(BW)
    return ber_reynders(eb_no, sf)

def per(sf,bw,cr,rssi,pl):
    snr = rssi  +174 - 10*math.log10(bw) - 6
    return 1 - (1 - ber_reynders_snr(snr, sf, bw, cr))**(pl*8)
######################################################à
#Compute collisions
#time
#frequency
#capture effect + pesudo-orthognal SFs

def timingCollision(p1,p2):
        # assuming p1 is the freshly arrived packet and this is the last check
        # we've already determined that p1 is a weak packet, so the only
        # way we can win is by being late enough (only the first n - 5 preamble symbols overlap)

        # assuming 8 preamble symbols
        Npream = 8

        # we can lose at most (Npream - 5) * Tsym of our preamble
        Tpreamb = 2**p1.sf/(1.0*p1.bw) * (Npream - 5)

        # check whether p2 ends in p1's critical section
        p2_end = p2.arriveTime + p2.rectime
        p1_cs = p1.arriveTime + (Tpreamb/1000.0)  # to sec
        
        if p1_cs < p2_end:
            return True
            
        return False
#
# frequencyCollision, conditions
#
#        |f1-f2| <= 120 kHz if f1 or f2 has bw 500
#        |f1-f2| <= 60 kHz if f1 or f2 has bw 250
#        |f1-f2| <= 30 kHz if f1 or f2 has bw 125
def frequencyCollision(p1,p2):
    if (abs(p1.freq-p2.freq)<=120 and (p1.bw==500 or p2.bw==500)):
        #print "frequency coll 500"
        return True
    elif (abs(p1.freq-p2.freq)<=60 and (p1.bw==250 or p2.bw==250)):
        #print "frequency coll 250"
        return True
    else:
        if (abs(p1.freq-p2.freq)<=30):
            #print "frequency coll 125"
            return True
        #else:
    #print "no frequency coll"
    return False

# check the capture effect and checking the effect of pesudo-orthognal SFs
def powerCollision_2(p1, p2):
    #powerThreshold = 6
    #print "SF: node {0.nodeid} {0.sf} node {1.nodeid} {1.sf}".format(p1, p2)
    #print "pwr: node {0.nodeid} {0.rssi:3.2f} dBm node {1.nodeid} {1.rssi:3.2f} dBm; diff {2:3.2f} dBm".format(p1, p2, round(p1.rssi - p2.rssi,2))
    if p1.sf == p2.sf:
        if abs(p1.rssi - p2.rssi) < IsoThresholds[p1.sf-7][p2.sf-7]:
            #print "collision pwr both node {} and node {}".format(p1.nodeid, p2.nodeid)
            # packets are too close to each other, both collide
            # return both packets as casualties
            return (p1, p2)
        elif p1.rssi - p2.rssi < IsoThresholds[p1.sf-7][p2.sf-7]:
            # p2 overpowered p1, return p1 as casualty
            #print "collision pwr node {} overpowered node {}".format(p2.nodeid, p1.nodeid)
            #print "capture - p2 wins, p1 lost"
            return (p1,)
        #print "capture - p1 wins, p2 lost"
        # p2 was the weaker packet, return it as a casualty
        return (p2,)
    else:
        if p1.rssi-p2.rssi > IsoThresholds[p1.sf-7][p2.sf-7]:
            #print "P1 is OK"
            if p2.rssi-p1.rssi > IsoThresholds[p2.sf-7][p1.sf-7]:
                #print "p2 is OK"
                return ()
            else:
                #print "p2 is lost"
                return (p2,)
        else:
            #print "p1 is lost"
            if p2.rssi-p1.rssi > IsoThresholds[p2.sf-7][p1.sf-7]:
                #print "p2 is OK"
                return (p1,)
            else:
                #print "p2 is lost"
                return (p1,p2)

#################################################################

def delay():
    return random.uniform(.2,.8)


class Packet():
    def __init__(self, sender, dest, freq, sf, bw, cr, txpow, distance, msg, ext, payload):
        global gamma
        global d0
        global var
        global Lpld0
        global GL
        self.payload=payload
        self.ext=ext
        self.msg = msg
        self.sender = sender
        self.dest = dest
        self.freq = freq
        self.sf = sf
        self.bw = bw
        self.cr = cr
        self.txpow = txpow
        self.pl = LorawanHeader+PcktLength_SF[self.sf-7]
        self.symTime = (2.0**self.sf)/self.bw
        self.arriveTime = 0
        if var == 0: 
            Lpl = Lpld0 + 10*gamma*math.log10(distance/d0)
            #print(Lpl)
        else:
            Lpl = Lpld0 + 10*gamma*math.log10(distance/d0) + np.random.normal(-var, var)

        self.rssi = self.txpow - GL - Lpl
        #print(self.rssi)
        self.rectime = airtime(self.sf,self.cr,self.pl,self.bw)
        self.collided=0

class LoraWANNode(wsp.Node):

    def init(self):
        self.cr = CodingRate
        self.txpow = Ptx
        self.set_tx_params()
        self.set_tx_range()
        #self.lastPacket=None
        self.dutyCycle1p = [0,0,0] # 3 channels with 1% duty cycle
        self.dutyCycle10p = 0 # one channel with 10% duty cycle
        self.sendedNo=0 #number of sended packages

    def set_tx_params(self):
        self.bw = BLOW#np.random.choice([BLOW,BMED,BHIGH])
        self.sf = np.random.choice(np.arange(7,12))
        self.freq = np.random.choice([868100000, 868300000, 868500000])

    def send(self,packet,**kwargs):
        #print(self.tx_range)
        for (dist,node) in self.neighbor_distance_list:
            if dist <= self.tx_range:
                sensitivity = sensi[packet.sf - 7, [125,250,500].index(packet.bw) + 1]
                if packet.rssi < sensitivity:
                    #packet not lost
                    #can send after last check
                    if (per(packet.sf,packet.bw,packet.cr,packet.rssi,packet.pl) < random.uniform(0,1)):
                            
                        chanlindex=[868100000, 868300000, 868500000].index(packet.freq)
                        if packet.msg=='ack' or packet.msg=='jacc':
                            acked=0
                            #if duty cycle permits
                            sendtime = self.now + 1  # one sec after receiving the packet
                            if (sendtime >= self.dutyCycle1p[chanlindex]):
                                self.delayed_exec(
                                    (1+packet.rectime),node.on_receive,packet,**kwargs)
                                self.dutyCycle1p[chanlindex]+=packet.rectime*1000
                                acked=1
                            if acked==0:
                                sendtime = self.now + 2
                                if (sendtime >= self.dutyCycle10p):
                                    self.delayed_exec(
                                        2+packet.rectime,node.on_receive,packet,**kwargs)
                                    self.dutyCycle10p+=packet.rectime*100
                        else:
                            if (self.now >= self.dutyCycle1p[chanlindex]):
                                self.delayed_exec(0.00000000001,node.on_receive,packet,**kwargs)
                                    #graphics
                                self.dutyCycle1p[chanlindex]+=packet.rectime*1000  
                                obj_id = self.scene.circle(
                                            self.vispos[0], self.vispos[1],
                                            self.tx_range*scale,
                                            line="wsnsimpy:tx")
                                    #
                                
                                self.delayed_exec(1,self.scene.delshape,obj_id)


    def set_tx_range(self):
        minsensi = np.amin(sensi[:,[125,250,500].index(self.bw) + 1])
        Lpl = Ptx - minsensi
        self.tx_range = d0*(10**((Lpl-Lpld0)/(10.0*gamma)))
        sim.update_neighbor_list(self.id)


class LoraEDNode(LoraWANNode):


    def init(self):
        super().init()
        self.sons = set()
        self.parent = None
        self.parents = {}
        self.EXT = None
        self.receivedP = None
        self.stickToTxParams=0
        self.jreqNo=0
        self.ackedNo=0
        self.lastHash=hash("0")

    def run(self):
        yield self.timeout(1)
        self.scene.nodecolor(self.id,.7,.7,.7)
        while True:
            if self.stickToTxParams==0:
                self.set_tx_params()
                self.set_tx_range()
            if self.dutyCycle1p[([868100000, 868300000, 868500000].index(self.freq))]==0:
                yield self.timeout(delay()*100)
            else:
                yield self.timeout(self.dutyCycle1p[([868100000, 868300000, 868500000].index(self.freq))]*np.random.randint(5,10))
            if self.parent is None:
                self.send_jreq()
            else:
                #print(f"{args.simTime}:{self.now}")
                #if args.simTime-self.now>args.simTime*.5:
                self.send_data()
    
    def send_jreq(self):
        if self.EXT is None:
            self.jreqNo+=1
            packet = Packet(self.id,None,self.freq,self.sf,self.bw,self.cr,self.txpow,self.tx_range,'jreq',None,None)
            self.send(packet)
        
    def send_data(self):
        payload= hash(str(self.id)+str(self.now))
        #self.sendedNo+=1
        self.lastHash=payload
        packet = Packet(self.id,self.parent,self.freq,self.sf,self.bw,self.cr,self.txpow,self.tx_range,'data',self.EXT,payload)
        self.send(packet)

    def on_receive(self, packet, **kwargs):
        if packet.dest==self.id:
            if packet.msg == 'jacc':
                self.parents[packet.sender] = packet.ext
                if self.EXT is None:
                    self.joinTime=self.now
                    self.parent, self.EXT = min(self.parents.items(), key=lambda x: x[1]) 
                    self.EXT+=1
                    self.parent_s = self.scene.addlink(packet.sender,self.id,"parent")
                    self.stickToTxParams=1
            if packet.msg == 'ack':
                if packet.payload == self.lastHash:
                    self.sendedNo+=1
                self.ackedNo+=1
                pass   


class RelayNode(LoraWANNode):


    def init(self):
        super().init()
        self.sons = set()
        self.parents = {}
        self.parent = None
        self.EXT = None
        self.receivedP = deque([])
        self.bandwidth=BLOW
        self.dataQ = deque([])
        self.maxBSReceives=2
        self.collisionWLine=deque([])
        self.stickToTxParams=0
        self.jreqNo=0
        self.ackedNo=0
        self.recId=0
    ###################
    def run(self):
        yield self.timeout(1)
        self.scene.nodecolor(self.id,0,0,1)
        while True:
            if self.stickToTxParams==0:
                self.set_tx_params()
                self.set_tx_range()
            if self.dutyCycle1p[([868100000, 868300000, 868500000].index(self.freq))]==0:
                yield self.timeout(delay()*100)
            else:
                yield self.timeout(self.dutyCycle1p[([868100000, 868300000, 868500000].index(self.freq))])
            if self.parent is None:
                self.send_jreq()
            else:
                self.send_data()
        

    def send_jreq(self):
        if self.EXT is None:
            self.jreqNo+=1
            if self.id == debugNodeId:
                self.log("Sending jreq\n")
            packet = Packet(self.id,None,self.freq,self.sf,self.bw,self.cr,self.txpow,self.tx_range,'jreq',None,None)
            self.send(packet)

    def send_data(self):
        if len(self.dataQ)>0:
            #print(str(self.id)+":"+str(len(self.dataQ)))
            payload=self.dataQ[0]#.popleft()
        #payload= "n:"+str(self.id)+"t:"+str(self.now)
            if self.id == debugNodeId:
                self.log("Sending data to: {}\n".format(self.parent))
            packet = Packet(self.id,self.parent,self.freq,self.sf,self.bw,self.cr,self.txpow,self.tx_range,'data',self.EXT,payload)
            self.send(packet)


    def send_jacc(self,packetIn,**kwargs):
        if self.id == debugNodeId:
            self.log("Sending jacc to: {}\n".format(packetIn.sender))
        packet = Packet(self.id,packetIn.sender,packetIn.freq,packetIn.sf,packetIn.bw,self.cr,packetIn.txpow,self.tx_range,'jacc',self.EXT,None)
        self.send(packet)

    def send_ack(self,packetIn,**kwargs):
        if self.id == debugNodeId:
            self.log("Sending ack to: {}\n".format(packetIn.sender))
        packet = Packet(self.id,packetIn.sender,packetIn.freq,packetIn.sf,packetIn.bw,self.cr,packetIn.txpow,self.tx_range,'ack',self.EXT,packetIn.payload)
        self.send(packet)
    
    def rx_end(self):
        for packet in self.collisionWLine:
            if packet.collided==1:
                x,y = self.vispos
                line1 = self.scene.line(x-5,y-5,x+5,y+5,line="wsnsimpy:collision")
                line2 = self.scene.line(x+5,y-5,x-5,y+5,line="wsnsimpy:collision")
                self.delayed_exec(1,self.scene.delshape,line1)
                self.delayed_exec(1,self.scene.delshape,line2)
            if packet.collided==0:
                if packet.msg == 'jreq':
                    if self.id == debugNodeId:
                        self.log("Received jreq from:{}\n".format(packet.sender))
                    if self.EXT is not None:
                        self.send_jacc(packet)
                        self.sons.add(packet.sender)
                        self.log("Sons:{}".format(self.sons))
                if packet.dest==self.id:
                    if packet.msg == 'data':
                        if self.id == debugNodeId:
                            self.log("Received data from:{}\n".format(packet.sender))
                        self.dataQ.append(packet.payload)
                        self.send_ack(packet)
                    if packet.msg == 'jacc':
                        if self.id == debugNodeId:
                            self.log("Received jacc from: {}\n".format(packet.sender))
                        self.parents[packet.sender] = packet.ext
                        if self.EXT is None:
                            self.joinTime=self.now
                            self.parent, self.EXT = min(self.parents.items(), key=lambda x: x[1]) 
                            self.EXT+=1
                            self.parent_s = self.scene.addlink(packet.sender,self.id,"parent")
                            self.stickToTxParams=1
                    if packet.msg == 'ack':
                        #self.ackedNo+=1
                        #if packet.sender is self.parent:
                        if packet.payload == self.dataQ[0]:
                            self.dataQ.popleft()
                            self.ackedNo+=1
                        if self.id == debugNodeId:
                            self.log("Received ack from:{}\n".format(packet.sender))
        self.collisionWLine=deque([])
        self.recId=0

    def on_receive(self, packet, **kwargs):
        #recid++
        #append packet
        #if recid = 1
        #delay (exec(toa),rx_end)
        #if recid > maxBSReceive 
        #check collision
        #set collision flag

        #def rx_end:
        #for p in cwl
        #if p not collided
        #decode
        #clear cwl
        self.recId+=1
        packet.arriveTime = self.now
        #self.log(len(self.collisionWLine))
        #if len(self.collisionWLine) < self.maxBSReceives:
        self.collisionWLine.append(packet)
        if self.recId==1:
            self.delayed_exec(packet.rectime,self.rx_end)
        if self.recId>self.maxBSReceives:
            #print("COLLISION!!!!")
            for other in self.collisionWLine:
                if other is not packet:
                    if frequencyCollision(packet, other) or timingCollision(packet, other):
                            # Capture + Non-orthognalitiy SFs effects
                        c = powerCollision_2(packet, other) #returns the lost packet
                        for a in c:
                            a.collided=1
                        #collision visualization
                        #x,y = self.vispos
                        #line1 = self.scene.line(x-5,y-5,x+5,y+5,line="wsnsimpy:collision")
                        #line2 = self.scene.line(x+5,y-5,x-5,y+5,line="wsnsimpy:collision")
                        #self.delayed_exec(1,self.scene.delshape,line1)
                        #self.delayed_exec(1,self.scene.delshape,line2)
                        #
"""
    def on_receive(self, packet, **kwargs):
        packet.arriveTime = self.now
        #if len(self.collisionWLine) < self.maxBSReceives:
        self.collisionWLine.append(packet)
        time.sleep(packet.rectime)
        #else:
            #return
        if len(self.collisionWLine)>self.maxBSReceives:
            #print("collision")
            print("COLLISION!!!!")
            for other in self.collisionWLine:
                if other is not packet:
                    if frequencyCollision(packet, other) and timingCollision(packet, other):
                            # Capture + Non-orthognalitiy SFs effects
                        c = powerCollision_2(packet, other) #returns the lost packet
                        for a in c:
                            a.collided=1
                            #collision visualization
                        x,y = self.vispos
                        line1 = self.scene.line(x-5,y-5,x+5,y+5,line="wsnsimpy:collision")
                        line2 = self.scene.line(x+5,y-5,x-5,y+5,line="wsnsimpy:collision")
                        self.delayed_exec(1,self.scene.delshape,line1)
                        self.delayed_exec(1,self.scene.delshape,line2)
                            #

        self.collisionWLine.remove(packet)
        if packet.collided==0:
            if packet.msg == 'jreq':
                if self.id == debugNodeId:
                    self.log("Received jreq from:{}\n".format(packet.sender))
                if self.EXT is not None:
                    self.send_jacc(packet)
                    self.sons.add(packet.sender)
                    self.log("Sons:{}".format(self.sons))
            if packet.dest==self.id:
                if packet.msg == 'data':
                    if self.id == debugNodeId:
                        self.log("Received data from:{}\n".format(packet.sender))
                    self.dataQ.append(packet.payload)
                    self.send_ack(packet)
                if packet.msg == 'jacc':
                    if self.id == debugNodeId:
                        self.log("Received jacc from: {}\n".format(packet.sender))
                    self.parents[packet.sender] = packet.ext
                    if self.EXT is None:
                        self.joinTime=self.now
                        self.parent, self.EXT = min(self.parents.items(), key=lambda x: x[1]) 
                        self.EXT+=1
                        self.parent_s = self.scene.addlink(packet.sender,self.id,"parent")
                        self.stickToTxParams=1
                if packet.msg == 'ack':
                    #self.ackedNo+=1
                    #if packet.sender is self.parent:
                    if packet.payload == self.dataQ[0]:
                        self.dataQ.popleft()
                        self.ackedNo+=1
                    if self.id == debugNodeId:
                        self.log("Received ack from:{}\n".format(packet.sender))
"""

class LoraSink(LoraWANNode):

    ###################
    def init(self):
        super().init()
        self.sons = set()
        self.parent = None
        self.EXT = 0
        self.dataQ = deque([])
        self.receivedP = deque([])
        self.maxBSReceives=8
        self.collisionWLine=deque([])
        self.recId=0
        self.collNo=0

    ###################
    def run(self):
        yield self.timeout(1)
        self.scene.nodecolor(self.id,1,0,0)
        self.scene.nodewidth(self.id,2)

    ###################
    def send_jacc(self,packetIn,**kwargs):
        packet = Packet(self.id,packetIn.sender,packetIn.freq,packetIn.sf,packetIn.bw,self.cr,self.txpow,self.tx_range,'jacc',self.EXT,None)
        self.send(packet)

    def send_ack(self,packetIn,**kwargs):
        packet = Packet(self.id,packetIn.sender,packetIn.freq,packetIn.sf,packetIn.bw,self.cr,packetIn.txpow,self.tx_range,'ack',self.EXT,packetIn.payload)        
        self.send(packet)

    ###################
    def rx_end(self):
        for p in self.collisionWLine:
            if p.collided==1:
                self.collNo+=1
                x,y = self.vispos
                line1 = self.scene.line(x-5,y-5,x+5,y+5,line="wsnsimpy:collision")
                line2 = self.scene.line(x+5,y-5,x-5,y+5,line="wsnsimpy:collision")
                self.delayed_exec(1,self.scene.delshape,line1)
                self.delayed_exec(1,self.scene.delshape,line2)
            if p.collided==0:
                if p.msg == 'jreq':
                    if self.EXT is not None:
                        self.send_jacc(p)
                        self.sons.add(p.sender)
                if p.msg == 'data':
                    if p.dest==self.id:
                        self.dataQ.append(p.payload)
                        self.send_ack(p)
        self.collisionWLine=deque([])
        self.recId=0


    def on_receive(self, packet, **kwargs):
        #recid++
        #append packet
        #if recid = 1
        #delay (exec(toa),rx_end)
        #if recid > maxBSReceive 
        #check collision
        #set collision flag

        #def rx_end:
        #for p in cwl
        #if p not collided
        #decode
        #clear cwl
        self.recId+=1
        packet.arriveTime = self.now
        #self.log(len(self.collisionWLine))
        #if len(self.collisionWLine) < self.maxBSReceives:
        self.collisionWLine.append(packet)
        if self.recId==1:
            self.delayed_exec(packet.rectime,self.rx_end)
        if self.recId>self.maxBSReceives:
            #print("COLLISION!!!!")
            for other in self.collisionWLine:
                if other is not packet:
                    if frequencyCollision(packet, other) or timingCollision(packet, other):
                            # Capture + Non-orthognalitiy SFs effects
                        c = powerCollision_2(packet, other) #returns the lost packet
                        for a in c:
                            a.collided=1
                        #collision visualization
                        #x,y = self.vispos
                        #line1 = self.scene.line(x-5,y-5,x+5,y+5,line="wsnsimpy:collision")
                        #line2 = self.scene.line(x+5,y-5,x-5,y+5,line="wsnsimpy:collision")
                        #self.delayed_exec(1,self.scene.delshape,line1)
                        #self.delayed_exec(1,self.scene.delshape,line2)
                        #

       

#not sure if can remove
debugNodeId=3
scale=100
tsx,tsy=(500,500)
tsx*=2
tsy*=2
terrSize=(tsx,tsy)
sim = wsp.Simulator(
        until=100000,
        timescale=0,
        visual=True,
        terrain_size=terrSize,
        title="LoraTree")
# place nodes
sink=LoraSink
sim.add_node(sink,(250,250),1)
sim.add_node(LoraEDNode,(300,200),1)
sim.add_node(LoraEDNode,(200,300),1)
sim.add_node(LoraEDNode,(200,200),1)
sim.add_node(LoraEDNode,(400,400),1)
sim.add_node(LoraEDNode,(250,200),1)
sim.add_node(LoraEDNode,(250,400),1)
sim.add_node(LoraEDNode,(300,250),1)
sim.add_node(LoraEDNode,(100,220),1)
sim.add_node(LoraEDNode,(400,220),1)

sim.scene.linestyle("parent", color=(0,.8,0), arrow="tail", width=2)
sim.run()
sen=0
rec=0
coll=0
for n in sim.nodes:
    if n.id == 0:
        rec=len(n.dataQ)
        coll=n.collNo
    else:
        sen+=n.sendedNo
print(f"{sen}:{rec}:{coll}")
"""
relNo=args.relay
ran=list(range(xGrid*yGrid-1))
sinksIds=np.random.choice(range(xGrid*yGrid-1), sinksNo, replace=False)
rel=[nid for nid in ran if nid not in sinksIds]
relayIds=np.random.choice(rel, relNo, replace=False)
nodeID=0
for x in range(xGrid): 
    for y in range(yGrid):
        px =  x*(terrSize[0]/xGrid) + random.uniform(-(terrSize[0]/(xGrid*6)),(terrSize[0]/(xGrid*6)))
        py =  y*(terrSize[1]/yGrid) + random.uniform(-(terrSize[1]/(xGrid*6)),(terrSize[1]/(xGrid*6)))
        if nodeID in relayIds:
            newNode = RelayNode
        else:
            newNode = LoraEDNode
        ##loraNodes+=1
        if nodeID in sinksIds:
            newNode = LoraSink
            #px = 50 + xGrid*100/2
            #py = 50 + yGrid*100/2    
        #if newNode.id is SINK:
        #    newNode = MyNode
        node = sim.add_node(newNode, (px,py), scale)
        #node.tx_range = 40
        #node.logging = False
        #if id==debugNodeId:
        node.logging = False
        nodeID+=1


# define a line style for parent links
sim.scene.linestyle("parent", color=(0,.8,0), arrow="tail", width=2)

# start the simulation
sim.run()

datas = map(lambda gw: (gw.id,gw.dataQ),filter(lambda node: node.id in sinksIds,sim.nodes))
gws=list(filter(lambda node: node.id in sinksIds,sim.nodes))
rels=list(filter(lambda node: node.id in relayIds,sim.nodes))
eds=list(filter(lambda node: node.id not in relayIds and node.id not in sinksIds,sim.nodes))
received=sum([len(gw.dataQ) for gw in gws])
#print(list(map(lambda gw: (gw.id,len(gw.dataQ)),gws)))
sended=sum([ed.sendedNo for ed in eds])
sonsMean=sum([len(rel.sons) for rel in rels])/len(rels)+sum([len(gw.sons) for gw in gws])/len(gws)
notConnected=0
for rel in rels:
    if rel.parent is None:
        notConnected+=1
for ed in eds:
    if ed.parent is None:
        notConnected+=1
print(xGrid*yGrid,":",args.gw,":",relNo,":",sended,":",received,":",sonsMean,":",notConnected,end="")

#print(list(datas))"""