#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 09:49:03 2021

@author: nav
"""

from qiskit import*
from qiskit.tools.visualization import*
from numpy import*
from matplotlib.pyplot import*
qr=QuantumRegister(2)
cr=ClassicalRegister(2)
qc=QuantumCircuit(qr,cr)
m12=7.39*10**-5
m13=2.525*10**-3
m23=m13
th12=0.5902
th23=0.8656
th13=0.1502
LE=linspace(0,1000,11)#(L/E)
le=linspace(0,1000,101)
muon_counts=[]
electron_counts=[]
tau_counts=[]
def phi_1(E):
    return(2*(1.27*m12*E))
def phi_2(E):
    return(2*(1.27*m13*E))
def prob_t(E):
    return(cos(th13)**2*(sin(2*th23))**2*sin(1.27*m23*E))
def prob_e(E):
    return(sin(2*th13)**2*(sin(th23))**2*sin(1.27*m23*E))  
def prob_m(E):
    return(1-cos(th13)**2*(sin(2*th23))**2*sin(1.27*m23*E)-sin(2*th13)**2*(sin(th23))**2*sin(1.27*m23*E))

#circuit
backend=BasicAer.get_backend('qasm_simulator')
for E in LE:
    p=phi_1(E)
    q=phi_2(E)
    qc.x(qr[1])
    qc.u(-0.7053,0,0,qr[0])
    qc.u(-1.3599,0,0,qr[1])
    qc.cx(qr[1],qr[0])
    qc.u(0.7966,0,0,qr[0])
    qc.u(-1.0139,0,0,qr[1])
    qc.cx(qr[1],qr[0])
    qc.u(0.6031,0,0,qr[0])
    qc.u(2.0125,0,0,qr[1])
    qc.u(0,0,p,qr[0])
    qc.u(0,0,q,qr[1])
    qc.u(-0.6031,0,0,qr[0])
    qc.u(-2.0125,0,0,qr[1])
    qc.cx(qr[1],qr[0])
    qc.u(-0.7966,0,0,qr[0])
    qc.u(1.0139,0,0,qr[1])
    qc.cx(qr[1],qr[0])
    qc.u(0.7053,0,0,qr[0])
    qc.u(1.3599,0,0,qr[1])

    qc.measure(qr[0],cr[0])
    qc.measure(qr[1],cr[1])
    counts=execute(qc,backend,shots=1024).result().get_counts()
    electron_counts.append(counts.get('00',0)/1024)
    muon_counts.append(counts.get('01',0)/1024)
    tau_counts.append(counts.get('10',0)/1024)
    qc.reset(qr[0])
    qc.reset(qr[1])
       
AC_t=prob_t(le)
AC_e=prob_e(le)
AC_m=prob_m(le)

plot(le,AC_t,label='tau')
plot(le,AC_e,label='electron')
plot(le,AC_m,label='muon')
plot(LE,electron_counts,'o',label='elec-Qasm')
plot(LE,muon_counts,'o',label='muon-Qasm')
plot(LE,tau_counts,'o',label='tau-Qasm')
legend()