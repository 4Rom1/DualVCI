#!/usr/bin/python3.5

'''
 Title of the program:  Convert.py

 Author : Romain Garnier

 Date : 11/12/2017

  There is no copyright or responsibility accepted for the use
  of this piece of code.

'''

import sys

if len(sys.argv) != 3:
 print("Should have 2 arguments")
 print("Usage ./Convert.py input.mop NMode")
 sys.exit()

if sys.argv[1] is None:
  print("Problem reading file")
  sys.exit()

InFile = sys.argv[1]
NMode = 0
NMode = int(sys.argv[2])

File = open(InFile, "r")
if File is not None:
 Txt=File.read()

if(NMode==0) :
 print("Number of modes should be indicated in second argument")
 sys.exit()

File.close()

Lines = Txt.split("\n") 

#NLines=len(Lines)
#print ("Number of lines {}".format(NLines))

import re

Exp1 = r"(^[ -][0-9].[0-9]{22}E[-+][0-9]{2})   Q\^([0-9][0-9]?)\(Q([0-9][0-9]?)\)"
Exp2 = r"(^[ -][0-9].[0-9]{22}E[-+][0-9]{2})   Q\^([0-9][0-9]?)\(Q([0-9][0-9]?)\)   Q\^([0-9][0-9]?)\(Q([0-9][0-9]?)\)"
Exp3 = r"(^[ -][0-9].[0-9]{22}E[-+][0-9]{2})   Q\^([0-9][0-9]?)\(Q([0-9][0-9]?)\)   Q\^([0-9][0-9]?)\(Q([0-9][0-9]?)\)   Q\^([0-9][0-9]?)\(Q([0-9][0-9]?)\)"

TabMode=[]

def InitMode(TabMode,NMode):
 for mm in range(0, NMode):
  TabMode.append(0)

def ZeroMode(TabMode,NMode):
 for mm in range(0, NMode):
  TabMode[mm]=0

def PrintMode(TabMode,NMode,Val):
 Chaine=" "
 for mm in range(0, NMode):
  Chaine+=str(TabMode[mm])+" " 
 Chaine+=", "
 Chaine+=str(Val)
 print(Chaine) 

HA_TO_CM=219474.6435
 
InitMode(TabMode,NMode)

print("FORCEFIELD")

for ii in range(0, len(Lines)):
 
 SubChaine1=re.match(Exp1,Lines[ii])
 SubChaine2=re.match(Exp2,Lines[ii])
 SubChaine3=re.match(Exp3,Lines[ii]) 
 
 ZeroMode(TabMode,NMode)

 if SubChaine1 and not(SubChaine2) is not None:
  NMode1=int(SubChaine1.group(3))
  Pow1=int(SubChaine1.group(2))
  TabMode[NMode1]=Pow1
  Val=float(SubChaine1.group(1))*HA_TO_CM
  if(abs(Val)>1e-20):
   PrintMode(TabMode,NMode,Val)  
#print("Q({})^{}".format(NMode1,Pow1))
  
 elif SubChaine1 and SubChaine2 and not(SubChaine3) is not None:
  NMode1=int(SubChaine1.group(3))
  Pow1=int(SubChaine1.group(2))
  NMode2=int(SubChaine2.group(5))
  Pow2=int(SubChaine2.group(4))
  TabMode[NMode1]=Pow1
  TabMode[NMode2]=Pow2
  Val=float(SubChaine1.group(1))*HA_TO_CM
  if(abs(Val)>1e-20):
   PrintMode(TabMode,NMode,Val) 
#  print("Q({})^{}, Q({})^{}".format(NMode1,Pow1,NMode2,Pow2))
 elif SubChaine1 and SubChaine2 and SubChaine3 is not None:
  NMode1=int(SubChaine1.group(3))
  Pow1=int(SubChaine1.group(2))
  NMode2=int(SubChaine2.group(5))
  Pow2=int(SubChaine2.group(4))
  NMode3=int(SubChaine3.group(7))
  Pow3=int(SubChaine3.group(6))
  TabMode[NMode1]=Pow1
  TabMode[NMode2]=Pow2
  TabMode[NMode3]=Pow3
  Val=float(SubChaine1.group(1))*HA_TO_CM
  if(abs(Val)>1e-20):
   PrintMode(TabMode,NMode,Val) 
#  print("Q({})^{}, Q({})^{}, Q({})^{}".format(NMode1,Pow1,NMode2,Pow2,NMode3,Pow3))
# else:
#  print("No polynomial expression")

print("ENDFF")

