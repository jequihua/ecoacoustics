#!/usr/bin/python
import sys
string=sys.argv[1]
out=[]
string=string.replace('(','').replace(')','')
splitstr=string.split(',')
for item in splitstr:
 if "x" in item:
  m=item.split('x')
  num=int(m[0])
  count=int(m[1])
  for z in range(0,count):
   out.append(str(num))
 else:
  out.append(item)
newstr=",".join(out)
print newstr
