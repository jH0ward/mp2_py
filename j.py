#!/usr/bin/python

alist=['a','b','c','d','e','f']

orb_pairs=[(i,j) for i in range(len(alist)) for j in range(len(alist)) if i<=j]
print orb_pairs
count=0
for pindex1 in range(len(orb_pairs)):
    for pindex2 in range(0,pindex1+1):
        print "Left side of intrgral = ", orb_pairs[pindex1]
        print "Right side of intrgral = ", orb_pairs[pindex2]
        print "----"
        count+=1

print "total = ",count

