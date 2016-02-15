#! /usr/bin/env python
#
#  pgex1: freely taken after PGDEMO1.F
#
import ppgplot, Numeric
import sys

# create an array 
xs=Numeric.array([1.,2.,3.,4.,5.])
ys=Numeric.array([1.,4.,9.,16.,25.])

# creat another array
yr = 0.1*Numeric.array(range(0,60))
xr = yr*yr


# pgplotting
if len(sys.argv) > 1: # if we got an argument use the argument as devicename
	ppgplot.pgopen(sys.argv[1])
else:
	ppgplot.pgopen('/xwin')
ppgplot.pgenv(0.,10.,0.,20.,0,1)
ppgplot.pglab('(x)', '(y)', 'PGPLOT Example 1:  y = x\u2')
ppgplot.pgpt(xs,ys,9)
ppgplot.pgline(xr,yr)
ppgplot.pgclos()
