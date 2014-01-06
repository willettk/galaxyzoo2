import numpy as np
import os
import Image,ImageMath


testdir = ('/Users/willettk/Astronomy/Research/GalaxyZoo/kaggle')
original = Image.open(os.path.join(testdir,'benchtest.jpg'))
rotated = original.rotate(180)

xsize,ysize = original.size
xsize_r,ysize_r = rotated.size

assert xsize == xsize_r, \
    'X sizes do not match, %i, %i' % (xsize,xsize_r)

assert ysize == ysize_r, \
    'Y sizes do not match, %i, %i' % (ysize,ysize_r)

r,g,b = original.split()
rr,gr,br = rotated.split()

imdiff_r = ImageMath.eval("convert(abs(a-b),'L')",a=r,b=rr)
imdiff_g = ImageMath.eval("convert(abs(a-b),'L')",a=g,b=gr)
imdiff_b = ImageMath.eval("convert(abs(a-b),'L')",a=b,b=br)

#imdiff_r = Image.fromarray(diff_r)
#imdiff_g = Image.fromarray(diff_g)
#imdiff_b = Image.fromarray(diff_b)

imnew = Image.merge('RGB', (imdiff_r,imdiff_g,imdiff_b))

imnew.save(os.path.join(testdir,'benchtest_asymmetry.jpg'))
rotated.save(os.path.join(testdir,'benchtest_rotated.jpg'))

imdiff_r.save(os.path.join(testdir,'benchtest_ra.jpg'))
imdiff_g.save(os.path.join(testdir,'benchtest_ga.jpg'))
imdiff_b.save(os.path.join(testdir,'benchtest_ba.jpg'))
r.save(os.path.join(testdir,'benchtest_r.jpg'))
g.save(os.path.join(testdir,'benchtest_g.jpg'))
b.save(os.path.join(testdir,'benchtest_b.jpg'))
rr.save(os.path.join(testdir,'benchtest_rr.jpg'))
gr.save(os.path.join(testdir,'benchtest_gr.jpg'))
br.save(os.path.join(testdir,'benchtest_br.jpg'))

# Isolate the inner 100x100 pixels and use that as a hash?
