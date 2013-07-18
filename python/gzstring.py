import json
import operator
import re

gzdir = '/Users/willettk/Astronomy/Research/GalaxyZoo/'

def load_data():

    filename = 'subjects_for_chris.csv'

    f = open(gzdir+filename,'rb')
    a = f.readlines()
    f.close()

    return a

def parse_type(glist):

    sloan = []
    candels = []
    for g in glist:
        prefix = json.loads(g.split('\t')[3].split('\n')[0]).keys()[0]
        if prefix[:7] == 'candels':
            candels.append(g)
        if prefix[:5] == 'sloan':
            sloan.append(g)

    return sloan, candels

def max_item(jdict):

    mi = max(jdict.iteritems(), key=operator.itemgetter(1))[0]

    return mi

def write_results():

    filename = 'subjects_for_chris_wstrings.csv'   

    glist = load_data()
    glist_new = []

    f = open(gzdir+filename,'wb')

    for g in glist:
        prefix = json.loads(g.split('\t')[3].split('\n')[0]).keys()[0]
        splititem = g.split('\t')
        if prefix[:7] == 'candels':
            galstring = candels_tree(g)
        if prefix[:5] == 'sloan':
            galstring = sloan_tree(g)
        splititem[1] = galstring
        glist_new.append(galstring)

        writeitem = '\t'.join(splititem)

        f.write(writeitem)

    f.close()

    return glist_new

def sloan_tree(gal):

    j = json.loads(gal.split('\t')[3].split('\n')[0])
    keys = j.keys()
    objname = gal.split('\t')[0]
    assert 'sloan-0' in keys, \
        'Cannot find sloan-0 in keys for %s' % objname

    char = ''

    # First answer

    s0_max = max_item(j['sloan-0'])
    char += 's0%s;' % re.sub('-','',str(s0_max))

    # Star/artifact

    if s0_max != 'a-2':

        # Smooth galaxies

        if s0_max == 'a-0':

            s7_max = max_item(j['sloan-7'])
            char += 's7%s;' % re.sub('-','',str(s7_max))

        # Features/disk

        if s0_max == 'a-1':
            s1_max = max_item(j['sloan-1'])
            char += 's1%s;' % re.sub('-','',str(s1_max))
            # Edge-on disk
            if s1_max == 'a-0':
                # Bulge shape
                s8_max = max_item(j['sloan-8'])
                char += 's8%s;' % re.sub('-','',str(s8_max))
            # Not edge-on disk
            else:
                # Bulge prominence
                s4_max = max_item(j['sloan-4'])
                char += 's4%s;' % re.sub('-','',str(s4_max))
                # Bar
                s2_max = max_item(j['sloan-2'])
                char += 's2%s;' % re.sub('-','',str(s2_max))
                # Spiral
                s3_max = max_item(j['sloan-3'])
                char += 's3%s;' % re.sub('-','',str(s3_max))
                if s3_max == 'a-0':
                    # Arm tightness
                    s9_max = max_item(j['sloan-9'])
                    char += 's9%s;' % re.sub('-','',str(s9_max))
                    # Arm multiplicity
                    s10_max = max_item(j['sloan-10'])
                    char += 's10%s;' % re.sub('-','',str(s10_max))

        # Anything odd?

        s5_max = max_item(j['sloan-5'])
        char += 's5%s;' % re.sub('-','',str(s5_max))
        if s5_max == 'a-0':
            s6_max = max_item(j['sloan-6'])
            char += 's6%s;' % re.sub('[", -]','',str(s6_max[1:-1])) 

    return char

def candels_tree(gal):

    j = json.loads(gal.split('\t')[3].split('\n')[0])
    keys = j.keys()
    objname = gal.split('\t')[0]
    assert 'candels-0' in keys, \
        'Cannot find candels-0 in keys for %s' % objname
 
    char = ''
 
    # First answer
 
    c0_max = max_item(j['candels-0'])
    char += 'c0%s;' % re.sub('-','',str(c0_max))
 
    # Star/artifact
 
    if c0_max != 'a-2':
 
        # Smooth galaxies
 
        if c0_max == 'a-0':
 
            c1_max = max_item(j['candels-1'])
            char += 'c1%s;' % re.sub('-','',str(c1_max))
 
        # Features/disk
 
        if c0_max == 'a-1':
            c2_max = max_item(j['candels-2'])
            char += 'c2%s;' % re.sub('-','',str(c2_max))
            # Clumpy galaxies
            if c2_max == 'a-0':
 
                # How many clumps?
                c3_max = max_item(j['candels-3'])
                char += 'c3%s;' % re.sub('-','',str(c3_max))
 
                # 1 clump
                if c3_max == 'a-0':
                    # Does galaxy appear symmetrical?
                    c7_max = max_item(j['candels-7'])
                    char += 'c7%s;' % re.sub('-','',str(c7_max))
                    # Are clumps embedded?
                    c8_max = max_item(j['candels-8'])
                    char += 'c8%s;' % re.sub('-','',str(c8_max))
 
                else:
                    # 3+ clumps
                    if c3_max != 'a-1':
                        c4_max = max_item(j['candels-4'])
                        char += 'c4%s;' % re.sub('-','',str(c4_max))
 
                    # Is one clump clearly bigger than others?
                    c5_max = max_item(j['candels-5'])
                    char += 'c5%s;' % re.sub('-','',str(c5_max))
 
                    if c5_max == 'a-0':
                        # Is brightest clump central?
                        c6_max = max_item(j['candels-6'])
                        char += 'c6%s;' % re.sub('-','',str(c6_max))
 
                        if c6_max == 'a-0':
                            # Is galaxy symmetrical?
                            c7_max = max_item(j['candels-7'])
                            char += 'c7%s;' % re.sub('-','',str(c7_max))
                            # Are clumps embedded?
                            c8_max = max_item(j['candels-8'])
                            char += 'c8%s;' % re.sub('-','',str(c8_max))
 
                    else:
                        # Is galaxy symmetrical?
                        c7_max = max_item(j['candels-7'])
                        char += 'c7%s;' % re.sub('-','',str(c7_max))
                        # Are clumps embedded?
                        c8_max = max_item(j['candels-8'])
                        char += 'c8%s;' % re.sub('-','',str(c8_max))
 
 
            else:
 
                # Not clumpy galaxies
 
                c9_max = max_item(j['candels-9'])
                char += 'c9%s;' % re.sub('-','',str(c9_max))
                # Edge-on disks
                if c9_max == 'a-0':
                    # Edge-on bulge?
                    c10_max = max_item(j['candels-10'])
                    char += 'c10%s;' % re.sub('-','',str(c10_max))
                # Not edge-on disks
                else:
                    # Bulge prominence
                    c15_max = max_item(j['candels-15'])
                    char += 'c15%s;' % re.sub('-','',str(c15_max))
                    # Bars
                    c11_max = max_item(j['candels-11'])
                    char += 'c11%s;' % re.sub('-','',str(c11_max))
                    # Spiral structure
                    c12_max = max_item(j['candels-12'])
                    char += 'c12%s;' % re.sub('-','',str(c12_max))
                    if c12_max == 'a-0':
                        # Arm tightness
                        c13_max = max_item(j['candels-13'])
                        char += 'c13%s;' % re.sub('-','',str(c13_max))
                        # Arm multiplicity
                        c14_max = max_item(j['candels-14'])
                        char += 'c14%s;' % re.sub('-','',str(c14_max))
 
        # Merging/tidal debris
 
        c16 = j['candels-16']
        char += 'c16%s;' % re.sub('-','',str(c16.keys()[0]))

    return char

