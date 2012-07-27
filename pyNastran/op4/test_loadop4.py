#!/usr/bin/env python

import numpy as np
from cop4 import OP4
from pyNastran.general.generalMath import printMatrix

for in_file in [ #'test/mat_b_dn.op4' ,
                 #'test/mat_b_s1.op4' ,
                 #'test/mat_b_s2.op4' ,
                 #'test/mat_t_dn.op4' ,
                 #'test/mat_t_s1.op4' ,
                 #'test/mat_t_s2.op4' , 
                 'test/b_sample.op4',
                 ]:
    try:
        op4fh = OP4(in_file, 'r')
    except:
        print('Failed to get header of %s, ignoring.' % (in_file))
        continue

    print('%s\n%s' % ('=' * 61, in_file))
    op4fh.print_header()

    for i in range(op4fh.nmat):
        a = op4fh.Load(nmat=1, skip=i)
        if a is None: 
            print('Failed to get %d-th matrix' % i)
        else:
            print('%s:' % op4fh.name[i])
            print printMatrix(a)
