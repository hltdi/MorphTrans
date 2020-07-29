"""
This file is part of the MorphTrans project: https://github.com/hltdi/MorphTrans

    Copyleft 2020. MorphTrans Collaborative.

    MorphTrans is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MorphTrans is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MorphTrans.  If not, see <http://www.gnu.org/licenses/>.

Author: Michael Gasser <gasser@indiana.edu>

-- 29.1.2020
   Created.
"""

import mt

A = mt.Aligner('am', 'ti', datafiles=["A04.tr", "A24.tr", "A34.tr"])
B = mt.Aligner('amG', 'tiG', datafiles=["A04ግ.tr", "A24ግ.tr", "A34ግ.tr"])
C = mt.Aligner('am', 'ti', datafiles=["A04G.tr", "A24G.tr", "A34G.tr"])
D = mt.Aligner('am_', 'ti_', datafiles=["A04_.tr", "A24_.tr", "A34_.tr"])
a = A.make_alignment(data_i=1, explicit=[0, 1, 3, -1, 6, 7, 8, -1, 9],
                     dscat=mt.Dataset.TEST)

def runD():
    D.run(restore="A_20200728204220.npy")

def test(aligner):
    for di in [1, 3, 25, 40, 59, 60, 75, 80, 100]:
        print(aligner.test.alignments[di])
        print()
#a.constraints = [0, 1, ['a', 2, 3], 'a', ['b', 4, 6], 7, ['c', 8, 8], 'c', 9]
#c1 = [((0, 0), (1, 1), (5, 7), (8, 9)),
#      [( (2, 3), (3, 3) ), ( (4, 4), (4, 6) ), ( (6, 8), (7, 8) )]
#      ]
# a
# ' a l : a     k a t m
# ' a y l e ' a K a   n
# b
# ' a   l :     a k a t m
# ' a y l   e ' a K a   n
b = A.make_alignment(data_i=1, explicit=[0, -1, 1, 2, 3, -1, 6, 8, 9],
                     dscat=mt.Dataset.TEST, direction=2)
#b = A.make_alignment(data_i=1, explicit=[0, 1, 3, -1, -1, 7, 8, -1, 9],
#                     dscat=mt.Aligner.TEST)
#c = A.make_alignment(data_i=1, explicit=[0, 1, 3, -1, 6, 7, -1, 8, 9],
#                     dscat=mt.Aligner.TEST)
#a.constraints = c1
#b.constraints = c1
#c.constraints = c1
#b = A.make_alignment(0)
#k = A.make_alignment(17)
##' a l g e n : e n k u m
##' a y g e n   e n k u n
#d = A.make_alignment(data_i=5, explicit=[0, 1, 2, 3, 4, -1, 5, 6, 7, 8, 9, 10])
#e = A.make_alignment(data_i=8, explicit=[0, 1, 2, 5, 6, 7]) #, left=False)
## y e q e d : e s kW      a t
## z   Q e d : e s k u w : a
#e = A.make_alignment(9, explicit=[-1, 0, 4, 5, 6, 7, 8, 9, 10, -1, 11]) #, left=False)
#f = A.make_alignment(9, explicit=[0, -1, 1, 2, 3, 4, 5, 6, 7, 11, -1])
                 # , explicit=[0, 1, 2, 3, 5, -1, -1])
#B = mt.Aligner('ti', 'am', datafiles=["0_4.tr", "2_4.tr", "3_4.tr"])
#c = mt.Alignment(B.data_indices[0], B.data[0], B, explicit=[0, 1, 2, 3, -1, 4])
#B = mt.align.Aligner('am', 'ti', datafiles=["test1.tr"],
#                     test=False, validate=False)

### Test of new data format and Data.read()
##def get_data(vc=0):
##    return mt.Data.read("0_4.tr", None)

### A simple test problem for alignment
##
##TG = [mt.TGroup([mt.Word(["a", "b", "c", "d"]),
##                 mt.Word(["b", "c", "d", "e"]),
##                 mt.Word(["a", "c", "d"]),
##                 mt.Word(["a", "b", "f", "c"]),
##                 mt.Word(["b", "c"]),
##                 mt.Word(["g", "a", "b", "d", "e", "d"]),
##                 mt.Word(["f", "a", "b", "d", "e", "g"])])  
##      ]
##P = mt.Problem(data=TG)
##TG0 = TG[0]
##S = TG0[0]
##T = TG0[1:]

# Problem with 3 translation groups in 5 languages.
# P = mt.Problem('1.tr', name="taTkc;1")
