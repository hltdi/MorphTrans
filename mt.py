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

A = mt.Aligner('am', 'ti', datafiles=["0_4.tr", "2_4.tr", "3_4.tr"])
a = A.make_alignment(1)
                 #, explicit=[0, 1, 3, -1, 6, 7, 8, -1, 9])
b = A.make_alignment(0)
                 # , explicit=[0, 1, 2, 3, 5, -1, -1])
#B = mt.Aligner('ti', 'am', datafiles=["0_4.tr", "2_4.tr", "3_4.tr"])
#c = mt.Alignment(B.data_indices[0], B.data[0], B, explicit=[0, 1, 2, 3, -1, 4])
#B = mt.align.Aligner('am', 'ti', datafiles=["test1.tr"],
#                     test=False, validate=False)

# Test of new data format and Data.read()
def get_data(vc=0):
    return mt.Data.read("0_4.tr", None)

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
