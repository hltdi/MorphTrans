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

TG = [mt.TGroup([mt.Word(["a", "b", "c", "d"]),
                 mt.Word(["b", "c", "d", "e"]),
                 mt.Word(["a", "c", "d"]),
                 mt.Word(["a", "b", "f", "c"]),
                 mt.Word(["b", "c"]),
                 mt.Word(["g", "a", "b", "d", "e", "d"]),
                 mt.Word(["f", "a", "b", "d", "e", "g"])])  
      ]

P = mt.Problem(data=TG)

TG0 = TG[0]

S = TG0[0]
T = TG0[1:]