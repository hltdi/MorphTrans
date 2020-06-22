#! /usr/bin/env python3

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

_version = 0.5

__all__ = ['problem', 'rule', 'word', 'segment', 'align']

from .problem import *

print('\n@@@@ This is MorphTrans, version {} @@@@\n'.format(_version))

def exit(save=True):
    """Exit the program."""
    print("Quitting...")

