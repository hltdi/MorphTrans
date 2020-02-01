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

class Phone(str):
    """A character or sequence of characters representing a phone type (not an instance) within a Word.
    As of 31.1.2020, not used."""

    # tentative list of phonetic features for phones
    FEATS = ['vowel', 'high', 'low', 'back', 'round', 'palatal', 'stop', 'ejective',
             'voiced', 'fricative', 'liquid', 'nasal', 'lateral', 'labial', 'coronal', 'velar', 'pharyngeal']

    def __init__(self, feats):
        """feats is a list of feature strings for which this segment is True."""
        str.__init__(self)
        self.feats = feats

    def feat_value(self, feat):
        if feat in self.feats:
            return True
        return False

    def distance(self, phone):
        total = 0
        for feat in Phone.FEATS:
            this_value = self.feat_value(feat)
            other_value = phone.feat_value(feat)
            if this_value != other_value:
                total += 1
        return total

