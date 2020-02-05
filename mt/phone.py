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

    # characters or character sequences used to represent vowels
    VOWELS = ['a', 'A', 'e', 'A', 'i', 'I', 'o', 'O', 'u', 'U']

    # miscellaneous characters, such as gemination, vowel length, and stress,
    # with low cost for deletion or substitution with vowels
    MISC = ['_', ':']

    DISTANCE = {frozenset({'C', 'C'}): 0.5,
                frozenset({'V', 'V'}): 0.1,
                frozenset({'V', 'C'}): 0.9,
                frozenset({'V', 'M'}): 0.2,
                frozenset({'C', 'M'}): 0.4,
                frozenset({'M', 'M'}): 0.1}

    # cost of inserting or deleting phone
    INSDEL = {'V': 0.5, 'C': 0.8, 'M': 0.1}

    def __init__(self, feats):
        """feats is a list of feature strings for which this segment is True."""
        str.__init__(self)
        self.feats = feats

    def feat_value(self, feat):
        if feat in self.feats:
            return True
        return False

    ## Insert-delete costs

    @staticmethod
    # Cost of inserting or deleting phone.
    def insdel_cost(phone):
        cat = Phone.get_cat(phone)
        return Phone.INSDEL[cat]

    ## Distance between phones

    # Detailed distance, by feature overlap
    def feat_distance(self, phone):
        total = 0
        for feat in Phone.FEATS:
            this_value = self.feat_value(feat)
            other_value = phone.feat_value(feat)
            if this_value != other_value:
                total += 1
        return total

    # Gross distance, consonants, vowels, misc
    # 0 is highest match
    @staticmethod
    def distance(phone1, phone2):
        if phone1 == phone2:
            return 0
        cats = frozenset([Phone.get_cat(phone1), Phone.get_cat(phone2)])
        return Phone.DISTANCE.get(cats)

    # phone categories

    @staticmethod
    def get_cat(phone):
        if phone in Phone.VOWELS:
            return 'V'
        elif phone in Phone.MISC:
            return 'M'
        else:
            return 'C'

    @staticmethod
    def is_vowel(phone):
        return phone in Phone.VOWELS

    @staticmethod
    def is_misc(phone):
        return phone in Phone.MISC

    @staticmethod
    def is_consonant(phone):
        return not Phone.is_vowel(phone) and not Phone.is_misc(phone)

