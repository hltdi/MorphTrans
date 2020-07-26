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

    PHONES = \
      {'am':
        [[
        'a', 'e', 'E', 'i', 'o', 'u', \
        'w', 'y', "'", ':', 'b', 'bW', 'c', 'cW', 'C', 'CW', 'd', 'dW', \
        'f', 'fW', 'g', 'gW', 'h', 'hW', 'j', 'jW', 'k', 'kW', 'l', 'lW', \
        'm', 'mW', 'n', 'nW', 'N', 'NW', 'q', 'qW', \
        'r', 'rW', 's', 'sW', 'S', 'SW', 't', 'tW', 'T', 'TW', \
        'x', 'xW', 'z', 'zW', 'Z', 'ZW'
        ],
        [
        'I', 'p', 'pW', 'P', 'PW', 'v', 'vW', \
        '^s', '^S', '^h', '^hW', '^sW', '^SW'
        ]],
        'ti':
        [[
        'a', 'e', 'E', 'i', 'o', 'u', \
        'w', 'y', "'", '`', ':', 'b', 'bW', 'c', 'cW', 'C', 'CW', 'd', 'dW', \
        'f', 'fW', 'g', 'gW', 'h', 'hW', 'H', 'HW', 'j', 'jW', \
        'k', 'kW', 'K', 'KW', 'l', 'lW', 'm', 'mW', 'n', 'nW', 'N', 'NW', \
        'q', 'qW', 'Q', 'QW', 'r', 'rW', 's', 'sW', 'S', 'SW', 't', 'tW', \
        'T', 'TW', 'x', 'xW', 'z', 'zW', 'Z', 'ZW'
         ],
        [
        'I', 'p', 'pW', 'P', 'PW', 'v', 'vW', \
        '^s', '^S', '^h', '^hW', '^sW', '^SW'
        ]],
        'amG':
        [[
          "ሀ", "ሁ", "ሂ", "ሃ", "ሄ", "ህ", "ሆ",\
          "ለ", "ሉ", "ሊ", "ላ", "ሌ", "ል", "ሎ", "ሏ",\
          "ሐ", "ሑ", "ሒ", "ሓ", "ሔ", "ሕ", "ሖ",\
          "መ", "ሙ", "ሚ", "ማ", "ሜ", "ም", "ሞ", "ሟ",\
          "ሠ", "ሡ", "ሢ", "ሣ", "ሤ", "ሥ", "ሦ", "ሧ",\
          "ረ", "ሩ", "ሪ", "ራ", "ሬ", "ር", "ሮ", "ሯ",\
          "ሰ", "ሱ", "ሲ", "ሳ", "ሴ", "ስ", "ሶ", "ሷ",\
          "ሸ", "ሹ", "ሺ", "ሻ", "ሼ", "ሽ", "ሾ", "ሿ",\
          "ቀ", "ቁ", "ቂ", "ቃ", "ቄ", "ቅ", "ቆ", "ቈ", "ቋ", "ቍ",\
          "በ", "ቡ", "ቢ", "ባ", "ቤ", "ብ", "ቦ", "ቧ",\
          "ተ", "ቱ", "ቲ", "ታ", "ቴ", "ት", "ቶ", "ቷ",\
          "ቸ", "ቹ", "ቺ", "ቻ", "ቼ", "ች", "ቾ", "ቿ",\
          "ኀ", "ኁ", "ኂ", "ኃ", "ኄ", "ኅ", "ኆ", "ኈ", "ኋ","ኍ",\
          "ነ", "ኑ", "ኒ", "ና", "ኔ", "ን", "ኖ", "ኗ",\
          "ኘ", "ኙ", "ኚ", "ኛ", "ኜ", "ኝ", "ኞ", "ኟ",\
          "አ", "ኡ", "ኢ", "ኣ", "ኤ", "እ", "ኦ",\
          "ከ", "ኩ", "ኪ", "ካ", "ኬ", "ክ", "ኮ", "ኰ", "ኳ", "ኵ",\
          "ኸ", "ኹ", "ኺ", "ኻ", "ኼ", "ኽ", "ኾ", "ዀ",\
          "ወ", "ዉ", "ዊ", "ዋ", "ዌ", "ው", "ዎ",\
          "ዐ", "ዑ", "ዒ", "ዓ", "ዔ", "ዕ", "ዖ",\
          "ዘ", "ዙ", "ዚ", "ዛ", "ዜ", "ዝ", "ዞ", "ዟ",\
          "ዠ", "ዡ", "ዢ", "ዣ", "ዤ", "ዥ", "ዦ", "ዧ",\
          "የ", "ዩ", "ዪ", "ያ", "ዬ", "ይ", "ዮ",\
          "ደ", "ዱ", "ዲ", "ዳ", "ዴ", "ድ", "ዶ", "ዷ",\
          "ጀ", "ጁ", "ጂ", "ጃ", "ጄ", "ጅ", "ጆ", "ጇ",\
          "ገ", "ጉ", "ጊ", "ጋ", "ጌ", "ግ", "ጎ", "ጐ", "ጓ", "ጕ",\
          "ጠ", "ጡ", "ጢ", "ጣ", "ጤ", "ጥ", "ጦ", "ጧ",\
          "ጨ", "ጩ", "ጪ", "ጫ", "ጬ", "ጭ", "ጮ", "ጯ",\
          "ጸ", "ጹ", "ጺ", "ጻ", "ጼ", "ጽ", "ጾ", "ጿ",\
          "ፀ", "ፁ", "ፂ", "ፃ", "ፄ", "ፅ", "ፆ",\
          "ፈ", "ፉ", "ፊ", "ፋ", "ፌ", "ፍ", "ፎ", "ፏ"
      ],
      [
          "ሗ", "ቊ", "ኊ", "ኲ", "ዂ", "ዃ", "ዄ", "ዅ", "ጒ", "ቌ", "ኌ", "ኴ", "ጔ",\
          "ቨ", "ቩ", "ቪ", "ቫ", "ቬ", "ቭ", "ቮ", "ቯ",\
           "ኧ",\
          "ጰ", "ጱ", "ጲ", "ጳ", "ጴ", "ጵ", "ጶ", "ጷ",\
          "ፐ", "ፑ", "ፒ", "ፓ", "ፔ", "ፕ", "ፖ", "ፗ"
      ]],

      'tiG':
      [[
          "ሀ", "ሁ", "ሂ", "ሃ", "ሄ", "ህ", "ሆ", "ኋ",\
          "ለ", "ሉ", "ሊ", "ላ", "ሌ", "ል", "ሎ", "ሏ",\
          "ሐ", "ሑ", "ሒ", "ሓ", "ሔ", "ሕ", "ሖ", "ሗ",\
          "መ", "ሙ", "ሚ", "ማ", "ሜ", "ም", "ሞ", "ሟ",\
          "ሠ", "ሡ", "ሢ", "ሣ", "ሤ", "ሥ", "ሦ", "ሧ",\
          "ረ", "ሩ", "ሪ", "ራ", "ሬ", "ር", "ሮ", "ሯ",\
          "ሰ", "ሱ", "ሲ", "ሳ", "ሴ", "ስ", "ሶ", "ሷ",\
          "ሸ", "ሹ", "ሺ", "ሻ", "ሼ", "ሽ", "ሾ", "ሿ",\
          "ቀ", "ቁ", "ቂ", "ቃ", "ቄ", "ቅ", "ቆ", "ቈ", "ቊ", "ቋ", "ቌ", "ቍ",\
          "ቐ", "ቑ", "ቒ", "ቓ", "ቔ", "ቕ", "ቖ", "ቘ", "ቚ", "ቛ", "ቜ", "ቝ",\
          "በ", "ቡ", "ቢ", "ባ", "ቤ", "ብ", "ቦ", "ቧ",\
          "ተ", "ቱ", "ቲ", "ታ", "ቴ", "ት", "ቶ", "ቷ",\
          "ቸ", "ቹ", "ቺ", "ቻ", "ቼ", "ች", "ቾ", "ቿ",\
          "ነ", "ኑ", "ኒ", "ና", "ኔ", "ን", "ኖ", "ኗ",\
          "ኘ", "ኙ", "ኚ", "ኛ", "ኜ", "ኝ", "ኞ", "ኟ",\
          "አ", "ኡ", "ኢ", "ኣ", "ኤ", "እ", "ኦ", "ኧ",\
          "ከ", "ኩ", "ኪ", "ካ", "ኬ", "ክ", "ኮ", "ኰ", "ኲ", "ኳ", "ኴ", "ኵ",\
          "ኸ", "ኹ", "ኺ", "ኻ", "ኼ", "ኽ", "ኾ", "ዀ", "ዂ", "ዃ", "ዄ", "ዅ",\
          "ወ", "ዉ", "ዊ", "ዋ", "ዌ", "ው", "ዎ",\
          "ዐ", "ዑ", "ዒ", "ዓ", "ዔ", "ዕ", "ዖ",\
          "ዘ", "ዙ", "ዚ", "ዛ", "ዜ", "ዝ", "ዞ", "ዟ",\
          "ዠ", "ዡ", "ዢ", "ዣ", "ዤ", "ዥ", "ዦ", "ዧ",\
          "የ", "ዩ", "ዪ", "ያ", "ዬ", "ይ", "ዮ",\
          "ደ", "ዱ", "ዲ", "ዳ", "ዴ", "ድ", "ዶ", "ዷ",\
          "ጀ", "ጁ", "ጂ", "ጃ", "ጄ", "ጅ", "ጆ", "ጇ",\
          "ገ", "ጉ", "ጊ", "ጋ", "ጌ", "ግ", "ጎ", "ጐ", "ጒ", "ጓ", "ጔ", "ጕ",\
          "ጠ", "ጡ", "ጢ", "ጣ", "ጤ", "ጥ", "ጦ", "ጧ",\
          "ጨ", "ጩ", "ጪ", "ጫ", "ጬ", "ጭ", "ጮ", "ጯ",\
          "ጸ", "ጹ", "ጺ", "ጻ", "ጼ", "ጽ", "ጾ", "ጿ",\
          "ፀ", "ፁ", "ፂ", "ፃ", "ፄ", "ፅ", "ፆ",\
          "ፈ", "ፉ", "ፊ", "ፋ", "ፌ", "ፍ", "ፎ", "ፏ"
          ],
      [
          "ኀ", "ኁ", "ኂ", "ኃ", "ኄ", "ኅ", "ኆ", "ኈ", "ኊ", "ኌ", "ኍ",\
          "ቨ", "ቩ", "ቪ", "ቫ", "ቬ", "ቭ", "ቮ", "ቯ",\
          "ጰ", "ጱ", "ጲ", "ጳ", "ጴ", "ጵ", "ጶ", "ጷ",\
          "ፐ", "ፑ", "ፒ", "ፓ", "ፔ", "ፕ", "ፖ", "ፗ"
      ]]
      }

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
