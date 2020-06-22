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

-- 20.1.2020
   Created.
   Functions to get data from a bilingual data file, to extract word pairs from the
   resulting data.
"""

### USAGE
### To get the data from a translation data file
### >>> data = get_data(languages, file), where
###     languages is a string representing one of the language pair
###     directories, and file is a filename string
### To get pairs of words (list of segment strings), given a data
### dict returned by get_data(),
### >>> words = get_words(data)

import os

def get_data(languages="AmKs", file="1.tr"):
    """
    Given a translation data file, return a dict containing the data, with
    source language root/aspect-voice as keys,
    to be handed to other functions in this module.
    """
    data = {}
    path = os.path.join(languages, file)
    try:
        with open(os.path.join(languages, file), encoding='utf8') as datafile:
            for line in datafile:
                words, misc = line.split(';;')
                misc_split = misc.split(';')
                if len(misc_split) != 3:
                    print("Something wrong with {}".format(line.strip()))
                source_ra, target_ra, features = misc.strip().split(';')
                features = features.split(',')
                sword, twords = words.split(';')
                twords = twords.split(',')
                if source_ra in data:
                    entry = data[source_ra]
                    entry.append((sword, twords, target_ra, features))
                else:
                    data[source_ra] = [(sword, twords, target_ra, features)]
        return data
    except IOError:
        print('No such data file as {}'.format(path))

def get_words(data):
    """
    data is a dict like that returned by get_data().
    Returns a list of source word, target word pairs.
    source word is a list of phones. target word is a list of lists
    of phones (because there can be more than one target word).
    Note: a phone can be represented by more than one character.
    """
    words = []
    for source_ra, entries in data.items():
        for sword, twords, target_ra, features in entries:
            sword = segment(sword, units=SEGUNITS)
            for i, tword in enumerate(twords):
                twords[i] = segment(tword, units=SEGUNITS)
            words.append((sword, twords))
    return words

def segment(word, units=None, correct=True):
    """Given a string representing a word in an ES language, segment it into
    a list of phones. Note: list(word) won't work because some phones are
    represented by two characters."""
    res = []
    pos = 0
    while pos < len(word):
        ch = word[pos]
        if ch in units[0]:
            res.append(ch)
            pos += 1
        else:
            sublists = units[1]
            if ch in sublists:
                sublist = sublists[ch]
                if pos < len(word) - 1 and ch + word[pos + 1] in sublist:
                    if pos < len(word) - 2 and ch + word[pos + 1:pos + 3] in sublist:
                        res.append(ch + word[pos + 1:pos + 3])
                        pos += 3
                    else:
                        res.append(ch + word[pos + 1])
                        pos += 2
                else:
                    res.append(ch)
                    pos += 1
            elif correct:
                print(ch, 'in', word, 'is not an acceptable character')
                return
            else:
                res.append(ch)
                pos += 1
    return res

### Segmentation units for ES languages.

SEGUNITS =  [["a", "e", "E", "i", "I", "o", "u", "@", "A", "w", "y", "'", "`", "|", "_"],
             {"b": ["b", "bW"], "c": ["c", "cW"], "C": ["C", "CW"],
              "d": ["d", "dW"], "f": ["f", "fW"], "g": ["g", "gW", "gY"],
              "h": ["h", "hW", "hY"], "H": ["H", "HW", "HY"], "j": ["j", "jW"],
              "k": ["k", "kW", "kY"], "K": ["K", "KW", "KY"], "l": ["l", "lW"],
              "m": ["m", "mW"], "n": ["n", "nW"], "p": ["p", "pW"], "P": ["P", "PW"],
              "N": ["N", "NW"], "q": ["q", "qW", "qY"], "Q": ["Q", "QW"],
              "r": ["r", "rW"], "s": ["s", "sW"], "S": ["S", "SW"],
              "t": ["t", "tW"], "T": ["T", "TW"], "v": ["v", "vW"],
              "x": ["x", "xW"], "z": ["z", "zW"], "Z": ["Z", "ZW"],
              "^": ["^s", "^S", "^h", "^hW", "^sW", "^SW"]}]

def main():
    pass

if __name__ == "__main__": main()

