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

-- 19.6.2020
   Created.
   Aligner uses the EM algorithm and word translation pairs in two languages
   (for convenience, referred to as 'source' and 'target') to learn alignments
   of characters in the words and character-to-character probabilities.
"""

from .data import *
from .word import *
import numpy

class Aligner:

    # probability of the same character in both languages
    eqprob = 0.5
    # other probabilities
    otherprob = 0.01

    def __init__(self, source, target, datafiles=None):
        self.source = source
        self.target = target
        self.schars = self.get_schars()
        self.tchars = self.get_tchars()
        self.l = len(self.schars)
        self.m = len(self.tchars)
        self.probs = self.make_probs()
        self.init_probs()
        self.tnorms = numpy.zeros(self.m+1)
        # if bidir:
        #     self.tsprobs = self.make_tsprobs()
        #     self.init_probs(False)
        #     self.tscounts = self.make_tscounts()
        self.data = []
        if datafiles:
            for datafile in datafiles:
                d = Data(datafile)
                wpairs = d.get_words([self.source, self.target])
                self.data.extend(wpairs)
            self.counts = self.make_counts()

    def __repr__(self):
        return "⤩:{},{}".format(self.source, self.target)

    def get_schars(self, every=False):
        phones = Phone.PHONES.get(self.source)
        if every:
            return phones[0] + phones[1]
        return phones[0]

    def get_tchars(self, every=False):
        phones = Phone.PHONES.get(self.target)
        if every:
            return phones[0] + phones[1]
        return phones[0]

    def get_tchar_index(self, char, check=False):
        if check:
            if char not in self.tchars:
                return -1
        return self.tchars.index(char)

    def get_schar_index(self, char, check=False):
        if check:
            if char not in self.schars:
                return -1
        return self.schars.index(char)

    def get_schar(self, index):
        if index >= self.l:
            return '0'
        return self.schars[index]

    def get_tchar(self, index):
        if index >= self.m:
            return '0'
        return self.tchars[index]

    def get_prob(self, sindex, tindex):
        return self.probs[sindex][tindex]

    #def get_tsprob(self, sindex, tindex):
    #    return self.tsprobs[tindex][sindex]

    def make_probs(self):
        """
        Make the 'st' table of probabilities: p(s|t) for s/t character
        combinations. source is rows, target columns.
        """
        array = numpy.full((self.l, self.m+1), Aligner.otherprob)
        return array

    # def make_tsprobs(self):
    #     """
    #     Make the 'ts' table of probabilities: p(t|s) for t/s characters
    #     combinations. target is rows, source columns.
    #     """
    #     array = numpy.full((self.m, self.l+1), Aligner.otherprob)
    #     return array

    def init_probs(self):
        """
        Initialize the probabilities in the table.
        """
        array = self.probs
        getrowchar = self.get_schar
        getcolindex = self.get_tchar_index
        for i, row in enumerate(array):
            char = getrowchar(i)
            colindex = getcolindex(char, True)
            if colindex >= 0:
                array[i,colindex] = Aligner.eqprob
        # Normalize
        self.norm_cols()

    def norm_cols(self):
        """
        Normalize each column of the probs array.
        """
        array = self.probs
        for c in range(array.shape[1]):
            self.norm_col(c, array)

    def norm_col(self, ci, array):
        """
        Normalize probabilities in column ci so that they sum to 1.0.
        """
        # sum of pross in column ci
        colsum = array[:,ci].sum()
        # divide each prob in c1 by colsum
        if colsum == 0.0:
            print("Can't normalize 0 column")
        else:
            array[:,ci] /= colsum

    def make_counts(self):
        """
        Make the 'st' table of normalized counts: c(s|t) for s/t character
        combinations. source is rows, target columns.
        """
        array = numpy.zeros((self.l, self.m+1, len(self.data)))
        return array

    # def make_tscounts(self):
    #     """
    #     Make the 'ts' table of normalized counts: c(t|s) for s/t character
    #     combinations. target is rows, source columns.
    #     """
    #     array = numpy.zeros((self.m, self.l+1))
    #     return array

    def update_counts(self):
        """
        Update all of the counts, based on current probs and occurrences
        of chars in word pairs.
        """
        for wpindex, wpair in enumerate(self.data):
            for si in range(self.l):
                denom = self.calc_prob_sum(si)
                for ti in range(self.m+1):
                    self.update_count(si, ti, wpair, wpindex, denom)

    def update_count(self, si, ti, wpair, wpindex, denom=None):
        """
        Calculate the current estimated weighted count of the character
        at index si associated with the character at index ti, given
        ti.
        """
        if not denom:
            denom = self.calc_prob_sum(si)
        sword, tword = wpair
        prob = self.probs[si,ti]
        schar = self.get_schar(si)
        tchar = self.get_tchar(ti)
        stoccs = self.calc_st_occurrences(sword, schar, tword, tchar)
        if stoccs:
            self.counts[si, ti, wpindex] = (prob * stoccs) / denom

    def calc_prob_sum(self, si):
        """
        Calculate the sum of all probs with s = si.
        Sum the si-th row of the probs matrix.
        """
        return self.probs[si].sum()

    def calc_st_occurrences(self, sword, schar, tword, tchar):
        """
        Product of number of schars in sword and number of tchars in tword.
        """
        return sword.count(schar) * tword.count(tchar)

    def update_tnorms(self):
        """
        Using current counts, calculate the normalizing factor
        for each column.
        """
        for tindex in range(self.m+1):
            total = self.counts[:,tindex,:].sum()
            self.tnorms[tindex] = total
