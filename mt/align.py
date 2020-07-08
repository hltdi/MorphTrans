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
    # null character
    NULL = '0'
    # number of length differences for parameters
    nlendiffs = 16
    # smallest length difference for parameters
    lendiff0 = -4
    # arbitrary constant to generate seed so that datasets are constant
    seedgen = 4
    # proportion of data to use for test set
    test_frac = 0.1
    # proprtion of data to use for validation set
    validation_frac = 0.1

    # parameter controlling effect of relative distance between
    # source and target indices on connection probability
    #delta = 0.9
    #sigma = 0.5

    def __init__(self, source, target, datafiles=None, null=True):
        self.source = source
        self.target = target
        # whether there is a target NULL position
        self.null = null
        self.schars = self.get_schars()
        self.tchars = self.get_tchars()
        self.ns = len(self.schars)
        self.nt = len(self.tchars)
        self.probs = self.make_probs()
        self.init_probs()
        self.tnorms = numpy.zeros(self.nt+1)
        # history of sum of squares of changes in probs during M-step
        self.changes = []
        # Deletion parameters
        self.deltas = self.make_del_counts()
        # Justification/alignment parameters
        self.lambdas = self.make_ljust_probs()
        # List of word pairs, consisting of character strings
        self.data = []
        # list of word pairs, consisting of lists of character indices
        self.data_indices = []
        if datafiles:
            for datafile in datafiles:
                d = Data(datafile)
                wpairs = d.get_words([self.source, self.target])
                self.data.extend(wpairs)
            # Shuffle data so that it can be split into training, testing
            # and validation subsets
            self.shuffle(self.data)
            self.test, self.validation, self.training = self.split_data(self.data)
            # Or only work with indexed data
            self.data_indices = [self.wpair2indices(wp) for wp in self.data]
            self.test_indices = [self.wpair2indices(wp) for wp in self.test]
            self.validation_indices = [self.wpair2indices(wp) for wp in self.validation]
            self.training_indices = [self.wpair2indices(wp) for wp in self.training]
            self.counts = self.make_counts()
            self.ljustcounts = self.make_just_counts()
            self.rjustcounts = self.make_just_counts()
            self.diffcounts = self.make_diff_counts()

    def __repr__(self):
        return "⤩:{},{}".format(self.source, self.target)

    ### datasets

    def shuffle(self, data, reproduce=True):
        """
        Shuffle items in dataset so they can be split into training, test,
        validation sets.
        """
        if reproduce:
            numpy.random.seed(Aligner.seedgen)
        numpy.random.shuffle(data)

    def split_data(self, data):
        """Split data into test, validation, and training sets."""
        n = len(data)
        test_i = round(n * Aligner.test_frac)
        validation_i = test_i + round(n * Aligner.validation_frac)
        return data[:test_i], data[test_i:validation_i], data[validation_i:]

    ### Accessing characters, probabilities, indices

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
        if index >= self.ns:
            return Aligner.NULL
        return self.schars[index]

    def get_tchar(self, index):
        if index >= self.nt:
            return Aligner.NULL
        return self.tchars[index]

    def get_prob(self, sindex, tindex):
        return self.probs[sindex][tindex]

    def get_schar_probs(self, char):
        charindex = self.get_schar_index(char)
        return self.probs[charindex]

    def get_tchar_probs(self, char):
        charindex = self.get_tchar_index(char)
        return self.probs[:,charindex]

    def get_char_prob(self, schar, tchar):
        sindex = self.get_schar_index(schar)
        tindex = self.get_tchar_index(tchar)
        return self.get_prob(sindex, tindex)

    def get_target_probs(self, tword, sindex):
        """
        Given a sindex, find the st probabilities of all characters
        in tword, with NULL in last position.
        """
        tprobs = [self.probs[sindex,tindex] for tindex in tword]
        if self.null:
            tprobs += [self.probs[sindex,self.nt]]
        return numpy.array(tprobs)

    def get_del_prob(self, sword, tword):
        """
        Get the probability that a character in source word will be
        deleted (connected to the NULL character).
        """
        m = len(sword)
        l = len(tword)
        count = self.deltas[m-l]
        return count / m

    def get_ljust_prob(self, sword, tword):
        m = len(sword)
        l = len(tword)
        index = self.get_lendiff_index(m-l)
        return self.lambdas[index]

    def get_NULL_prob(self, si):
        return self.probs[si,self.nt]

    def get_lendiff(self, index):
        """Return the len(s) - len(t) difference associated with the
        count or parameter in the index position."""
        return Aligner.lendiff0 + index

    def get_lendiff_index(self, diff):
        """
        Return the index associated with the len(s) - len(t) difference
        in the count or parameter table.
        """
        return diff - Aligner.lendiff0

    def get_wpair_lendiff(self, wpair):
        """
        Return the difference in lengths between the source and target
        words in wpair.
        """
        return len(wpair[0]) - len(wpair[1])

    def wpair2indices(self, wpair):
        """
        Convert words in wpair to lists of indices into arrays.
        """
        sindices = [self.get_schar_index(c) for c in wpair[0]]
        tindices = [self.get_tchar_index(c) for c in wpair[1]]
        return sindices, tindices

    def targetN(self):
        """
        Number of target probs and counts, depending on whether
        there is a NULL position with targets.
        """
        if self.null:
            return self.nt+1
        return self.nt

    #def get_tsprob(self, sindex, tindex):
    #    return self.tsprobs[tindex][sindex]

    ### Character positions within words

    ### Arrays and parameters

    def make_probs(self):
        """
        Make the 'st' table of probabilities: p(s|t) for s/t character
        combinations. source is rows, target columns.
        """
        array = numpy.full((self.ns, self.targetN()), Aligner.otherprob)
        return array

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
        array = numpy.zeros((self.ns, self.targetN(), len(self.training)))
        return array

    def make_just_counts(self):
        """
        Make the left- or right-justified table of normalized counts, one for each
        combination of len(s) - len(t) distance."""
        array = numpy.zeros((Aligner.nlendiffs))
        return array

    def make_diff_counts(self):
        """
        Count the number of pairs for each len(s) - len(t) difference.
        """
        array = numpy.zeros((Aligner.nlendiffs))
        for wpair in self.training:
            diff = self.get_wpair_lendiff(wpair)
            index = self.get_lendiff_index(diff)
            array[index] += 1
        return array

    def make_del_counts(self):
        """
        Make the parameters governing how many characters in a source
        sentence will be deleted (connected to the NULL character), given
        different differences between source and target word lengths.
        """
        return {-5: 0.0, -4: 0.0, -3: 0.01, -2: 0.05, -1: 0.5,
                0: 0.0,
                1: 1.3, 2: 2.5, 3: 3.7, 4: 5.0, 5: 6.3,
                6: 7.5, 7: 9.0, 8: 10.3, 9: 11.5, 10: 13.0, 11: 14.3}

    def make_ljust_probs(self):
        """
        Make the probabilities of left alignment, given different
        differences between source and target word lengths.
        """
        array = numpy.full((Aligner.nlendiffs), 0.5)
        return array

        # return {-5: 0.5, -4: 0.5, -3: 0.5, -2: 0.5, -1: 0.5,
        #         0: 0.5,
        #         1: 0.5, 2: 0.5, 3: 0.5, 4: 0.5,
        #         5: 0.5, 6: 0.5, 7: 0.5, 8: 0.5,
        #         9: 0.5, 10: 0.5, 11: 0.5}

    def init_lambdas(self):
        """
        Initialize the lambdas (left justification parameters).
        """
        for key in self.lambdas.keys():
            self.lambdas[key] = 0.5

    def init_deltas(self):
        """
        Initialize the deltas (character deletion parameters).
        """
        for key in self.deltas.keys():
            self.deltas[key] = 0.5

    ### EM

    def EM(self, iter_cutoff=10, error_cutoff=-10.0):
        iteration = 0
        changes = 1.0
        while iteration < iter_cutoff and changes > error_cutoff:
            print("ITERATION {}".format(iteration+1))
            print("E step")
            self.e_init()
            self.update_counts()
            print("M step")
            self.update_tnorms()
            self.update_probs()
            self.update_lambdas()
            iteration += 1
            changes = self.changes[-1]

    def reinit(self):
        """
        Reinitialize parameters for a new EM run.
        """
        self.init_probs()
        self.init_lambdas()

    ## E-step

    def e_init(self):
        """
        Do whatever initialization is required before the E-step on each
        EM iteration.
        """
        # Reset ljust and rjust counts to 0
        self.ljustcounts[:] = 0.0
        self.rjustcounts[:] = 0.0

    def update_counts(self):
        """
        On the E-step, update all of the expected counts, based on current probs
        and cooccs of chars in word pairs.
        """
        for wpindex, wpair in enumerate(self.training_indices):
            lendiff = self.get_wpair_lendiff(wpair)
            diffindex = self.get_lendiff_index(lendiff)
            #left = 0.0
            #right = 0.0
            for si in range(self.ns):
                denom = self.calc_prob_sum(si, wpair[1])
                #left0 = 0.0
                #right0 = 0.0
                for ti in range(self.nt):
                    self.update_count(si, ti, wpair, wpindex,
                                      denom, diffindex)
                    #left += left1
                    #right += right1
                #if left0:
                    #print("wp {} si {} left {} right {}".format(wpindex, si, left0, right0))
                if self.null:
                    # Figure the deletion count separately
                    self.update_del_count(si, wpair, wpindex, denom)
            #self.ljustcounts[diffindex] += left
            #self.rjustcounts[diffindex] += right
            #print("wp {} left {} right {} index {}".format(wpindex, left, right, diffindex))

    def update_count(self, si, ti, wpair, wpindex, denom=None,
                     diffindex=0):
        """
        Calculate the current expected weighted count of the character
        index si associated with the character index ti, given
        ti.
        And increment the left and right justification counts.
        """
        sword, tword = wpair
        if not denom:
            # Denominator is the s|t probability normalized by the
            # s|Σt_i for all i in target word.
            denom = self.calc_prob_sum(si, tword)
        prob = self.probs[si,ti]
        #schar = self.get_schar(si)
        #tchar = self.get_tchar(ti)
        leftprox, rightprox = self.calc_cooccs(sword, si, tword, ti)
        if leftprox or rightprox:
            cooccs = (leftprox + rightprox) / 2.0
            self.counts[si, ti, wpindex] = (prob * cooccs) / denom
            # also normalize these?
            self.ljustcounts[diffindex] += (prob * leftprox)
            self.rjustcounts[diffindex] += (prob * rightprox)
        #return leftprox, rightprox

    def calc_prob_sum(self, si, tword):
        """
        Calculate the sum of all probs with s = si and ti in tword.
        """
        probs = [self.probs[si,ti] for ti in tword]
        return sum(probs)

    def update_del_count(self, si, wpair, wpindex, denom=None):
        """
        Calculate the current expected weighted count of deletions
        of character si for word pair wpair.
        """
        sword, tword = wpair
        if not denom:
            denom = self.calc_prob_sum(si, tword)
        prob = self.get_NULL_prob(si)
        scount = sword.count(si)
        delcount = self.get_del_count(len(sword), len(tword))
        cooccs = scount * delcount
        if cooccs:
            self.counts[si, self.nt, wpindex] = (prob * cooccs) / denom

    def get_del_count(self, slength, tlength):
        """
        Get the expected number of deleted characters in source word of
        length slength associated with target word of length tlength.
        """
        diff = slength - tlength
        param = self.deltas[diff]
        return param / slength

    def calc_cooccs(self, sword, schar, tword, tchar):
        """
        Product of number of schars in sword and number of tchars in tword
        weighted by left justification probabilities applied to characters'
        positions.
        Return both left and right estimates.
        """
        if schar not in sword or tchar not in tword:
            return 0.0, 0.0
        total = 0.0
        lbd = self.get_ljust_prob(sword, tword)
        slength = len(sword)
        tlength = len(tword)
        maxsep = max([slength, tlength])
        leftproxsum = 0.0
        rightproxsum = 0.0
        for j, s in enumerate(sword):
            if s != schar:
                continue
            for i, t in enumerate(tword):
                if t != tchar:
                    continue
                # characters at positions j, j match schar, tchars
                sep = j - i
                leftsep = abs(sep)
                rightsep = abs(j-slength - (i - tlength))
                leftsepnorm = leftsep / maxsep
                rightsepnorm = rightsep / maxsep
                leftprox = lbd * (1.0 - leftsepnorm)
                rightprox = (1.0 - lbd) * (1.0 - rightsepnorm)
                total += leftprox + rightprox
                # update justification make_counts
                leftproxsum += leftprox
                rightproxsum += rightprox
        # self.ljustcounts[diffindex] += leftproxsum
        # self.rjustcounts[diffindex] += rightproxsum
        return leftproxsum, rightproxsum
#        return sword.count(schar) * tword.count(tchar)

    ## M-step

    def update_tnorms(self):
        """
        Using current counts, calculate the normalizing factor
        for each column (λ_e in Brown et al.).
        """
        for tindex in range(self.targetN()):
            total = self.counts[:,tindex,:].sum()
            self.tnorms[tindex] = total

    def update_probs(self):
        """
        Update all of the probabilities in the M-step.
        """
        # This might not need to be
        changes = 0.0
        for sindex in range(self.ns):
            for tindex in range(self.targetN()):
                oldprob = self.probs[sindex,tindex]
                total = self.counts[sindex,tindex,:].sum()
                if total:
                    self.probs[sindex,tindex] = total / self.tnorms[tindex]
                else:
                    self.probs[sindex,tindex] = 0.0
                change = numpy.power(self.probs[sindex,tindex] - oldprob, 2)
                changes += change
        changes /= self.probs.size
        changes = numpy.log(changes)
        self.changes.append(changes)
        print("Changes: {}".format(changes))

    def update_lambdas(self):
        """
        Update all left justification probabilities.
        """
        for lendiff in range(Aligner.nlendiffs):
            # Find the left and right justification counts for this lendiff
            leftcount = self.ljustcounts[lendiff]
            rightcount = self.rjustcounts[lendiff]
            if leftcount or rightcount:
                self.lambdas[lendiff] = leftcount / (leftcount + rightcount)

    ### Relative positions within source and target Words

    # @staticmethod
    # def pdf(x, mu=0.0, sigma=0.0):
    #     sigma = sigma or Aligner.sigma
    #     x = float(x - mu) / sigma
    #     return numpy.exp(-x*x/2.0) / numpy.sqrt(2.0*numpy.pi) / sigma

    @staticmethod
    def rel_index(index, length):
        """
        Return the relative position of index within word.
        IF index = 0 => 0.0
        IF index = len(word)-1 => 1.0
        """
        # if index < 0:
        #     # convert negative index to positive
        #     index = length + index
        # n = length-1
        return index / (length - 1)

    @staticmethod
    def rel_sep(length1, index1, length2, index2):
        """
        The relative left-justified separation of index1 in word of length1
        and index2 in word of length2.
        0.0 <= rvalue <= 1.0
        """
        relindex1 = Aligner.rel_index(index1, length1)
        relindex2 = Aligner.rel_index(index2, length2)
        return abs(relindex1 - relindex2)

    @staticmethod
    def right_rel_sep(length1, index1, length2, index2):
        """
        The relative right-justified separation of index1 in word of length1
        and index2 in word of length2.
        0.0 <= rvalue <= 1.0
        """
        relindex1 = Aligner.rel_index(length1 - index1, length1)
        relindex2 = Aligner.rel_index(length2 - index2, length2)
        return abs(relindex1 - relindex2)

    # def scaled_rel_proximity(length1, index1, length2, index2):
    #     proximity = Aligner.rel_proximity(length1, index1, length2, index2)
    #     return Aligner.pdf(proximity)

    @staticmethod
    def rel_proximity_sum(sindex, sword, tword):
        l = len(sword)
        m = len(tword)
        return sum([Aligner.scaled_rel_proximity(l, sindex, m, tindex) \
        for tindex in range(m)])

    @staticmethod
    def rel_proximity_prob(sindex, sword, tindex, tword, total=0):
        if not total:
            total = Aligner.rel_proximity_sum(sindex, sword, tword)
        l = len(sword)
        m = len(tword)
        proximity = Aligner.scaled_rel_proximity(l, sindex, m, tindex)
        return proximity/l

    # @staticmethod
    # def rel_proximity(word1, index1, word2, index2, length1=0, length2=0):
    #     """
    #     Relative proximity of position index1 in word1 and index2
    #     in word2, scaled by delta.
    #     """
    #     length1 = length1 or len(word1)
    #     length2 = length2 or len(word2)
    #     distance = \
    #     Aligner.rel_dist(length1, index1, length2, index2)
    #     return 1.0 - Aligner.delta * distance

    # @staticmethod
    # def rel_proximity_sum(sindex, sword, tword):
    #     l = len(sword)
    #     m = len(tword)
    #     return sum([Aligner.rel_proximity(sword, sindex, tword, tindex, l, m) \
    #     for tindex in range(m)])

    ### Alignments for word pairs.

    # @staticmethod
    # def get_NULL_tprob(tprobs):
    #     """
    #     Return the probability associated with tchar = NULL
    #     within the vector of probs associated with a given sindex.
    #     """
    #     # or tprobs[-1]
    #     return tprobs[self.nt]

    @staticmethod
    def argmaxes(array):
        """
        Utility function. Like argmax except that it returns a list of
        positions in case there's a tie.
        """
        indices = [0]
        maximum = array[0]
        for i, a in enumerate(array[1:]):
            if a > maximum:
                maximum = a
                indices = [i+1]
            elif a == maximum:
                indices.append(i+1)
        return indices

    def find_best_alignment(self, wpair):
        sword, tword = wpair
        alignment = []
        conflicts = []
        for scindex, sindex in enumerate(sword):
            tprobs = self.get_target_probs(tword, sindex)
            maxtcindices = Aligner.argmaxes(tprobs)
            if len(maxtcindices) > 1:
                # the best target char appears more than once;
                # pick the closest one
                print("Multiple maxtcindices for {}|{}: {}".format(scindex, sindex, maxtcindices))
                stdists = [(abs(scindex - tcindex), tcindex) for tcindex in maxtcindices]
                stdists.sort(key=lambda x: x[0])
                besttcindex = stdists[0][1]
            else:
                besttcindex = maxtcindices[0]
            print("Best tc index for {}|{}: {}".format(scindex, sindex, besttcindex))
            if besttcindex in alignment:
                prevscindex = alignment.index(besttcindex)
                conflicts.append(((prevscindex, scindex), besttcindex))
#                if self.keep_old_stconn(sword, tword, scindex, prevscindex, besttcindex):
                    # Previous connection is better; connect this scindex to NONE
#                    alignment.append(len(tword))
#                else:
                    # New connection is better; replace old one with NONE
#                    alignment[prevscindex] = len(tword)
                    #alignment.append(besttcindex)
            alignment.append(besttcindex)
        return alignment, conflicts

    def keep_old_stconn(self, sword, tword, newsci, oldsci, tci):
        """
        When two source char indices are to be connected to the same
        target char index (other than NONE), pick the scindex that is closer
        in relative index to the target char index.
        """
        olddist = Aligner.get_rel_index_distance(sword, oldsci, tword, tci)
        newdist = Aligner.get_rel_index_distance(sword, newsci, tword, tci)
        return olddist < newdist
