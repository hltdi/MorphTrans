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

from .phone import *

class TGroup(list):
    """A list of words that are translations of each other.
    The first word is normally the source, the others targets."""

    def __init__(self, words, languages=None, problem=None):
        list.__init__(self, words)
        # list of language abbreviations
        self.languages = languages
        # instance of Problem for which this TGroup is one piece of data
        self.problem = problem
        self.source = self[0]
        self.targets = self[1:]
        # list of positions corresponding to source positions, one for each target word
        self.alignments = [self.init_alignment(l) for l in range(len(self)-1)]
        # set the number of phone matches for each source word position; to be update
        # after each re-alignment
        self.set_matches()

    def __repr__(self):
        """Print name for the TGroup."""
        return "<<{}...>>".format(self.source.name)

    def init_alignment(self, twordindex):
        """Create an initial alignment for the word at twordindex based on its
        length and similarity of ends to ends of source word."""
        slength = len(self.source)
        target = self.targets[twordindex]
        tlength = len(target)
        a = list(range(slength))
        diff = slength - tlength
        if diff > 0:
            if self.source[0] == target[0]:
                for i in range(tlength-1, slength):
                    a[i] = tlength-1
            elif diff == 1 or self.source[-1] == target[-1]:
                for i in range(0, diff):
                    a[i] = -1
                for i in range(0, tlength):
                    a[i+diff] = i
            else:
                offset = diff // 2
                a = [-1] * slength
                for i in range(tlength):
                    a[i+offset] = i
        elif diff < 0:
            if self.source[-1] == target[-1]:
                for i in range(slength):
                    a[i] = i - diff
            elif diff < -1:
                offset = -diff // 2
                for i in range(slength):
                    a[i] = i + offset
        return a

    def set_matches(self):
        """Set the number of matching phones in target words for each source word position,
        given the current alignments."""
        self.matches = [self.get_matches1(i) for i in range(len(self.source))]

    def get_matches1(self, sindex):
        """Get the number of phones in target words that match the source phone in
        position sindex, given the current alignments."""
        sphone = self.source[sindex]
        count = 0
        for i, (target, alignment) in enumerate(zip(self.targets, self.alignments)):
            alignindex = alignment[sindex]
            if alignindex < 0:
                continue
            else:
                tphone = target[alignindex]
                if tphone == sphone:
                    # a place for degrees of similarity
                    count += 1
        return count

    def align(self, verbosity=0):
        """Repeatedly align source to target words, minimizing edit costs,
        until there are no changes to alignments."""
        change = True
        while change:
            cost, change = self.minimize_all(verbosity=verbosity)

    def minimize_all(self, verbosity=0):
        """Minimize costs once for all target words, adjusting alignments accordingly."""
        total_cost = 0
        change = False
        print("ALIGNING SOURCE TO TARGETS")
        for tword_i in range(len(self.targets)):
            # make a copy of the current alignment for this target word
            current_alignment = self.alignments[tword_i][:]
            cost = self.minimize(tword_i, verbosity=verbosity)
            print("Cost for target {}: {}".format(tword_i, cost))
            total_cost += cost
            if self.alignments[tword_i] != current_alignment:
                change = True
        # Update the matches vector.
        self.set_matches()
        print("TOTAL COST: {}".format(total_cost))
        if change:
            print("SOME ALIGNMENT CHANGED")
        else:
            print("CONVERGED; NO ALIGNMENT CHANGED")
        return total_cost, change

    def minimize(self, tword_index, verbosity=0):
        """Minimize edit costs for the source word and target word at tword_index,
        returning the cost and changing the alignment if necessary."""
        source = self.source
        target = self.targets[tword_index]
        alignment = self.alignments[tword_index]
        sourcecopy = source.copy()
#        if verbosity:
        print()
        print("Aligning {} with {}".format(source, target))
        return self.minimize1(source, target, len(source)-1, len(target)-1,
                              alignment, sourcecopy, verbosity=verbosity)

    def minimize1(self, source, target, sindex, tindex, alignment,
                  sourcecopy=None, verbosity=0):
        """Return the cost of editing source word to produce target word,
        ending in positions sindex in source, tindex in target. Change the current
        alignment where there are deletions or substitutions.
        Start at the end of source and target words, recursively figuring the
        cost, using the expression on p. 1187 of Durrett & DeNero."""
        if sindex < 0 and tindex < 0:
            return 0
        if not sourcecopy:
            sourcecopy = source.copy()
        costs = self.costs(source, target, sindex, tindex)
        if verbosity:
            print(" Found costs {} at si {}, ti {}".format(costs, sindex, tindex))
        mincost = min(costs, key=lambda d: d[0])
        cost, (newsindex, newtindex), oper = mincost
        # Change sourcecopy to reflect operation selected (not needed?)
        # and update alignment.
        if oper == 1:
            # Deletion
            if verbosity > 1:
                print("  Deleting {}".format(source[sindex]))
            alignment[sindex] = -1
            sourcecopy.delete(sindex)
        else:
            phone = target[tindex]
            if oper == 0:
                # Insertion; this does *not* change the alignment, it seems
                if verbosity > 1:
                    print("  Inserting {}".format(phone))
                sourcecopy.insert(sindex, phone)
            else:
                # Substitution
                alignment[sindex] = tindex
                sphone = source[sindex]
                if sphone == phone:
                    if verbosity > 1:
                        print("  No change for {}".format(phone))
                else:
                    sourcecopy.substitute(sindex, phone)
                    if verbosity > 1:
                        print("  Substituting {} for {}".format(phone, sphone))
        if verbosity:
            print(" Cost: {}, si* {}, ti* {}".format(cost, newsindex, newtindex))
            print(" Current transformed word: {}".format(sourcecopy.chars()))
            print(" Current alignment: {}".format(alignment))
        return cost + self.minimize1(source, target, newsindex, newtindex,
                                     alignment, sourcecopy=sourcecopy, verbosity=verbosity)

    def insert_cost(self, sindex=-1, ldiff=0, tindex=-1, tphone='', scontext=None):
        """
        Return the cost of inserting tphone in position sindex in source,
        the new indices to check next,
        and 0 (id for insertion).
        """
        return (self.problem.score_insert(sindex, tphone, ldiff, scontext), (sindex, tindex-1), 0)

    def delete_cost(self, sindex=-1, sphone='', ldiff=0, tindex=-1, scontext=None):
        """
        Return the cost of deleting sphone in position sindex in source,
        the new indices to check next,
        and 1 (id for deletion).
        """
        return (self.problem.score_delete(sindex, sphone, ldiff), (sindex-1, tindex), 1)

    def compare_cost(self, sindex=-1, sphone='', tindex=-1, tphone='',
                     ldiff=0, scontext=None, tcontext=None, matches=0):
        """
        Return the cost of keeping sphone (if it's the same as tphone) or replacing it with tphone (if they're different),
        the new indices to check next,
        and 2 (id for substitution).
        """
        return (self.problem.score_compare(sindex, tphone, sphone, scontext, tcontext, matches), (sindex-1, tindex-1), 2)

    def costs(self, source, target, sindex, tindex):
        """
        Return the costs associated with each operation, along with updated source and target indices, and indices
        to distinguish the operations.
        """
        # difference in length of source and target words
        ldiff = len(source) - len(target)
        # current source and target phones in positions sindex and tindex
        sphone = source[sindex] if sindex >= 0 else ''
        tphone = target[tindex] if tindex >= 0 else ''
        # preceding contexts of in source and target words
        scontext = source[max([0,sindex-2]):sindex] if sindex >= 1 else []
        tcontext = target[max([0,tindex-2]):tindex] if tindex >= 1 else []
        # number of other target words with the current source phone in positions aligned with sindex
        matches = self.matches[sindex] - 1 if sindex >= 0 else 0
        if sindex < 0:
            # insertion is the only possibility because we're reached the beginning of the source word
            return [self.insert_cost(sindex=sindex, ldiff=ldiff, tindex=tindex, tphone=tphone, scontext=scontext)]
        if tindex < 0:
            # deletion is the only possibility because we've reached the beginning of the target word
            return [self.delete_cost(sindex=sindex, sphone=sphone, ldiff=ldiff, tindex=tindex, scontext=scontext)]
        # all three edit operations are possible in these source and target positions
        return [self.insert_cost(sindex=sindex, ldiff=ldiff, tindex=tindex, tphone=tphone, scontext=scontext),
                self.delete_cost(sindex=sindex, sphone=sphone, ldiff=ldiff, tindex=tindex, scontext=scontext),
                self.compare_cost(sindex=sindex, sphone=sphone, tindex=tindex, tphone=tphone,
                                  ldiff=ldiff, scontext=scontext, tcontext=tcontext, matches=matches)]

class Word(list):
    """A list of characters (strings) or Phone instances."""

    def __init__(self, phones, name=''):
        list.__init__(self, phones)
        self.name = name or ''.join(self)

    def __repr__(self):
        return "<{}>".format(self.name)

    def chars(self):
        return ''.join(self)

    # String editing

    def copy(self):
        """A copy of the Word."""
        return Word(self, name="*" + self.name)

    def insert(self, index, phone):
        self[index+1:index+1] = [phone]

    def delete(self, index):
        del self[index]

    def substitute(self, index, phone):
        self[index] = phone
        
