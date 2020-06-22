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
    Any word in the list can be treated as the source, the others
    as targets."""

    # Number of positions before current index in source or target word to check
    # for presence of current source or target phone
    context_size = 3

    def __init__(self, words, languages=None, problem=None, comments=None):
        list.__init__(self, words)
        # list of language abbreviations
        self.languages = languages
        # instance of Problem for which this TGroup is one piece of data
        self.problem = problem
        self.source = self[0]
        self.targets = self[1:]
        self.comments = comments
        # list of positions corresponding to source positions, one for each target word
        self.alignments = [self.init_alignment(l) for l in range(len(self)-1)]
        # set the number of phone matches for each source word position; to be updated
        # after each re-alignment
        self.set_matches()

    def __repr__(self):
        """Print name for the TGroup."""
        return "<<{}...>>".format(self.source.name)

    def get_source(self, index):
        return self[index]

    def get_targets(self, index):
        return [w for (i, w) in enumerate(self) if i != index]

    def copy_alignments(self):
        """Make a copy of the current alignments so it doesn't change as the algorithm
        is tried in one or the other direction."""
        alignments = []
        for alignment in self.alignments:
            alignments.append(alignment[:])
        return alignments

    def split(self, index, assign=False):
        """Separate the TGroup into source and target Words."""
        source = None
        targets = []
        for i, w in enumerate(self):
            if i == index:
                source = w
            else:
                targets.append(w)
        if assign:
            self.source = source
            self.targets = targets
        return source, targets

    def print_alignments(self, alignments=None, verbosity=0):
        """Pretty print alignments. If alignments isn't passed, the currently
        saved alignments are printed."""
        alignments = alignments or self.alignments
#        print("source {}".format(self.source))
        current_positions = [-1] * len(self.targets)
        widths = [0] * len(self.source)
        t_strings = [[] for x in range(len(self.targets))]
        s_string = []
        for si, sphone in enumerate(self.source):
            if verbosity:
                print()
                print(" source phone {}:{}".format(si, sphone))
            current_max = len(sphone) + 1
            string = " " + sphone
            strings = []
            for ti, (alignment, target) in enumerate(zip(alignments, self.targets)):
                align = alignment[si]
                current_tpos = current_positions[ti]
                if verbosity:
                    print()
                    print("  current max {}".format(current_max))
                    print("  target {}, alignment {}".format(target, align))
                    print("  current tpos {}".format(current_tpos))
                if align == -1:
                    # sphone corresponds to nothing in target
#                    print("  not aligned")
                    strings.append('  ')
                elif align > current_tpos + 1:
#                    print("  gap since last alignment")
                    # positions in target since last aligned phone
                    new_phones = target[current_tpos+1:align+1]
                    new_length = sum([len(p)+1 for p in new_phones])
#                    print("  new phones {}, new length {}".format(new_phones, new_length))
                    current_max = max(current_max, new_length)
                    strings.append(''.join([' ' + p for p in new_phones]))
                    current_positions[ti] = align
                else:
#                    print("  alignment agrees with current position")
                    new_phone = target[align]
                    new_length = len(new_phone)+1
#                    print("  new phone {}, new length {}".format(new_phone, new_length))
                    current_max = max(current_max, new_length)
                    if si > 0 and align == alignment[si-1]:
                        # This character already part of string
                        strings.append(' ' + ' ' * len(new_phone))
                    else:
                        strings.append(' ' + new_phone)
                    current_positions[ti] = align
                if verbosity:
                    print("  updated max {}".format(current_max))
            if verbosity:
                print(" final max width for position {}: {}".format(si, current_max))
            widths[si] = current_max
            string = string.rjust(current_max)
            s_string.append(string)
            for i, s in enumerate(strings):
                strings[i] = s.rjust(current_max)
            for i, s in enumerate(strings):
                t_strings[i].append(s)
            if verbosity:
                print(" target strings {}".format(strings))
                print(" source string {}".format(string))
        # Targets that have positions beyond position aligned with end of source
        for ti, (target, alignment, t_string) in enumerate(zip(self.targets, self.alignments, t_strings)):
            # last aligned position that is not -1
            last_position = [p for p in alignment if p >= 0][-1]
            if last_position < len(target) -1:
                t_string.append(''.join([' ' + s for s in target[last_position+1:]]))
        s_string = ''.join(s_string)
        for i, s in enumerate(t_strings):
            t_strings[i] = ''.join(t_strings[i])
        if verbosity:
            print("STRINGS")
        print(s_string)
        for s in t_strings:
            print(s)
        return widths

    def get_spans(self, source, target, alignment):
        """Given an alignment of source with target, return the spans."""
        si, ti = 0, 0
        spans = []
        current_span = []
        print("SOURCE {}, len {}; TARGET {}, len {}".format(source, len(source), target, len(target)))
        while si < len(source) or ti < len(target):
#            print("SI {}, TI {}".format(si, ti))
            current = None
            if ti >= len(target):
                # DELETE
                s = source[si]
#                print("Deleting {} at end of source".format(s))
                current = (si, ti, s, 'd')
                si += 1
            elif si >= len(source):
                # INSERT
                t = target[ti]
#                print("Inserting {} at end of source".format(t))
                current = (si, ti, t, 'i')
                ti += 1
            else:
                a = alignment[si]
                if a < 0:
                    # DELETE
                    s = source[si]
#                    print("Deleting {}".format(s))
                    current = (si, ti, s, 'd')
                    si += 1
                elif ti < a:
                    # INSERT
                    t = target[ti]
#                    print("Inserting {}".format(t))
                    current = (si, ti, t, 'i')
                    ti += 1
                else:
                    # either substitute or nochange
                    s = source[si]
                    t = target[ti]
                    if s == t:
#                        print("No change: {}".format(s))
                        current = (si, ti, s, 'n')
                    else:
#                        print("Substituting {} for {}".format(t, s))
                        current = (si, ti, (s, t), 's')
                    si += 1
                    ti += 1
            # Decide whether to include with current span or close off
#            print("current span {}".format(current_span))
#            print("new {}".format(current))
#            print("spans {}".format(spans))
            if current[-1] == 'n':
                if current_span and current_span[-1][-1] == 'n':
                    # continue with nochange span
                    current_span.append(current)
                else:
                    # start new nochange span
                    spans.append(current_span)
                    current_span = [current]
            elif current_span and current_span[-1][-1] != 'n':
                # continue change span
                current_span.append(current)
            else:
                # start new change span
                if current_span:
                    spans.append(current_span)
                current_span = [current]
        if current_span:
            spans.append(current_span)
        return spans

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

    def align(self, direction=None, alignments=None, verbosity=0):
        """Repeatedly align source to target words, minimizing edit costs,
        until there are no changes to alignments."""
        alignments = alignments or self.alignments
        change = True
        cost = 0
        if verbosity:
            print()
            print("DIRECTION: {}, INITIAL ALIGNMENTS:".format(direction))
            self.print_alignments(alignments)
        iteration = 1
        while change:
            if verbosity:
                print()
                print("ALIGNING SOURCE TO TARGETS in {}, DIRECTION {}, ITERATION {}".format(self, direction, iteration))
            cost, change = self.minimize_all(direction=direction, alignments=alignments,
                                             verbosity=verbosity)
            iteration += 1
        if verbosity:
            print()
            print("DIRECTION: {}, FINAL ALIGNMENTS:".format(direction))
            self.print_alignments(alignments)
        return cost

    def minimize_all(self, direction=None, alignments=None, verbosity=0):
        """Minimize costs once for all target words, adjusting alignments accordingly."""
        total_cost = 0
        change = False
        for tword_i in range(len(self.targets)):
            # make 3 copies of the current alignment for this target word
            alignment = alignments[tword_i]
            current_alignment = alignment[:]
            cost = self.minimize(tword_i, forward=direction=='forward', alignment=alignment,
                                 verbosity=verbosity)
            # Update the matches vector.
            self.set_matches()
            if verbosity:
                print(" Cost for target {}: {}".format(tword_i, cost))
            total_cost += cost
            if self.alignments[tword_i] != current_alignment:
                change = True
#        # Update the matches vector.
#        self.set_matches()
        if change:
            print("COST: {}; SOME ALIGNMENT CHANGED, REALIGNING".format(total_cost))
        else:
            print("COST: {}; CONVERGED, NO ALIGNMENT CHANGED".format(total_cost))
        return total_cost, change

    def minimize(self, tword_index, forward=False, alignment=None,
                 verbosity=0):
        """Minimize edit costs for the source word and target word at tword_index,
        returning the cost and changing the alignment if necessary."""
        source = self.source
        target = self.targets[tword_index]
        # if alignment is specified it's a copy of the current alignment
        alignment = alignment or self.alignments[tword_index]
        sourcecopy = source.copy()
        if verbosity:
            print(" Aligning {} with {}".format(source, target))
        if forward:
            sstart = 0; tstart = 0
        else:
            sstart = len(source)-1; tstart = len(target)-1
        return self.minimize1(source, target, sstart, tstart, alignment,
                              forward=forward, sourcecopy=sourcecopy, verbosity=verbosity)

    def minimize1(self, source, target, sindex, tindex, alignment,
                  forward=False, sourcecopy=None, copyoffset=0, verbosity=0):
        """Return the cost of editing source word to produce target word,
        ending in positions sindex in source, tindex in target. Change the current
        alignment where there are deletions or substitutions.
        Start at the *beginning* or *end* of source and target words, recursively figuring the
        cost, using the expression on p. 1187 of Durrett & DeNero.
        forward option determines direction of recursion."""
        if forward:
            sout = sindex >= len(source)
            tout = tindex >= len(target)
        else:
            sout = sindex < 0
            tout = tindex < 0
        if sout and tout:
            return 0
#        if sindex >= len(source) and tindex >= len(target):
#            return 0
        if forward:
            copyindex = sindex + copyoffset
        else:
            copyindex = sindex
        if verbosity:
            print("Checking positions {}:{}, current phones: {}:{}".format(sindex, tindex,
                                                                           source[sindex] if not sout else '',
                                                                           target[tindex] if not tout else ''))
        if not sourcecopy:
            sourcecopy = source.copy()
        costs = self.costs(source, target, sindex, tindex, forward=forward, verbosity=verbosity)
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
            sourcecopy.delete(copyindex)
            copyoffset += -1
        else:
            phone = target[tindex]
            if oper == 0:
                # Insertion; this does *not* change the alignment, it seems
                if verbosity > 1:
                    print("  Inserting {}".format(phone))
                # for forward algorithm, we insert before the current index
                if forward:
                    sourcecopy.insert(copyindex-1, phone)
                else:
                    sourcecopy.insert(copyindex, phone)
                copyoffset += 1
            else:
                # Substitution
                alignment[sindex] = tindex
                sphone = source[sindex]
                if sphone == phone:
                    if verbosity > 1:
                        print("  No change for {}".format(phone))
                else:
                    if verbosity > 1:
                        print("  Substituting {} for {} in position {}".format(phone, sphone, sindex))
                    sourcecopy.substitute(copyindex, phone)
        if verbosity:
            print(" Cost: {}, si* {}, ti* {}".format(cost, newsindex, newtindex))
            print(" Current transformed word: {}".format(sourcecopy.chars()))
            print(" Current alignment: {}".format(alignment))
        return cost + self.minimize1(source, target, newsindex, newtindex,
                                     alignment, sourcecopy=sourcecopy,
                                     forward=forward, copyoffset=copyoffset, verbosity=verbosity)

    def insert_cost(self, sindex=-1, ldiff=0, tindex=-1, tphone='', scontext=None, tcontext=None, forward=False):
        """
        Return the cost of inserting tphone in position sindex in source,
        the new indices to check next,
        and 0 (id for insertion).
        """
        if forward:
            indices = (sindex, tindex+1)
        else:
            indices = (sindex, tindex-1)
        return (self.problem.score_insert(sindex, tphone, ldiff, scontext, forward=forward),
                indices, 0)

    def delete_cost(self, sindex=-1, sphone='', ldiff=0, tindex=-1, scontext=None, forward=False):
        """
        Return the cost of deleting sphone in position sindex in source,
        the new indices to check next,
        and 1 (id for deletion).
        """
        if forward:
            indices = (sindex+1, tindex)
        else:
            indices = (sindex-1, tindex)
        return (self.problem.score_delete(sindex, sphone, ldiff), indices, 1)

    def compare_cost(self, sindex=-1, sphone='', tindex=-1, tphone='',
                     ldiff=0, scontext=None, tcontext=None, matches=0, forward=False):
        """
        Return the cost of keeping sphone (if it's the same as tphone) or replacing it with tphone (if they're different),
        the new indices to check next,
        and 2 (id for substitution).
        """
        if forward:
            indices = (sindex+1, tindex+1)
        else:
            indices = (sindex-1, tindex-1)
        return (self.problem.score_compare(sindex, tphone, sphone, scontext, tcontext, matches, forward=forward),
                indices, 2)

    def costs(self, source, target, sindex, tindex, forward=False, verbosity=0):
        """
        Return the costs associated with each operation, along with updated source and target indices, and indices
        to distinguish the operations.
        """
        # difference in length of portions of source and target words remaining
        if forward:
            ldiff = len(source) - len(target) - sindex + tindex
        else:
            ldiff = sindex - tindex
        # current source and target phones in positions sindex and tindex
        if forward:
            sphone = source[sindex] if sindex < len(source) else ''
            tphone = target[tindex] if tindex < len(target) else ''
        else:
            sphone = source[sindex] if sindex >= 0 else ''
            tphone = target[tindex] if tindex >= 0 else ''
        # preceding contexts of in source and target words
        if forward:
            scontext = source[sindex+1:sindex+1+TGroup.context_size] if sindex+1 < len(source) else []
            tcontext = target[tindex+1:tindex+1+TGroup.context_size] if tindex+1 < len(target) else []
        else:
            scontext = source[max([0,sindex-TGroup.context_size]):sindex] if sindex >= 1 else []
            tcontext = target[max([0,tindex-TGroup.context_size]):tindex] if tindex >= 1 else []
        # number of other target words with the current source phone in positions aligned with sindex
        if forward:
            matches = self.matches[sindex] if sindex < len(source) else 0
        else:
            matches = self.matches[sindex] if sindex >= 0 else 0
        if sindex < 0 or sindex >= len(source):
            # insertion is the only possibility because we're reached the beginning of the source word
            return [self.insert_cost(sindex=sindex, ldiff=ldiff, tindex=tindex, tphone=tphone, scontext=scontext,
                                     tcontext=tcontext, forward=forward)]
        if tindex < 0 or tindex >= len(target):
            # deletion is the only possibility because we've reached the beginning of the target word
            return [self.delete_cost(sindex=sindex, sphone=sphone, ldiff=ldiff, tindex=tindex, scontext=scontext, forward=forward)]
        if verbosity > 1:
            print("  Calculating costs for {}:{}; {}:{}".format(sindex, tindex, sphone, tphone))
        # all three edit operations are possible in these source and target positions
        return [self.insert_cost(sindex=sindex, ldiff=ldiff, tindex=tindex, tphone=tphone,
                                 scontext=scontext, tcontext=tcontext, forward=forward),
                self.delete_cost(sindex=sindex, sphone=sphone, ldiff=ldiff, tindex=tindex,
                                 scontext=scontext, forward=forward),
                self.compare_cost(sindex=sindex, sphone=sphone, tindex=tindex, tphone=tphone,
                                  ldiff=ldiff, scontext=scontext, tcontext=tcontext, matches=matches, forward=forward)]

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
        
