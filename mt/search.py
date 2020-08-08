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

-- 13.7.2020
   Created.
"""

#import queue
import heapq

class Searcher:
    """
    Simple unconstrained search using the queue module.
    """

    cutoff = 2000
    qmax = 100

    def __init__(self, name='', goal_test=None,
                 extend=None, make_start=None):
        self.name = name
        # Function taking a state testing whether it's a goal state
        self.goal_test = goal_test
        # Function tatking a state returning a list of new new states
        # reached from it
        self.extend = extend
        # Function of no arguments creating the start state
        self.make_start = make_start
        # Queue of states
        self.queue = self.make_queue()

    def make_queue(self):
        """Create the queue."""
#        return queue.SimpleQueue()
        return []

    def empty_queue(self):
        """Is the queue empty?"""
        return not self.queue

    def add_state(self, state):
        """Add a new state to the queue."""
        heapq.heappush(self.queue, state)
#        self.queue.put(state)

    def add_states(self, states):
        """
        Add new states (already sorted) to the queue.
        """
        for state in states:
            self.add_state(state)

    def truncate(self, qmax=0):
        """
        Shorten the queue to qmax
        """
        qmax = qmax or Searcher.qmax
        if self.qsize() > qmax:
            self.queue[qmax:] = []

    def qsize(self):
        """Number of states in the queue."""
#        return self.queue.qsize()
        return len(self.queue)

    def run(self, cutoff=0, truncate=True, verbosity=0):
        """
        Run search, stopping after cutoff iterations if no goal
        is found.
        """
        cutoff = cutoff or Searcher.cutoff
        # Create start state
        self.add_state(self.make_start())
        iter = 0
        success = False
        while (not success and iter < cutoff and not self.empty_queue()):
            if verbosity:
                print("Iteration {}, q length {}".format(iter, self.qsize()))
            success = self.expand(verbosity=verbosity)
            if truncate:
                self.truncate()
            iter += 1
        # Final (goal) state
        if not success:
            print("{} failed!".format(self))
        return success

    def expand(self, verbosity=0):
        """
        Expand the queue by taking a state from it and replacing
        it with new states reached by extending that state.
        """
        if self.empty_queue():
            if verbosity:
                print("No more states")
            return False
#        next_state = self.queue.get()
        next_state = heapq.heappop(self.queue)
        if verbosity:
            print(" Next {}".format(next_state))
        if self.goal_test(next_state):
            if verbosity:
                print("goal state")
            return next_state
        new_states = self.extend(next_state)
        for ns in new_states:
            self.add_state(ns)
        return False

class BestFirst(Searcher):
    """
    Best first search using queue.PriorityQueue to rank
    states in the queue by their values.
    """

    def __init__(self, name='', goal_test=None, make_start=None,
                 extend=None, evaluate=None):
        Searcher.__init__(self, name=name, goal_test=goal_test,
                          make_start=make_start, extend=extend)
        self.evaluate=evaluate

    def __repr__(self):
        return "BFS:{}".format(self.name)

    def make_queue(self):
        """Create a PriorityQueue."""
        return []
#        return queue.PriorityQueue()

    def add_state(self, state):
        """
        Add a state to the queue, along with its
        value (cost+distance to goal); low values go to the front.
        """
#        self.queue.put((self.evaluate(state), state))
        heapq.heappush(self.queue, (self.evaluate(state), state))

    def expand(self, verbosity=0):
        """
        Expand the queue by removing the first (current best)
        state and replacing it in the queue by the states
        reached by extending it.
        """
        if self.empty_queue():
            if verbosity:
                print("No more states")
            return False
        # This line is different from Searcher
#        next_value, next_state = self.queue.get()
        next_value, next_state = heapq.heappop(self.queue)
        if verbosity:
            print("{}\nVALUE {}".format(next_state, next_value))
        if self.goal_test(next_state):
            if verbosity:
                print(" Goal state")
            return next_state, next_value
        new_states = self.extend(next_state)
        for ns in new_states:
            self.add_state(ns)
        return False

# # An example to test the algorithm
# MAP = {'S': {'a': 1.5, 'd': 2},
#        'a': {'b': 2},
#        'd': {'e': 3},
#        'b': {'c': 3},
#        'e': {'G': 2},
#        'c': {'G': 4}}
#
# CROW = {'a': 4, 'b': 2, 'c': 4, 'd': 4.5, 'e': 2, 'S': 5.5, 'G':0}
#
# def dests_from(origin):
#     if origin not in MAP:
#         return {}
#     return MAP[origin]
#
# def dist_to(distances, dest):
#     if dest not in distances:
#         return -1
#     return distances[dest]
#
# def extend(val_place):
#     value, place = val_place
#     print("Extending {}, {}".format(value, place))
#     if place not in MAP:
#         return []
#     dests = dests_from(place)
#     dist_dests = [[v, k] for k, v in dests.items()]
#     dist_dests = [[d + CROW[s] + value, s] for d, s in dist_dests]
#     print("New states {}".format(dist_dests))
#     return dist_dests
#
# S = Searcher('map', goal_test=lambda s: s[1] == 'G',
#              make_start=lambda: [0, 'S'],
#              extend=lambda s: extend(s))
#
# BF = BestFirst('mapBF',
#                goal_test=lambda s: s[1] == 'G',
#                make_start=lambda: [0, 'S'],
#                extend=lambda s: extend(s),
#                evaluate=lambda s: s[0])
