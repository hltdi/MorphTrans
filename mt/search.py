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

import queue

# class State:
#
#     def __init__(self, obj=0, parent=None, cost=0):
#         self.parent = parent
#         self.cost = cost
#         self.obj = obj
#
#     def __repr__(self):
#         return "{}:{}".format(self.position, self.cost)

class Searcher:
    """
    Simple unconstrained search using the queue module.
    """

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
        return queue.SimpleQueue()

    def empty_queue(self):
        """Is the queue empty?"""
        return not self.queue

    def add_state(self, state):
        """Add a new state to the queue."""
        self.queue.put(state)

    def qsize(self):
        """Number of states in the queue."""
        return self.queue.qsize()

    def run(self, cutoff=50):
        """
        Run search, stopping after cutoff iterations if no goal
        is found.
        """
        # Create start state
        self.add_state(self.make_start())
        iter = 0
        success = False
        while (not success and iter < 100 and not self.empty_queue()):
            print("Iteration {}".format(iter))
            success = self.expand()
            iter += 1

    def expand(self):
        """
        Expand the queue by taking a state from it and replacing
        it with new states reached by extending that state.
        """
        if self.empty_queue():
            print("No more states")
            return False
        next_state = self.queue.get()
        print(" Next {}".format(next_state))
        if self.goal_test(next_state):
            print("goal state")
            return True
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
        return queue.PriorityQueue()

    def add_state(self, state):
        """
        Add a state to the queue, along with its
        value (cost+distance to goal); low values go to the front.
        """
        self.queue.put((self.evaluate(state), state))

    def expand(self, verbosity=1):
        """
        Expand the queue by removing the first (current best)
        state and replacing it in the queue by the states
        reached by extending it.
        """
        if self.empty_queue():
            if verbosity:
                print("No more states")
            return False
        if verbosity:
            print("Queue length {}".format(self.qsize()))
        # This line is different from Searcher
        next_value, next_state = self.queue.get()
        if verbosity:
            print(" Next state\n{}\n VALUE {}".format(next_state, next_value))
        if self.goal_test(next_state):
            if verbosity:
                print(" Goal state")
            return True
        new_states = self.extend(next_state)
        for ns in new_states:
#            print(" Adding new state\n{}".format(ns))
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
