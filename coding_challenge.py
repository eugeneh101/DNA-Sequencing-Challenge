from collections import defaultdict
from difflib import SequenceMatcher


class GeneSequencer(object):


    def __init__(self):
        self.sequence_name_string = dict()
        self.match_to_left = defaultdict(set)
        self.match_to_right = defaultdict(set)
        self.match_to_left_confirmed = defaultdict(dict)
        self.match_to_right_confirmed = defaultdict(dict)
        self.final_frag_name_sequence = None
        self.correct_string_sequence = None

    def load_data(self, file_name):
        """Loads fragment name and fragment string into dictionary from a file name.
        
        Given a file name, will load the data and perform simple cleaning to 
        load into sequene_name_string dictionary attribute.
        """
        f = open(file_name, 'rb')
        lines = [line[:-1] for line in f.readlines()] # remove '\n'

        starter_indices = []
        for i, line in enumerate(lines): # assumes that a starter will not be a last line
            if line[0] == '>':
                starter_indices.append(i)
        starter_indices.append(len(lines)) # adds the terminating length for the zip below

        for start_index, end_index in zip(starter_indices[0:-1], starter_indices[1:]):
            self.sequence_name_string[lines[start_index][1:]] = ''.join(lines[
                    start_index + 1:end_index])

    def find_adjacent_pairs(self):
        """For any fragment string, find POSSIBLE fragments that can connect
        to the left or right
        
        Given a fragment name, find its respective DNA string and find what
        fragment strings that can be paired to the right and left of its string.
        Save possible right and left fragment names into dictionary. This method
        doesn't guarantee that the pairing is possible but runs very quickly.
        """
        for frag_name, frag_string in self.sequence_name_string.iteritems():
            frag_string_first_half = frag_string[ :len(frag_string) / 2 + 1]
            frag_string_second_half = frag_string[len(frag_string) / 2: ]

            for other_frag_name, other_frag_string in \
                self.sequence_name_string.iteritems():
                if frag_name != other_frag_name:
                    if frag_string_first_half in other_frag_string:
                        self.match_to_left[frag_name].add(other_frag_name)
                        self.match_to_right[other_frag_name].add(frag_name)
                    if frag_string_second_half in other_frag_string:
                        self.match_to_right[frag_name].add(other_frag_name)
                        self.match_to_left[other_frag_name].add(frag_name)
                                                
    def find_adjacent_confirmed_pairs(self):
        """For any fragment string, find GUARANTEED fragments that can connect
        to the left or right
        
        Builds on top of find_adjacent_pairs() method where it will find
        guaranteed fragments that can connect to the left or right. 
        find_adjacent_pairs() uses a simple heuristic that is quick and filters
        out possible candidates. This method confirms whether the fragment can
        be connected by checking if the overlapping characters are at the 
        beginning or end of the possible paired fragment string. Although there 
        are 3 'for' loops, each loop is iterating over a fairly small list.
        """
        
        def __confirmed_match_pairs(f_n, f_s, o_f_n, o_f_s, matches):
            """Adds confirmed adjacent string name into class dictionary
            f_n: fragment name
            f_s: fragment string
            o_f_n: other fragment name
            o_f_s: other fragment string
            matches: output object from SequenceMatcher
            """
            for match in [temp for temp in matches if (temp.size > 
                                                       len(f_s) / 2)]:
                if match.a == 0 and match.b + match.size == len(o_f_s):
                    #print "match to the left"
                    self.match_to_left_confirmed[f_n][o_f_n] = match.size
                    self.match_to_right_confirmed[o_f_n][f_n] = match.size
                if match.b == 0 and match.a + match.size == len(f_s):
                    #print"match to the right"
                    self.match_to_right_confirmed[f_n][o_f_n] = match.size
                    self.match_to_left_confirmed[o_f_n][f_n] = match.size

        for frag_name, other_frag_name_list in self.match_to_left.items() + \
            self.match_to_right.items():
            frag_string = self.sequence_name_string[frag_name]
            for other_frag_name in other_frag_name_list:
                other_frag_string = self.sequence_name_string[other_frag_name]
                matches = SequenceMatcher(None, frag_string, other_frag_string, 
                                          autojunk=False).get_matching_blocks()
                __confirmed_match_pairs(frag_name, frag_string, other_frag_name, 
                                    other_frag_string, matches)

    def build_final_string_from_left_side(self):
        """Build a list of fragment names where each adjacent fragment name
        is the correct fragment that connects to the left or right in the 
        final solution.
        
        Assumes there exists a fragment name such that it has no fragment that
        can connect the left. If such fragment exists, then it must be the
        first fragment name on the left. Then, connect fragment names to the
        right of this fragment name until all fragment names are used. Saves
        this list of fragment names.
        """

        left_starter = set(self.sequence_name_string.keys()) - \
            set(self.match_to_left_confirmed.keys())
        #print left_starter
        if len(left_starter) != 1:
            print """Either there exists no fragment name that is the first 
            fragment OR there are more than 1 fragment name. 0 fragment names
            occur when the leftmost fragment (call it frag0) CAN be connected 
            to another fragment on frag0's left side but is NOT connected in
            the final unique solution. If there exist more than 2 fragment 
            names, then there exists a problem in either the input file or 
            previous steps.
            """
            return
        # if the if statement is false, there should be a guaranteed solution
        left_starter = left_starter.pop()
        sequence_dict = {None: [[None]], left_starter: [[]]}
        current_keys = [left_starter]
        #print current_keys

        for _ in xrange(len(self.sequence_name_string) - 1): # has finite iterations
            next_keys = []
            for current_key in current_keys:
                for right_pair in self.match_to_right_confirmed[current_key]:
                    temp_sequence = [element + [current_key] for element in 
                                     sequence_dict[current_key]]
                    sequence_dict[right_pair] = temp_sequence
                    next_keys.append(right_pair)
            current_keys = next_keys

        for frag_name, frag_seq_list in sequence_dict.iteritems():
            for frag_seq in frag_seq_list:
                if len(set(frag_seq)) + 1 == len(self.sequence_name_string) \
                    and len(frag_seq) + 1 == len(self.sequence_name_string):
                    self.final_frag_name_sequence = frag_seq + [frag_name]
                    #print frag_seq
                    return

    def build_final_string_from_right_side(self):
        """Does the same thing as build_final_string_from_left_side except
        the starter fragment is on the right side of the final string. This
        method is ONLY called if build_final_string_from_left_side failed to
        find correct solution.

        This method is used when the leftmost fragment in the final solution
        (call it frag0) CAN be connected to another fragment on frag0's 
        left side but is NOT connected in the final unique solution. Then,
        I cannot used build_final_string_from_left_side. Hence, I try to
        build from right to left instead.
        
        Most of the code is copied and pasted from the previous method.
        However, there are modifications to build the fragment sequence from
        the right side to the left side. 
        """
        
        if self.final_frag_name_sequence:
            return # only excute if build_final_string_from_left_side failed

        right_starter = set(self.sequence_name_string.keys()) - \
            set(self.match_to_right_confirmed.keys())
        #print right_starter
        if len(right_starter) != 1:
            print """Either there exists no fragment name that is the last 
            fragment OR there are more than 1 fragment name. 0 fragment names
            occur when the rightmost fragment (call it frag_last) CAN be 
            connected to another fragment on frag_last's right side but is NOT 
            connected in the final unique solution. If there exist more than 2 
            fragment names, then there exists a problem in either the input 
            file or previous steps.
            """
            return
        # if the if statement is not true, there should be a guaranteed solution
        right_starter = right_starter.pop()
        sequence_dict = {None: [[None]], right_starter: [[]]}
        current_keys = [right_starter]
        #print current_keys
        
        # has finite iterations even if code breaks
        for _ in xrange(len(self.sequence_name_string) - 1):
            next_keys = []
            for current_key in current_keys:
                for left_pair in self.match_to_left_confirmed[current_key]:
                    temp_sequence = [element + [current_key] for element in 
                                     sequence_dict[current_key]]
                    sequence_dict[left_pair] = temp_sequence
                    next_keys.append(left_pair)
            current_keys = next_keys

        for frag_name, frag_seq_list in sequence_dict.iteritems():
            for frag_seq in frag_seq_list:
                if len(set(frag_seq)) + 1 == len(self.sequence_name_string) \
                    and len(frag_seq) + 1 == len(self.sequence_name_string):
                    self.final_frag_name_sequence = (frag_seq + 
                         [frag_name])[::-1] # have to reverse the sequence
                    #print frag_seq
                    return

    def build_final_string_from_any_fragment(self):
        """Used as last resort to find correct sequence of fragment names.
        
        If both build_final_string_from_left_side and 
        build_final_string_from_right_side have failed (ie. no identifiable
        fragment to start building sequence from left- or right-most side), 
        this method can find final sequence solution from ANY node.
        It will guess what is the first fragment (call it frag_guess) 
        is and build from left to right. If frag_guess is not the correct
        starting fragment for the left side, then it will try another 
        fragment as the first fragment.
        This method is much slower than in terms of iterations compared to
        the other previous methods but nonetheless runtime is still quick. 
        Basically, it will force a solution out at the end. Actually, a 
        modified version of build_final_string_from_left_side as it builds
        from left to right.
        """
        
        if self.final_frag_name_sequence:
            return # only excute if both previous methods have failed

        for starter in self.sequence_name_string:
            #print starter
            sequence_dict = {None: [[None]], starter: [[]]}
            current_keys = [starter]
            # finite number of iterations even if code breaks
            for _ in xrange(len(self.sequence_name_string) - 1):
                next_keys = []
                for current_key in current_keys:
                    for right_pair in self.match_to_right_confirmed[current_key]:
                        temp_sequence = [element + [current_key] for element in 
                                         sequence_dict[current_key]]
                        sequence_dict[right_pair] = temp_sequence
                        next_keys.append(right_pair)
                current_keys = next_keys

            for frag_name, frag_seq_list in sequence_dict.iteritems():
                for frag_seq in frag_seq_list:
                    if len(set(frag_seq)) + 1 == len(self.sequence_name_string)\
                        and len(frag_seq) + 1 == len(self.sequence_name_string):
                        self.final_frag_name_sequence = frag_seq + [frag_name]
                        #print frag_seq
                        return # when solution found, return from method

    def return_unique_string(self):
        """Extracts the final string output by chaining the adjacent
        fragment strings together while accounting for overlaps.

        This method is called after build_final_string_from_left_side
        OR build_final_string_from_right_side has successfully created
        the correct sequence. Nonetheless, the final string is constructed
        from left to right.
        """

        sequence_first = [None] + self.final_frag_name_sequence
        sequence_second = self.final_frag_name_sequence + [None]
        for current_frag, next_frag in zip(sequence_first, sequence_second):
            if current_frag is None:
                self.correct_string_sequence = self.sequence_name_string[
                    next_frag]
                continue
            elif next_frag is None:
                return self.correct_string_sequence
            num_overlapping_chars = self.match_to_right_confirmed[
                current_frag][next_frag]
            self.correct_string_sequence += self.sequence_name_string[
                next_frag][num_overlapping_chars:]


def main(file_name):
    gene_sequencer = GeneSequencer()
    gene_sequencer.load_data(file_name)
    gene_sequencer.find_adjacent_pairs()
    gene_sequencer.find_adjacent_confirmed_pairs()
    gene_sequencer.build_final_string_from_any_fragment()
    return gene_sequencer.return_unique_string()