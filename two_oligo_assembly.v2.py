# This script needs codontable.tab in the same directory
#
# Authors: Eva Maria Strauch, Benjamin Basanta 2018
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.

import os
import re
from random import randint
from random import shuffle
import sys
from Bio import SeqIO
import optparse
from Bio.Alphabet.IUPAC import unambiguous_dna, ambiguous_dna
from Bio.SeqUtils import  MeltingTemp
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from multiprocessing import Pool
import argparse

########## Option system: ###########
argparser = argparse.ArgumentParser(description='Split genes in orthogonal pieces that can be used in multiplex assembly')
argparser.add_argument('-input_list', type=str,help='Name of file containing: mygenename DNAsequence AAsequence')
argparser.add_argument('-adaptor_fname', type=str,default='./pool_adaptors_short_list.txt',help='Name of file containing adaptor sequences. Format: First line: column names, followed by lines: adaptor_name fiveprime_5 fiveprime_3 threeprime_5 threeprime_3')
argparser.add_argument('-codontable_fname', type=str,default='./codontable.tab',help='Codon table to use')
argparser.add_argument('-adaptor_number', type=int,help='What adaptor to use? starting at 1')
argparser.add_argument('-nproc', type=int,default=1,help='Number of processors to use. Must be in the same node (DIGs: -cX -N1 AFAIK)')
argparser.add_argument('-dNTPsmM', type=float, default=0.3, help='Concentration of dNTPs in mM')
argparser.add_argument('-MgmM', type=float, default=2.5, help='Concentration of Mg++ in mM')
argparser.add_argument('-min_melt_temp', type=int,default=62,help='minimum melting temperature of the overlapping region')
argparser.add_argument('-maxoverlaplen', type=int,default=30,help='What is the max overlap length?')
argparser.add_argument('-max_oligo_size', type=int,default=230,help='Absolute max length of orderable oligo')
argparser.add_argument('-alignthreshold', default=15, type=int,help='Number of positions over which an alignment will be considered as such')
args = argparser.parse_args()

####################################
minmelt = args.min_melt_temp
# Maximum number of consecutively paired bases when looking for mismatched priming:
alignthreshold = args.alignthreshold
#final oligos length
#max_oligo_size = 160
max_oligo_size = args.max_oligo_size
#
yeast_primer5 = "GGGTCGGCTTCGCATATG"
yeast_primer3 = "CTCGAGGGTGGAGGTTCC"

#files to write:

# 1. sequences
final_sequences = open('final_sequences_pool_%d.tab'%args.adaptor_number,'w')

# 2. filtered out sequences
trashed_sequences = open('ChimeraForming_pool_%d.tab'%args.adaptor_number,'w')
impossible_sequences = open('Impossible_pool_%d.tab'%args.adaptor_number,'w')
# 4. sequences with adapters ready for order
final_order_large = open('final_order_large_pool_%d.tab'%args.adaptor_number, 'w')

def reverse_complement(seq):
    rc = Seq(seq, ambiguous_dna ).reverse_complement()
    return rc

def read_adapters(infile):
    handle = open(infile, 'r')
    # get all lines, skip final newline and drop first line that has titles:
    lines = [ line[:-1] for line in handle.readlines() ][1:]
    print("Your input adaptors:")
    for line in lines:
        print(line)
    fiveprime_5 = [ line.split()[1] for line in lines ]
    fiveprime_3 = [ line.split()[2] for line in lines ]
    threeprime_5 = [ line.split()[3] for line in lines ]
    threeprime_3 = [ line.split()[4] for line in lines ]
    return fiveprime_5, fiveprime_3, threeprime_5, threeprime_3

def number_to_base(num):
    if num is 1:
        return 'A'
    elif num is 2:
        return 'C'
    elif num is 3:
        return 'T'
    elif num is 4:
        return 'G'
    else:
        print "Input number is not 1, 2, 3 or 4, review code"
        sys.exit(0)

def get_random_bases(length):
    s = ""
    for i in range(length):
        s += number_to_base(randint(1,4))
    return s


def isComplement(x, y):
    cs = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}
    return cs[x] == y


def complement_positions(top, bottom):
    positions = []
    for i, (t, b) in enumerate(zip(top, bottom)):
        if isComplement(t, b):
            positions.append(i)
    return positions

def get_score(positions):
    # for this function see Nguyen-Dumon, BioTechniques, Vol 55, No 2, August 2013, pp 69
    '''
    Score is the number of consecutive positions that paired in the input
    alignment positions.
    '''
    score = 0
    maxscore = 0
    if len(positions) == 1:
        score = 1
    elif len(positions) > 1:
        score = 1
        lastposition = positions[0]
        for position in positions[1:]:
            if position == lastposition + 1:
                score += 1
            else:
                if score > maxscore:
                    maxscore = score
                score = 1
            lastposition = position
    if score > maxscore:
        maxscore = score
    return maxscore


class Dimer(object):
    # for this class see Nguyen-Dumon, BioTechniques, Vol. 55, No. 2, August 2013, pp. 69
    def __init__(self, top, bottom=None):
        self._top = top
        # If bottom not available: do it against top itself
        self._bottom = bottom if bottom is not None else top[::-1]
        self._maxscore = None  # (score, direction, startIndex,  positions)

    def score(self, align_threshold):

        def detect_dimers():
            matches = []
            if (len(self._top) < align_threshold) or (len(self._bottom) < align_threshold):
                # No need to check for continous matches in a dimer where one of the
                # sequences doesn't even reach align_threshold in length.
                print("The priming region is too short to mismatch > align_threshold")
                return matches
            #for i in range(0, len(self._top)):
            for i in range(0, len(self._top)-align_threshold):
                top = self._top[i:]
                positions = complement_positions(top, self._bottom)
                matches.append(('forward', i, positions))
            for i in range(1, len(self._bottom)):
                bottom = self._bottom[i:]
                positions = complement_positions(self._top, bottom)
                matches.append(('backward', i, positions))
            return matches

        def find_maxscore(matches):
            maxscore = (0, None, None, None)
            for direction, index, positions in matches:
                score = get_score(positions)
                if score > maxscore[0]:
                    maxscore = (score, direction, index, positions)
            return maxscore

        def print_dimer(direction, index, positions):
            if positions:
                matchSigns = [s for s in ' '*min(len(self._top),
                                                 len(self._bottom))]
                for i in positions:
                    matchSigns[i] = '|'
                matchSigns = ''.join(matchSigns)
                indentation = index
                if direction == 'backward':
                    top = ' '*indentation + "5'> " + self._top + " >3'"
                    bottom = "3'< " + self._bottom + " <5'"
                    matchSigns = ' '*(indentation+4) + matchSigns
                    print(top)
                    print(matchSigns)
                    print(bottom)
                if direction == 'forward':
                    top = "5'> " + self._top + " >3'"
                    bottom = ' '*indentation + "3'< " + self._bottom + " <5'"
                    matchSigns = ' '*(indentation+4) + matchSigns
                    print(top)
                    print(matchSigns)
                    print(bottom)

        matches = detect_dimers()
        maxscore = find_maxscore(matches)
        # For debugging:
        #print_dimer(maxscore[1], maxscore[2], maxscore[3])
        return maxscore[0]
#end of Nguyen-Dumon et al., BioTechniques, Vol. 55, No. 2, August 2013, pp. 69

def is_part(array,idx):
    for i in array:
        if i == idx:
            return True
    else:
        return False

def fill_sequence_dictionary(seq_name, new_seq, current_sequences):
    updated_sequences = {}
    updated_sequences[seq_name] = new_seq
    for n in current_sequences.keys():
        updated_sequences[n] = current_sequences[n]

#all on base count level
def replace_codons(s,frame_end,overlap_len):

    # this is to keep track of the whole sequence:
    #remainder = s[:frame_end+1] #inclusive for frame_end
    #print "Remainder:"
    #print remainder

    #starting at the first position for new codon hence +1
    randomized_positions = [i for i in range((frame_end+1)/3, (frame_end+1+ overlap_len)/3 )]
    shuffle(randomized_positions)

    #change 6 codons at a time
    poss_to_change = min(6,len(randomized_positions)-1)
    changed = randomized_positions[:poss_to_change]
    a = sorted(changed)

    sh = ""
    for p in range(frame_end+1,frame_end+1+overlap_len, 3):
        if is_part(a,p/3):
            #Make a new codon
            codon = ''
            for i in range(3):
                codon += s[p + i]
            #grab a codon at random
            aa = str(Seq(codon).translate())
            r = randint(1,len(codons[aa]))
            sh += codons[aa][r-1]
        else:
            # Return original codon
            for i in range(3):
                sh += s[p + i]
    return sh

def filter_by_alignment_p(arg_list):
    rc = arg_list[0]
    s2 = arg_list[1]
    align_threshold = arg_list[2]
    sc1 = Dimer(rc, s2).score(align_threshold)
    if sc1 > align_threshold:
        return sc1
    else:
        return 0

def filter_by_alignment(rc,s2, align_threshold, container1):
    sc1 = Dimer(rc, s2).score(align_threshold)
    s2_revcomp = str(Seq(s2).reverse_complement())
    sc2 = Dimer(rc, s2_revcomp).score(align_threshold)
    if sc1 > align_threshold:
        container1.append(sc1)
    if sc2 > align_threshold:
        container1.append(sc2)

    return container1

def grow_overlap(startpoint, seq):
    '''
    Overlap grows from the middle outwards, regardless of the oligo length limit, that
    is checked afterwards, and if it's not good enough codons get swapped randomly until
    the overlap has the proper length and Tm. The returned stings are the three parts of
    sequence: five prime unique section, the overlap that will be later added to both,
    and the three-prime unique section.
    '''

    #seed for primer, start with min primer length
    #for i in range(startpoint, startpoint + min_primer_length):
    #    overlap += seq[i]

    #Minimum length of the overlap must be len = Tm/4, which assumes that it's %100 GC:
    min_len = minmelt/4
    overlap = seq[startpoint - min_len/2 : startpoint] + seq[startpoint : startpoint + min_len/2 ]
    counter = 0
    firsthalf = seq[:(startpoint - min_len/2 - counter)]
    tm = MeltingTemp.Tm_NN(overlap, Na=50, K=0, Tris=0, Mg=args.MgmM, dNTPs=args.dNTPsmM )
    gc = GC(overlap) 
    #while (startpoint + min_len/2 + counter < len(seq)) and ( (tm < minmelt) or ( (gc < 40) or (gc > 60) ) ):
    while (min_len/2 + counter < args.maxoverlaplen) and ( (tm < minmelt) or ( (gc < 40) or (gc > 60) ) ): # GC% must be between 40 and 60, Tm should be above mininimum, and the overlap should not be longer than 30bps
        counter += 1
        overlap = seq[startpoint - min_len/2 - counter] + overlap + seq[startpoint + min_len/2 + counter - 1]
        tm = MeltingTemp.Tm_NN(overlap, Na=50, K=0, Tris=0, Mg=args.MgmM, dNTPs=args.dNTPsmM )
        gc = GC(overlap)
    #print(counter)
    firsthalf = seq[:(startpoint - min_len/2 - counter)]
    secondhalf = seq[(startpoint + min_len/2 + counter):len(seq)]
    #print len(firsthalf + overlap + secondhalf), len(seq)
    #print(firsthalf)
    #print(overlap)
    #print(secondhalf)
    assert len(firsthalf + overlap + secondhalf) == len(seq)
    assert firsthalf + overlap + secondhalf == seq
    #assert str(Seq(seq, unambiguous_dna).translate()) == protein_seqs[design] # this is done later
    return firsthalf, overlap, secondhalf

def get_frame_end(startpoint):

    ## determine the frame, assuming input sequences given were all inframe
    frameend = startpoint

    if frameend % 3 == 0 :
        frameend -= 1
    elif frameend % 3 == 1:
        frameend += 1

    return frameend

def record_redesigned_sequence(mydesign,myseq,mynewseq,end,a1,a2):
    '''
    Record redesigned sequences to file and update sequence array
    (which will be evenutally formated with adapters)
    '''
    a1.append(myseq[0:end+1] + mynewseq)
    a2.append(mynewseq + myseq[(end +1) + len(mynewseq):len(myseq)])

    #print >> order_sequences, mydesign + '_1st' , myseq[0:end+1] + mynewseq
    #print >> order_sequences, mydesign + '_2nd' , mynewseq + myseq[(end +1) + len(mynewseq):len(myseq)]
    print >> final_sequences, mydesign , myseq[0:end+1] + mynewseq + seq[(end +1) + len(mynewseq):len(myseq)]

    #return a1,a2

'''
BEGGINING OF MAIN:
'''
# Load and process sequence files:
# for this verison, specify which line of the adapters should be used
# numbering starts from 1, not like python 0
lines=open( args.input_list).readlines()
adapter_idx = args.adaptor_number - 1
nproc = args.nproc

designs = {}
protein_seqs = {}
for line in lines:
    its = line.strip().split()
    if len(its) != 3:
        print("Skipping:")
        print(its)
        continue
    else:
    	designs[its[0]] = its[1]
        protein_seqs[its[0]] = its[2]

# Get adapters ready
five_prime_addition_5, five_prime_addition_3, three_prime_addition_5, three_prime_addition_3 = read_adapters(args.adaptor_fname)

#Load codons
codonfile = open(args.codontable_fname, 'r').readlines()
codons = {}
for line in codonfile:
    its = line.strip().split()
    try:
        if not codons[its[0]]:
            codons[its[0]].append(its[1])

        else:
            codons[its[0]].append(its[1])
    except KeyError:
        bases = []
        codons[its[0]] = bases
        codons[its[0]].append(its[1])

print "your condons\n", codons

newdesign_seqs = {}
final_seq = {}
trashed_seq = {}
order_seq = {}

curr_seq_stash = {}
for i in designs:
    curr_seq_stash[i] = designs[i]

bad_list = []

#container for final sequence to order
oligos5 = []
oligos3 = []
names = []

overlap_report_handle = open('overlap_report_%d.txt'%args.adaptor_number, 'w')

"""
End of setup, now go design by design and optimize overlaps:
"""
p = Pool(nproc)
for design in designs:
    print("Working on %s"%design)
    # Setup:
    seq = designs[design]
    accepted = False
    #startseq = len(seq)/2 - min_primer_length
    startseq = len(seq)/2
    # Checks:
    if 'K' in seq or 'R' in seq:
        print "it appears there are amino acids in the dna sequence, going to next design ..."
        continue

    # Generate initial overlap and flanking 3' and 5' sequences:
    fiveprime, overlapsegment, threeprime = grow_overlap(startseq, seq)
    # Calculate max possible length of oligo:
    five_prime_addition_3_len = len(five_prime_addition_3[adapter_idx])
    three_prime_addition_5_len = len(three_prime_addition_5[adapter_idx])
    yeast_primer5_len = len(yeast_primer5)
    yeast_primer3_len = len(yeast_primer3)
    current_len_5prime_oligo = yeast_primer5_len + len(fiveprime) + len(overlapsegment) + five_prime_addition_3_len
    current_len_3prime_oligo = three_prime_addition_5_len + len(overlapsegment) + len(threeprime) + yeast_primer3_len
    # Now make changes until the Tm is as high as desired and the complete oligo length is
    # in a viable range:
    for trial in range(1000):
	tm = MeltingTemp.Tm_NN(overlapsegment, Na=50, K=0, Tris=0, Mg=args.MgmM, dNTPs=args.dNTPsmM )
        gc = GC(overlapsegment) 
        if (tm >= minmelt) and ( (gc > 40) and (gc<60) ) and (len(overlapsegment) <= args.maxoverlaplen ):
            if not ((current_len_5prime_oligo > max_oligo_size ) or (current_len_3prime_oligo > max_oligo_size )):
                 break
        frame_end = get_frame_end(len(fiveprime))
        newseq = replace_codons(seq,frame_end,len(overlapsegment))
        seq = seq[:frame_end+1] + newseq + seq[frame_end+1+len(newseq):]
        #assert str(Seq(seq, unambiguous_dna).translate()) == protein_seqs[design]
        fiveprime, overlapsegment, threeprime = grow_overlap(startseq, seq)
        current_len_5prime_oligo = yeast_primer5_len + len(fiveprime) + len(overlapsegment) + five_prime_addition_3_len
        current_len_3prime_oligo = three_prime_addition_5_len + len(overlapsegment) + len(threeprime) + yeast_primer3_len
	seq = fiveprime+overlapsegment+threeprime
	#assert str(Seq(seq, unambiguous_dna).translate()) == protein_seqs[design]
    skip = False
    tm = MeltingTemp.Tm_NN(overlapsegment, Na=50, K=0, Tris=0, Mg=args.MgmM, dNTPs=args.dNTPsmM ); gc = GC(overlapsegment); overlapsegment_len = len(overlapsegment)
    if ((current_len_5prime_oligo > max_oligo_size ) or (current_len_3prime_oligo > max_oligo_size )) or (not ( (tm>= minmelt ) and ( (gc>40) and (gc<60) ) and (overlapsegment_len<= args.maxoverlaplen) ) ):
        print("Olifos do not pass requirements even after 1000 optimization attempts, discarding %s and placing on Impossible design list"%design)
        skip = True
        bad_list.append(design)
        trashed_seq[design] = seq
        print >> impossible_sequences, design, seq
    curr_seq_stash[design] = seq
    # Create reverse complement of the overlap to search for mispriming:
    overlapsegment_rc = str(Seq(overlapsegment).reverse_complement())
    missaligned1 = []
    
    result_fwd = p.map_async(filter_by_alignment_p,[ [overlapsegment_rc[::-1],curr_seq_stash[des],alignthreshold] for des in curr_seq_stash ])
    results_fwd = result_fwd.get()
    fwd = [ i for i in results_fwd if i > 0 ]
    result_rev = p.map_async(filter_by_alignment_p,[ [overlapsegment_rc[::-1],str(Seq(curr_seq_stash[des]).reverse_complement()),alignthreshold] for des in curr_seq_stash ])
    results_rev = result_rev.get()
    rev = [ i for i in results_rev if i > 0 ]
    missaligned1 = fwd + rev
    
    # If there was more than the one expected alignement of the overlap to the gene:
    if (len(missaligned1) > 1) and not skip:
        fixed = False
        print("Possible chimera detected, atempting to fix ...")
        for trial in range(20):
            frame_end = get_frame_end(len(fiveprime))
            newseq = replace_codons(seq,frame_end,len(overlapsegment))
            seq = seq[:frame_end+1] + newseq + seq[frame_end+1+len(newseq):]
            curr_seq_stash[design] = seq
            fiveprime, overlapsegment, threeprime = grow_overlap(startseq, seq)
            seq = fiveprime+overlapsegment+threeprime
            #assert str(Seq(seq, unambiguous_dna).translate()) == protein_seqs[design]
            overlapsegment_rc =  str(Seq(overlapsegment).reverse_complement())
            result_fwd = p.map_async(filter_by_alignment_p,[ [overlapsegment_rc[::-1],curr_seq_stash[des],alignthreshold] for des in curr_seq_stash ])
            results_fwd = result_fwd.get()
            fwd = [ i for i in results_fwd if i > 0 ]
            result_rev = p.map_async(filter_by_alignment_p,[ [overlapsegment_rc[::-1],str(Seq(curr_seq_stash[des]).reverse_complement()),alignthreshold] for des in curr_seq_stash ])
            results_rev = result_rev.get()
            rev = [ i for i in results_rev if i > 0 ]
            missaligned1 = fwd + rev

            current_len_5prime_oligo = yeast_primer5_len + len(fiveprime) + len(overlapsegment) + five_prime_addition_3_len
            current_len_3prime_oligo = three_prime_addition_5_len + len(overlapsegment) + len(threeprime) + yeast_primer3_len
            
            if (len(missaligned1) == 1) and not ((current_len_5prime_oligo > max_oligo_size ) or (current_len_3prime_oligo > max_oligo_size )):
                tm = MeltingTemp.Tm_NN(overlapsegment, Na=50, K=0, Tris=0, Mg=args.MgmM, dNTPs=args.dNTPsmM )
                gc = GC(overlapsegment)
                if (tm >= minmelt) and ( (gc > 40) and (gc<60) ) and (len(overlapsegment) <= args.maxoverlaplen ):
                    fixed = True
                    break

        if fixed:
            print("Success!")
            assert str(Seq(seq, unambiguous_dna).translate()) == protein_seqs[design]
            Tm = MeltingTemp.Tm_NN(overlapsegment, Na=50, K=0, Tris=0, Mg=args.MgmM, dNTPs=args.dNTPsmM )
            GC_cont = GC(overlapsegment)
            print("%s length: %d Tm: %02f GC: %d"%(overlapsegment,len(overlapsegment),Tm,GC_cont))
            print >> overlap_report_handle, "%s,%s,%d,%02f,%d"%(design,overlapsegment,len(overlapsegment),Tm,GC_cont)

            print >> final_sequences, design , seq
            names.append(design)
            final_seq[design] = seq
            oligos5.append(fiveprime + overlapsegment )
            oligos3.append(overlapsegment + threeprime)
        else:
            bad_list.append(design)
            print "Failed to remove all chimeras after 20 attempts, discarding %s and savign to ChimeraForming list" %design
            trashed_seq[design] = seq
            print >> trashed_sequences, design, seq

    elif not skip:
	assert str(Seq(seq, unambiguous_dna).translate()) == protein_seqs[design]
        Tm = MeltingTemp.Tm_NN(overlapsegment, Na=50, K=0, Tris=0, Mg=args.MgmM, dNTPs=args.dNTPsmM )
        GC_cont = GC(overlapsegment)
	print("%s length: %d Tm: %02f GC: %d"%(overlapsegment,len(overlapsegment),Tm,GC_cont))
        print >> overlap_report_handle, "%s,%s,%d,%02f,%d"%(design,overlapsegment,len(overlapsegment),Tm,GC_cont)
        final_seq[design] = seq
        oligos5.append(fiveprime + overlapsegment )
        oligos3.append(overlapsegment + threeprime)
        names.append(design)

        print >> final_sequences, design , seq

    ## update design pool after each processed design sequence
    curr_seq_stash.clear()
    for des in final_seq:
        curr_seq_stash[des] = final_seq[des]

    #then add all old designs minus the bad list
    for d in designs:
        if d in bad_list:
            continue
        else:
            curr_seq_stash[d] = designs[d]


trashed_sequences.close()
final_sequences.close()
impossible_sequences.close()
#order_sequences.close()

if len(oligos3) != len(oligos5):
    print "Number of 3' oligos (%i) is different from number of 5' oligos (%i)"%(len(oligos3),len(oligos5))


print len(oligos3),len(names)


###################### output formating section ... #######################


batch_small_id = adapter_idx
batch_big_id = adapter_idx
report = open('report_%d.txt'%args.adaptor_number, 'w')

distributions =[0] * len(five_prime_addition_5)

#gene sequence: up to 61 residues
#5' oligos:  (160bp)  17bp outside primer + (+ random bits) + 18bp yeast adaptor + gene half + overlap + AA + 17bp outside primer
#3' oligos:  (160bp)  17bp outside primer + TT + 30 overlap + 93 unique + 18bp yeast adaptor + ( random bits) +17bp outside primer

#gene sequence: up to 72 residues
#5' oligos:  (160bp)  18bp yeast adaptor +  unique + overlap + AA + 17bp outside primer (+ random bits)
#3' oligos:  (160bp)  17bp outside primer + TT + 30 overlap + 93 unique + 18bp yeast adaptor (+ random bits)

#using first 6 set for shorter genes and the last 6 for longer genes

## defining size limits -- asuming adapter are equal lengths, adjust as needed ####
upper_limit = max_oligo_size - len(yeast_primer5) - len(five_prime_addition_5[0])
lower_limit = max_oligo_size - len(yeast_primer5) - len(five_prime_addition_5[0]) -len(three_prime_addition_3[0])

#keeping track of small and big
type1 = 0
type2 = 0

for idx in range(len(oligos3)):

    fivePlen = len(oligos5[idx])
    threePlen = len(oligos3[idx])
    '''
    if fivePlen < lower_limit and threePlen  < lower_limit :
        type1 += 1

        distributions[batch_small_id] += 1

        ###  5' oligo  - adding random pieces between 5 addition and yeast adapter to adjust size to max length
        s1 = five_prime_addition_5[batch_small_id]
        addition = get_random_bases( max_oligo_size - len(five_prime_addition_5[batch_small_id]  + yeast_primer5 + oligos5[idx] + five_prime_addition_3[batch_small_id] ))

        s1 += addition
        s1 += yeast_primer5 + oligos5[idx] + five_prime_addition_3[batch_small_id]

        ### 3' oligo - adding after yeast adapter before 3' addition
        s2 = three_prime_addition_5[batch_small_id] + oligos3[idx]  + yeast_primer3
        addition = get_random_bases( max_oligo_size - len(three_prime_addition_5[batch_small_id] + oligos3[idx] + three_prime_addition_3[batch_small_id] + yeast_primer3))
        s2 += addition
        s2 += three_prime_addition_3[batch_small_id]

        print >>final_order, names[idx] + '_1st' , s1
        print >>final_order, names[idx] + '_2nd' , s2

	# update not necessary for this version, all one idx as specified by the user
        #batch_small_id += 1
        #if batch_small_id > 2 :#only the first 3 will be used for the smaller
        #    batch_small_id = 0
    '''
    if fivePlen < upper_limit and threePlen < upper_limit :

        type2 += 1

        #batch_big_id += 1

        #if batch_big_id > 9:# 5-10 , in python 4 - 9 of the adapters are used for the bigger designs
        #    batch_big_id = 5 # reset
        distributions[batch_big_id] += 1

        #adding to end random bases for customArray synthesis
        s1 = yeast_primer5 + oligos5[idx] + five_prime_addition_3[batch_big_id]
        #addition = get_random_bases( max_oligo_size - len(s1))
        #s1 = s1 + addition

        #if len(s2) < max_oligo_size+1:
        s2 = three_prime_addition_5[batch_big_id] + oligos3[idx] + yeast_primer3
        #addition = get_random_bases( max_oligo_size - len(s2))
        #s2 += addition
        #print len(s1), len(s2), max_oligo_size

        # If these asserts go off, you messed up something:
        assert len(s1) <= max_oligo_size # Final product len < max oligo len
        assert len(s2) <= max_oligo_size # Final product len < max oligo len
        print >>final_order_large, names[idx] + '_1st' , s1
        print >>final_order_large, names[idx] + '_2nd' , s2


    else:
        print names[idx], "is too big!!" , fivePlen, threePlen, "limits:" , lower_limit, upper_limit, oligos5[idx], oligos3[idx]
	print >> report , names[idx], "is too big!!" , fivePlen, threePlen, "limits:" , lower_limit, upper_limit, oligos5[idx], oligos3[idx]

print "double primer sequences: " , type1,"single", type2
print >> report ,"\nshort primers with double custom primers:", type1, "with single " , type2,"\n"


print distributions
print >> report,"adapter distributions:", distributions


report.close()
p.close()
p.join()

final_order_large.close()
overlap_report_handle.close()
exit()
