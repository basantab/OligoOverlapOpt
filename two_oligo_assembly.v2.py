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
#pETCON adaptors:
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
    """
    Shortcut to get reverse complement of DNA sequence
    """
    rc = Seq(seq, ambiguous_dna ).reverse_complement()
    return rc

def read_adapters(infile):
    """
    Reads in the adaptor file, obetaining the correct sequences from the expected format.
    """
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

def isComplement(x, y):
    """
    Returns true if the fist argument base is complentary to the base in the second argument.
    """
    cs = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}
    return cs[x] == y


def complement_positions(top, bottom):
    """
    Returns the position indexes (starts at 0) that are properly paired between the top and bottom sequences.
    """
    positions = []
    for i, (t, b) in enumerate(zip(top, bottom)):
        if isComplement(t, b):
            positions.append(i)
    return positions

def get_score(positions):
    """
    Score is the number of consecutive positions that paired in the input
    alignment positions.
    For this function see Nguyen-Dumon, BioTechniques, Vol 55, No 2, August 2013, pp 69
    """
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
    """
    For this class see Nguyen-Dumon, BioTechniques, Vol. 55, No. 2, August 2013, pp. 69
    This object is used in the chimera evaluation.
    """
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
    """
    Returns True if idx is in array
    """
    for i in array:
        if i == idx:
            return True
    else:
        return False


def replace_codons(s,frame_end,overlap_len):
    """
    Returns a new sequence encoding for the same peptide, but with different codons in the overlap region.
    """
    #starting at the first position for new codon hence +1
    randomized_positions = [i for i in range((frame_end+1)/3, (frame_end+1+ overlap_len)/3 )]
    shuffle(randomized_positions)

    #change 6 codons at a time
    poss_to_change = min(6,len(randomized_positions)-1)
    changed = randomized_positions[:poss_to_change]
    a = sorted(changed)

    sh = ""
    for p in range(frame_end+1,frame_end+1+overlap_len, 3):
        #if is_part(a,p/3):
        if p/3 in a:
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
    """
    Takes a senquence and the reverse complement of the current overlap, and searches for possible priming in other sequences.
    Returns the scores of the aligment that are above the input threshold. The score is just the number of bases properly aligned.
    """
    rc = arg_list[0]
    s2 = arg_list[1]
    align_threshold = arg_list[2]
    sc1 = Dimer(rc, s2).score(align_threshold)
    if sc1 > align_threshold:
        return sc1
    else:
        return 0

def grow_overlap(startpoint, seq):
    """
    Returns the sequence divided in three parts: 5', overlap and 3'.
    Overlap grows from the middle outwards, regardless of the oligo length limit, that
    is checked afterwards, and if it's not good enough codons get swapped randomly until
    the overlap has the proper length and Tm. The returned stings are the three parts of
    sequence: five prime unique section, the overlap that will be later added to both,
    and the three-prime unique section.
    """
    #Minimum length of the overlap must be len = Tm/4, which assumes that it's %100 GC:
    min_len = minmelt/4
    overlap = seq[startpoint - min_len/2 : startpoint] + seq[startpoint : startpoint + min_len/2 ]
    counter = 0
    firsthalf = seq[:(startpoint - min_len/2 - counter)]
    tm = MeltingTemp.Tm_NN(overlap, Na=50, K=0, Tris=0, Mg=args.MgmM, dNTPs=args.dNTPsmM )
    gc = GC(overlap) 
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
    """
    Returns the frame end given the that the sequence has length "startpoint"
    """
    frameend = startpoint

    if frameend % 3 == 0 :
        frameend -= 1
    elif frameend % 3 == 1:
        frameend += 1

    return frameend

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
    # Now make changes until the Tm is as high as desired, the GC% is between 40% and 60%, and the complete oligo length is
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
    # After 1000 attempts are done, check if sequence passes requierements:
    if ((current_len_5prime_oligo > max_oligo_size ) or (current_len_3prime_oligo > max_oligo_size )) or (not ( (tm>= minmelt ) and ( (gc>40) and (gc<60) ) and (overlapsegment_len<= args.maxoverlaplen) ) ):
        print("Oligos do not pass requirements even after 1000 optimization attempts, discarding %s and placing on Impossible design list"%design)
        skip = True
        bad_list.append(design)
        trashed_seq[design] = seq
        print >> impossible_sequences, design, seq
    curr_seq_stash[design] = seq
    # Create reverse complement of the overlap to search for mispriming:
    overlapsegment_rc = str(Seq(overlapsegment).reverse_complement())
    missaligned1 = []
    # Check for priming in other sequences:
    #result_fwd = p.map_async(filter_by_alignment_p,[ [overlapsegment_rc[::-1],curr_seq_stash[des],alignthreshold] for des in curr_seq_stash ])
    result_fwd = p.map_async(filter_by_alignment_p,[ [overlapsegment_rc[::-1],final_seq[des],alignthreshold] for des in final_seq ])
    results_fwd = result_fwd.get()
    fwd = [ i for i in results_fwd if i > 0 ]
    #result_rev = p.map_async(filter_by_alignment_p,[ [overlapsegment_rc[::-1],str(Seq(curr_seq_stash[des]).reverse_complement()),alignthreshold] for des in curr_seq_stash ])
    result_rev = p.map_async(filter_by_alignment_p,[ [overlapsegment_rc[::-1],str(Seq(final_seq[des]).reverse_complement()),alignthreshold] for des in final_seq ])
    results_rev = result_rev.get()
    rev = [ i for i in results_rev if i > 0 ]
    missaligned1 = fwd + rev
    
    # If there was more than the one expected alignement of the overlap to a gene:
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
            #result_fwd = p.map_async(filter_by_alignment_p,[ [overlapsegment_rc[::-1],curr_seq_stash[des],alignthreshold] for des in curr_seq_stash ])
            result_fwd = p.map_async(filter_by_alignment_p,[ [overlapsegment_rc[::-1],final_seq[des],alignthreshold] for des in final_seq ])
            results_fwd = result_fwd.get()
            fwd = [ i for i in results_fwd if i > 0 ]
            #result_rev = p.map_async(filter_by_alignment_p,[ [overlapsegment_rc[::-1],str(Seq(curr_seq_stash[des]).reverse_complement()),alignthreshold] for des in curr_seq_stash ])
            result_rev = p.map_async(filter_by_alignment_p,[ [overlapsegment_rc[::-1],str(Seq(final_seq[des]).reverse_complement()),alignthreshold] for des in final_seq ])
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

if len(oligos3) != len(oligos5):
    print "Number of 3' oligos (%i) is different from number of 5' oligos (%i)"%(len(oligos3),len(oligos5))


print len(oligos3),len(names)


batch_big_id = adapter_idx
report = open('report_%d.txt'%args.adaptor_number, 'w')
distributions =[0] * len(five_prime_addition_5)
upper_limit = max_oligo_size - len(yeast_primer5) - len(five_prime_addition_5[0])
lower_limit = max_oligo_size - len(yeast_primer5) - len(five_prime_addition_5[0]) -len(three_prime_addition_3[0])

#keeping track of small and big
type1 = 0
type2 = 0

for idx in range(len(oligos3)):

    fivePlen = len(oligos5[idx])
    threePlen = len(oligos3[idx])
    
    if fivePlen < upper_limit and threePlen < upper_limit :
        type2 += 1
        distributions[batch_big_id] += 1
        #adding to end random bases for customArray synthesis
        s1 = yeast_primer5 + oligos5[idx] + five_prime_addition_3[batch_big_id]
        s2 = three_prime_addition_5[batch_big_id] + oligos3[idx] + yeast_primer3

        # If these asserts go off, you messed up something:
        assert len(s1) <= max_oligo_size # Final product len < max oligo len
        assert len(s2) <= max_oligo_size # Final product len < max oligo len
        print >>final_order_large, names[idx] + '_1st' , s1
        print >>final_order_large, names[idx] + '_2nd' , s2


    else:
        print names[idx], "is too big!!" , fivePlen, threePlen, "limits:" , lower_limit, upper_limit, oligos5[idx], oligos3[idx]
	print >> report , names[idx], "is too big!!" , fivePlen, threePlen, "limits:" , lower_limit, upper_limit, oligos5[idx], oligos3[idx]

print >> report ,"\nshort primers with double custom primers:", type1, "with single " , type2,"\n"


print distributions
print >> report,"adapter distributions:", distributions


report.close()
p.close()
p.join()

final_order_large.close()
overlap_report_handle.close()
