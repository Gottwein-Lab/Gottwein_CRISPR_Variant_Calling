import Bio.Align as Align
from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq
from typing import Optional,Tuple, Dict, Union
from dataclasses import dataclass
from annotation_alignment_utils import get_var_cds, get_top_alignment,\
    VariantAlignmentError

CharacterizationReport = Dict[str, Optional[Union[str, int]]]

@dataclass
class FrameshiftResult:
    is_shifted: bool
    shift_amount: int
    shift_pos: Optional[int] = None

def is_indel(
    var_cds: Seq,
    ref_cds: Seq
    ) -> Tuple[bool, Align.PairwiseAlignment]:
    """Tests if the variant contains indels
    (if the alignment contains more than one segment)

    Parameters
    ----------
    var_cds : Seq
        predicted coding sequence for the variant
    ref_cds : Seq
        reference coding sequence

    Returns
    -------
    Tuple[bool, Align.PairwiseAlignment]
        Whether an indel was found and the alignment used to determine this
    """
    alignment = get_top_alignment(var_cds, ref_cds)
    var_aligned = alignment.aligned[0]
    var_aligned = [a for a in var_aligned if a != None]

    _is_indel =  len(var_aligned) > 1

    return (_is_indel, alignment)

def is_frameshifted(
    aligned: Tuple[tuple],
    max_codon_shift: int,
    ) -> FrameshiftResult:
    """
    Iterates over alignment regions and tests if the cumulative sum of indel
    regions is not divisible by three for more than 3*max_codon_shift
    nucleotides total

    Parameters
    ----------
    aligned : Tuple[tuple]
        Pair of alignment coordinates (1 or more); first set is variant,
        followed by reference
    max_codon_shift : int
        maximum number of codons to allow out-of-frame before the variant is
        called a frameshift

    Returns
    -------
    FrameshiftResult
        Information about whether the variant is frameshifted and if so by how
        much/where
    """
    var_aligned = [nt for nt in aligned[0] if nt != None]
    ref_aligned = [nt for nt in aligned[1] if nt != None]
    
    shift_amount = 0
    # These variables will be used to track if we exceed the maximum
    # allowed codon shifted
    left_shift_pos = None
    right_shift_pos = None
    is_shifted = False

    for var_nt, ref_nt in zip(var_aligned, ref_aligned):
        #Check if the alignments are in sync
        var_start = var_nt[0]
        ref_start = ref_nt[0]
        shift = var_start - ref_start
        shift_pos = ref_start
        if shift != 0:
        # If shifted, update the amount
            shift_amount += shift
            
            if shift_amount % 3 != 0:
                # If shift changes frame, record how long it is
                right_shift_pos = ref_start
                if left_shift_pos == None:
                    left_shift_pos = ref_start
            else:
                # If this shift does not change frame or stops changing the frame
                # Check if max_codon_shift has been exceeded
                if left_shift_pos != None:
                    # We only need to check this if this isn't the first shift
                    left_shift_codon = left_shift_pos // 3
                    right_shift_codon = right_shift_pos // 3
                    codons_shifted = right_shift_codon - left_shift_codon
                    # If we have exceeded the max codons, record and break
                    if codons_shifted >= max_codon_shift:
                        shift_amount = right_shift_pos - left_shift_pos
                        shift_pos = left_shift_pos
                        return FrameshiftResult(True, shift_amount%3, shift_pos//3)
                    # Otherwise, let's start tracking anew
                    else:
                        left_shift_pos = right_shift_pos        
        else:
        # If we have restored frame, we should reset the shift amount and pos
            shift_amount = 0
            shift_pos = None
    
    try:
        shift_pos = shift_pos // 3
    except TypeError:
        pass
    
    shift_amount = shift_amount % 3
    is_shifted = shift_amount != 0

    return FrameshiftResult(is_shifted, shift_amount, shift_pos)

def get_mut_type(
    var_peptide: Seq,
    ref_peptide: Seq
    ) -> Tuple[str, Optional[int]]:
    """Returns a string describing the type (or lackthereof) mutation and the
    position of the first predicted stop codon.

    Parameters
    ----------
    var_peptide : Seq
        Predicted variant peptide sequence
    ref_peptide : Seq
        Reference peptide sequence

    Returns
    -------
    Tuple[str, Optional[int]]
    mutation type and position of first stop codon (None if nonstop)
        
    """
    var_stop = var_peptide.find('*')

    if var_peptide == ref_peptide:
        mut_type = 'synonymous'
    else:
        if var_stop == -1:
            mut_type = 'nonstop'
            var_stop = None
        else:
            var_peptide = var_peptide[0:var_stop]
            if len(var_peptide) == len(ref_peptide):
                mut_type = 'missense'
            else:
                mut_type = 'nonsense'

    return (mut_type, var_stop)

def characterize_variant(
    variant: SeqRecord,
    ref_cds: Seq,
    ref_amplicon: SeqRecord,
    amp_left: Seq,
    amp_right: Seq,
    ref_peptide: Seq,
    max_codon_shift: int
    ) -> Tuple[CharacterizationReport, Align.PairwiseAlignment]:
    """Classifies a variant in terms of indel/frameshift/overall mutation

    Parameters
    ----------
    variant : SeqRecord
        Raw variant sequence, will be converted to CDS based on alignment
    ref_cds : Seq
        reference coding sequence
    ref_amplicon : SeqRecord
        reference amplicon sequence
    ref_peptide : Seq
        reference peptide sequenc

    Returns
    -------
    Tuple[CharacterizationReport, Align.PairwiseAlignment]
        Description of variant, including indel/frameshift status, mutation
        type (if any), stop codon location, number of read pairs, and frequecy
        of variant + the alignment generated between the variant and reference
    """
    # ref_peptide is passed instead of translated below to minimize the number
    # of times this is done
    results = {}
    try:
        var_cds, var_alignment = get_var_cds(
            variant = variant,
            amp_left = amp_left,
            amp_right = amp_right,
            ref_amplicon = ref_amplicon
            )
        
        results['variant_cds'] = str(var_cds)

        results['indel'], cds_alignment = is_indel(var_cds, ref_cds)
    
        frameshift_res = is_frameshifted(cds_alignment.aligned, max_codon_shift)
        results['frameshift'] = frameshift_res.is_shifted
        results['frameshift_amount'] = frameshift_res.shift_amount

        if (len(var_cds) % 3 != 0):
            padding = 'N'*frameshift_res.shift_amount
            padding = Seq(padding)
            var_peptide = var_cds+padding
            var_peptide = var_peptide.translate()
        else:
            var_peptide = var_cds.translate()

        results['peptide'] = str(var_peptide)

        results['mut_type'], results['stop_loc'] = get_mut_type(
            var_peptide,
            ref_peptide
            )
    except VariantAlignmentError:
        keys = [
            'indel',
            'frameshift',
            'frameshift_pos',
            'frameshift_amount',
            'stop_loc'
            ]
        results['mut_type'] = 'no_alignment'
        for k in keys:
            results[k] = None

    descr_str = variant.description.split(' ')
    descr_str = [i for i in descr_str if i != '']
    results['read_pairs'] = int(descr_str[1])
    results['frequency'] = float(descr_str[-1].strip('%'))/100

    return (results, var_alignment)