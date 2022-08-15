from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord
from typing import List, Tuple, Optional
import Bio.SeqFeature as SeqFeature
import Bio.Align as Align

class CannotAlignError(Exception):
    pass

class ReferenceError(Exception):
    pass

class VariantAlignmentError(Exception):
    pass

def get_top_alignment(
    seqA: Seq,
    seqB: Seq,
    mode:  str = 'global',
    match_score: float = 1.0,
    mismatch_score: float = 0.0,
    open_gap_score: float = -.99,
    extend_gap_score: float = 0,
    ) -> Align.PairwiseAlignment:
    """Gets the top alignment between two sequences

    Parameters
    ----------
    seqA : Seq
        first sequence to align
    seqB : Seq
        second sequence to align
    mode : str, optional
        local or global alignment, by default 'global'
    match_score : float, optional
        added to alignment score upon match, by default 1.0
    mismatch_score : float, optional
        added to alignment score upon mismatch, by default 0.0
    open_gap_score : float, optional
        added to alignment score upon opening a gap, by default -1.0
    extend_gap_score : float, optional
        added to alignment score upon extending a gap, by default -.1

    Returns
    -------
    Align.PairwiseAlignment
        Top-scored alignment object

    Raises
    ------
    CannotAlignError
        raised if biopython fails to return a valid alignment
    """
    if isinstance(seqA, SeqRecord):
        seqA = seqA.seq
    if isinstance(seqB, SeqRecord):
        seqB = seqB.seq

    aligner = Align.PairwiseAligner(
        mode = mode,
        match_score = match_score,
        mismatch_score = mismatch_score,
        open_gap_score = open_gap_score,
        extend_gap_score = extend_gap_score)

    for _alignment in aligner.align(seqA, seqB):
        alignment = _alignment
        break

    if isinstance(alignment, Align.PairwiseAlignment):
        return alignment
    else:
        raise CannotAlignError()

def get_exons_from_gb(
    genbank: SeqRecord,
    gene_name: str,
    cds_num: int = 0,
    ) -> List[Seq]:
    """Returns a list of exon sequences (only used if a CDS is not explicitly passed)

    Parameters
    ----------
    genbank : SeqRecord
        genbank record
    gene_name : str
        name of gene contained in genbank
    cds_num : int, optional
        the coding sequence (by order) to use, by default 0

    Returns
    -------
    List[Seq]
        List of exon sequences

    Raises
    ------
    ValueError
        Raised if the cds number is out of range
    """
    #Get list of valid CDS regions for gene
    features = [f for f in genbank.features if f.type == 'CDS']
    if any(['gene' in f.qualifiers for f in features]):
        features = [f for f in features if gene_name in f.qualifiers['gene']]
    
    num_feats = len(features)
    if cds_num > num_feats:
        raise ValueError(f'Invalid cds_num, should be between 1 & {num_feats}')

    cds_annot = features[cds_num-1]

    exons = []
    for part in cds_annot.location.parts:
        exons.append(genbank.seq[part.start:part.end])
    
    return exons

def get_exon_alignment_from_gb(
    ref_amplicon: Seq,
    genbank: SeqRecord,
    gene_name: str,
    cds_num: int = 0,
    ) -> Align.PairwiseAlignment:
    """Aligns reference amplicon sequence to exon contained in
    genbank entry

    Parameters
    ----------
    ref_amplicon : Seq
        sequence of reference amplicon
    genbank : SeqRecord
        genbank record
    gene_name : str
        name of gene in genbank file
    cds_num : int, optional
        the coding sequence (by order) to use, by default 0

    Returns
    -------
    Align.PairwiseAlignment
        Top-scored alignment object between an exon and the reference amplicon

    Raises
    ------
    ReferenceError
        Raised if the reference amplicon is not alignable to the CDS
    """
    exons = get_exons_from_gb(
        genbank = genbank,
        gene_name = gene_name,
        cds_num = cds_num
        )

    # Get the list of possible exon alignments
    exon_alignments = []
    scores = []
    for exon in exons:
        try:
            pcr_exon_align = get_top_alignment(exon, ref_amplicon)
            exon_alignments.append(pcr_exon_align)
            scores.append(pcr_exon_align.score)
        except CannotAlignError:
            scores.append(0)
    if not(any(scores)): # If all scores are 0 raise error.
        raise ReferenceError(f'The reference amplicon:\n{ref_amplicon.seq}\ndoes not align to {gene_name} CDS {cds_num}.')

    # Find the best exon alignment over the whole CDS
    top_alignment_ind = scores.index(max(scores))
    target_exon_alignment = exon_alignments[top_alignment_ind]

    return target_exon_alignment

def _generate_exon_feature(
    coords: Tuple[tuple]
    )-> SeqFeature.FeatureLocation:
    """Helper function for generating exon annotations on reference amplicon

    Parameters
    ----------
    coords : Tuple[tuple]
        pair of nucleotide numbers corresponding to feature position start/end

    Returns
    -------
    SeqFeature.FeatureLocation
        Biopython FeatureLocation
    """
    exon_start = coords[0]
    exon_end = coords[1]
    # Create exon feature
    exon_loc = SeqFeature.FeatureLocation(
        start = exon_start,
        end = exon_end
        )
        
    return exon_loc


def annotate_ref_amplicon(
    ref_amplicon: Seq,
    genbank: SeqRecord,
    gene_name: str,
    codons_spanned: int,
    cds_num: int = 0,
    cds_fragment: Optional[Seq] = None
    ) -> SeqRecord:
    """Annotate coding featurs on reference amplicon

    Parameters
    ----------
    ref_amplicon : Seq
        sequence of reference amplicon
    genbank : SeqRecord
        genbank record
    gene_name : str
        name of gene in genbank file
    cds_num : int, optional
        the coding sequence (by order) to use, by default 0
    cds_fragment : Optional[Seq], optional
        explicitly passed sequence of the coding sequence targeted by reference
        amplicon, by default None

    Returns
    -------
    SeqRecord
        Annotated copy of reference amplicon

    Raises
    ------
    ReferenceError
        Raised if the amplicon and reference CDS align disjointly
    """
    # Get the sequences of all exonic regions for CDS
    try:
        alignment = get_top_alignment(cds_fragment, ref_amplicon, mode='local')
        print("Using indicated cds_fragment as targeted coding region")
    except ValueError:
        print(f'Getting targeted coding region for {gene_name} CDS {cds_num} from GenBank')
        alignment = get_exon_alignment_from_gb(
            ref_amplicon = ref_amplicon,
            genbank = genbank,
            gene_name = gene_name,
            cds_num = cds_num
            )


    # Get just the coordinates of the PCR product that map to the exon
    ref_pcr_aligned_coords = alignment.aligned[1]
    ref_pcr_aligned_coords = [i for i in ref_pcr_aligned_coords if i != None]

    # Get the parameters for the exon feature
    # Case where the alignment to the exon is a single continuous piece.
    if len(ref_pcr_aligned_coords) == codons_spanned:
        segments = []
        for coords in ref_pcr_aligned_coords:
            feature = _generate_exon_feature(
                coords = coords
                )
            segments.append(feature)

        if len(segments) > 1:
            exon_loc = SeqFeature.CompoundLocation(segments)
        else:
            exon_loc = segments[0]

        exon_feature = SeqFeature.SeqFeature(
            location = exon_loc,
            type="exon",
            location_operator = '',
            )

        if not(isinstance(cds_fragment, Seq)) and not(isinstance(cds_fragment, SeqRecord)):
            print("Determining CDS fragment from genbank..")
            frags = []
            for coords in ref_pcr_aligned_coords:
                frag = ref_amplicon[coords[0]:coords[1]]
                frags.append(str(frag))
            cds_fragment = Seq(''.join(frags))

        record = SeqRecord(
            seq = ref_amplicon,
            id = f'ref_amplicon',
            features = [exon_feature],
            annotations={'cds_fragment': cds_fragment}
            )

        return record

    # If there is more than one aligned segment, it means the alignment is
    # disjoint
    else:
        err = f'The reference amplicon:\n{ref_amplicon}\naligns to'
        err += f' {gene_name} CDS {cds_num} in a disjoint manner:'
        err += f'\n{alignment}'
        raise ReferenceError(err)


def get_cds_flanks(
    cds_fragment: Seq,
    ref_cds: Seq,
    ) -> Tuple[Seq, Seq]:
    """Returns the sequence abutting/flanking amplicon-targeted CDS fragment

    Parameters
    ----------
    cds_fragment : Seq
        explicitly passed sequence of the coding sequence targeted by reference
        amplicon
    ref_cds : Seq
        Full length reference coding sequence

    Returns
    -------
    Tuple[Seq, Seq]
        Pair of sequences to the left, or right of the amplicon-targeted CDS
        fragment

    Raises
    ------
    ReferenceError
        Raised if the concatenation of the fragment + flanks is not the same as
        the reference CDS
    """
    alignment = get_top_alignment(ref_cds, cds_fragment, mode='local')

    ref_cds_aligned = alignment.aligned[0]
    ref_cds_aligned = [a for a in ref_cds_aligned if a != None]

    left_start = ref_cds_aligned[0][0]
    right_start = ref_cds_aligned[-1][1]
    left = ref_cds[:left_start]
    right = ref_cds[right_start:]

    cds_whole = left + cds_fragment + right
    
    if cds_whole != ref_cds:
        alignment = get_top_alignment(cds_whole, ref_cds)
        raise ReferenceError("Reference not properly mapped to CDS")


    return (left, right)


def get_var_cds(
    variant: SeqRecord,
    amp_left: Seq,
    amp_right: Seq,
    ref_amplicon: SeqRecord,
    ) -> Tuple[Seq, Align.PairwiseAlignment]:
    """Align variant to CDS fragment and then determine a full-length,
    hypothetical coding sequence from this

    Parameters
    ----------
    variant : SeqRecord
        Sequence from FASTA of assembled contigs
    ref_cds : Seq
        Full legnth reference CDS
    ref_amplicon : SeqRecord
        Reference amplicon sequence

    Returns
    -------
    Seq
        Predicted coding sequence for the provided variant
    """
    # Annotate exonic section of variant and extract
    cds_fragment = ref_amplicon.annotations['cds_fragment']
    alignment = get_top_alignment(
        cds_fragment,
        variant,
        mode='local',
        extend_gap_score = 0
        )

    var_aligned = [a for a in alignment.aligned[1] if a != None]
    var_start = var_aligned[0][0]
    var_end = var_aligned[-1][1]
    
    var_cds = amp_left + variant.seq[var_start:var_end] + amp_right

    return (var_cds, alignment)