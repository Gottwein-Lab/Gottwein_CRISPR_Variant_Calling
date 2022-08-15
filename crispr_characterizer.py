#!/usr/bin/env python
import os
import csv
import json
from Bio import SeqIO
import math
import Bio.Align as Align
from Bio.SeqIO import SeqRecord
from pathlib import Path
import pathlib
from Bio.Seq import Seq
from typing import  Dict, Union, Optional, Tuple
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
from characterization import characterize_variant, CharacterizationReport
from annotation_alignment_utils import annotate_ref_amplicon,\
    get_exons_from_gb, get_cds_flanks, ReferenceError

CHAR_KWARGS = [
    'var_fasta_path',
    'ref_amplicon',
    'amplicon_strand'
    ]

REF_KWARGS = [
    'genbank_path',
    'ref_amplicon',
    'gene_name',
    'cds_num',
    'amplicon_strand',
    'cds_strand',
    'cds_fragment',
    'codons_spanned'
    ]

MUT_TYPES = [
    'synonymous',
    'missense',
    'nonsense',
    'nonstop',
    'no_alignment'
    ]

CWD = pathlib.Path(__file__).parent.absolute()

class CodonError(Exception):
    pass

def run_crispresso(
    fq_path: Path,
    ref_amplicon: Seq,
    sgrna_seq: Seq,
    gene_name: str,
    locus: str,
    cds_fragment: Seq
    ):
    """Function to run CRISPResso2 by issuing a system terminal command

    Parameters
    ----------
    fq_path : Path
        path to fastq file
    ref_amplicon : Seq
        reference amplicon sequence
    sgrna_seq : Seq
        sequence of sgrna
    gene_name : str
        gene name, used for output name
    locus : str
        CRIPSR guide locus identifier, used for output name
    cds_fragment : Seq
        sequence of the fragment of the codoing sequence hit by amplicon

    Returns
    -------
    Tuple[float, float]
    Pair of values corresponding to frequency of editing and frameshifts
    """
    name = f'{gene_name}_{locus}'
    outdir = CWD/"Crispresso2"

    cmd = f'CRISPResso --fastq_r1 "{fq_path}" --split_interleaved_input --amplicon_seq "{ref_amplicon}" --guide_seq "{sgrna_seq}" -gn "{locus}" --quantification_window_size 10 --quantification_window_center -3 --base_editor_output --coding_seq "{cds_fragment}" --name "{name}" --output_folder {outdir} --exclude_bp_from_left 5 --exclude_bp_from_right 5 --no_rerun'
    print(cmd)
    os.system(cmd)

    outdir = outdir/f'CRISPResso_on_{name}'

    editing_freq = outdir/'CRISPResso_quantification_of_editing_frequency.txt'
    editing_freq = pd.read_csv(editing_freq, sep='\t', header=0)
    editing_freq = editing_freq['Modified']/editing_freq['Reads_aligned']

    frameshift_freq = outdir/'Frameshift_analysis.txt'
    with open(frameshift_freq, 'r') as f:
        frameshift_freq = f.read()
    
    frameshift_freq = frameshift_freq.split('\n')
    frameshift_freq = [frameshift_freq[2], frameshift_freq[3]]
    frameshift_freq = [i.split(':')[1] for i in frameshift_freq]
    frameshift_freq = [int(i.split(' ')[0]) for i in frameshift_freq]
    frameshift_freq = frameshift_freq[1]/(frameshift_freq[0]+frameshift_freq[1])

    return (editing_freq.values[0], frameshift_freq)

def strand_check(
    target: Union[Seq, SeqRecord],
    strand: str,
    rc_kwargs: dict = {}
    ) -> Union[Seq, SeqRecord]:
    """Helper function that corrects biopython sequence/sequence record so it
    is in the proper sense.

    Parameters
    ----------
    target : Union[Seq, SeqRecord]
        sequence to check and fix if not properly stranded
    strand : str
        strand which can be considered "sense" (coding)

    Returns
    -------
    Union[Seq, SeqRecord]
        Correctly stranded sequence

    Raises
    ------
    ValueError
        Raised if strand passed is not + or -
    """
    if strand not in (['+','-']):
        raise ValueError(f'Strand must be + or - but {strand} was provided')
    if strand == '-':
        return target.reverse_complement(**rc_kwargs)
    else:
        return target

def get_ref_data(
    genbank_path: Path,
    cds_strand: str,
    ref_amplicon: Seq,
    amplicon_strand: str,
    gene_name: str,
    codons_spanned: int,
    cds_num: int = 0,
    cds_fragment: Optional[Seq] = None
    ) -> Tuple[SeqRecord, SeqRecord]:
    """Get reference CDS sequence

    Parameters
    ----------
    genbank_path : Path
        path to genbank file
    cds_strand : str
        strand on which the cds is found
    ref_amplicon : Seq
        sequence of reference amplicon
    amplicon_strand : str
        strand on which the translation proceeds 5' to 3' for the amplicon
    gene_name : str
        name of the gene in the genbank entry
    codons_spanned : int
        number of codons which the amplicon spans, typically 1
    cds_num : int, optional
        the coding sequence (by order) to use, by default 0
    cds_fragment : Optional[Seq], optional
        explicitly passed sequence of the coding sequence targeted by reference
        amplicon, by default None

    Returns
    -------
    Tuple[SeqRecord, SeqRecord]
        Reference coding sequence and reference amplicon
    """
    ref_amplicon = strand_check(ref_amplicon, amplicon_strand)

    gb = SeqIO.read(genbank_path, format='genbank')
    gb = strand_check(gb, cds_strand)

    exons = get_exons_from_gb(gb, gene_name, cds_num)
    exons = [str(seq) for seq in exons]
    ref_cds = Seq(''.join(exons))
    ref_amplicon = annotate_ref_amplicon(
        ref_amplicon = ref_amplicon,
        genbank = gb,
        gene_name = gene_name,
        cds_num = cds_num,
        cds_fragment = cds_fragment,
        codons_spanned = codons_spanned
        )

    return (ref_cds, ref_amplicon)


def characterize_variants(
    var_fasta_path: Path,
    ref_cds: Seq,
    ref_amplicon: SeqRecord,
    ref_peptide: Seq,
    amplicon_strand: str,
    max_codon_shift: int
    ) -> Tuple[Dict[str, CharacterizationReport], Dict[str, Align.PairwiseAlignment]]:
    """Takes a FASTA file of assembled contigs and characterizes each one.

    Parameters
    ----------
    var_fasta_path : Path
        path to variant fasta file, should be the set of assembled contigs
    ref_amplicon : str
        sequence of the reference amplicon
    gene_name : str
        namne of gene in genbank entry
    cds_num : int, optional
        the coding sequence (by order) to use, by default 0
    cds_strand : str
        strand on which the cds is found
    amplicon_strand : str
        strand on which the translation proceeds 5' to 3' for the amplicon
    cds_fragment : Optional[Seq], optional
        explicitly passed sequence of the coding sequence targeted by reference
        amplicon, by default None

    Returns
    -------
    Tuple[dict, dict]
    Pair of dictionaries containign the characterization report and alignment
    for each variant
        
    """
    variants = {}

    amp_left, amp_right = get_cds_flanks(
        cds_fragment = ref_amplicon.annotations['cds_fragment'],
        ref_cds = ref_cds)
        
    alignments = {}
    with open(var_fasta_path) as filestream:
        for i,variant in enumerate(SeqIO.parse(filestream, 'fasta')):
            if i > 0:
                variant = variant.upper()
                if amplicon_strand == '-':
                    rc_kwargs = {'description': True, 'id': True, 'name': True}
                    variant = strand_check(variant, amplicon_strand, rc_kwargs)
                results, alignment = characterize_variant(
                    variant = variant,
                    ref_cds = ref_cds,
                    ref_amplicon = ref_amplicon,
                    amp_left = amp_left,
                    amp_right = amp_right,
                    ref_peptide = ref_peptide,
                    max_codon_shift = max_codon_shift
                    )
                variants[variant.id] = results

                alignments[variant.id] = alignment
        
    return [variants, alignments]

def summarize(
    run_id: str,
    variants: Dict[str, CharacterizationReport],
    alignments: Dict[str, Align.PairwiseAlignment],
    ref_peptide_len: int,
    out_path: Path
    ) -> Dict[str, Union[int, str]]:
    """Generates summary files/reports

    Parameters
    ----------
    run_id : str
        Identifier string for run
    variants : Dict[str, CharacterizationReport]
        Results dictionary for variants characterization
    alignments : Dict[str, Align.PairwiseAlignment]
        Results dictionary for variant alignments
    ref_peptide_len : int
        length of the reference peptide
    out_path : Path
        path to deposit output files

    Returns
    -------
    Dict[str, Union[int, str]]
        Dictionary containing summary of results
    """
    df = pd.DataFrame.from_dict(variants, orient='index')
    total_read_pairs = df['read_pairs'].sum()
    total_freq = df['frequency'].sum()

    print(f'{total_read_pairs} total reads characterized, representing {total_freq*100}% of NGS reads')

    indels = df['indel']
    indel_rate = df.loc[indels, 'read_pairs'].sum()  / total_read_pairs
    frameshifts = df['frameshift']
    frameshift_rate = df.loc[frameshifts, 'read_pairs'].sum()  / total_read_pairs
    
    mean_stop_loc = df['stop_loc'].mean(skipna=True)
    median_stop_loc = df['stop_loc'].median(skipna=True)
    mode_stop_loc = np.mean(df['stop_loc'].mode(dropna=True).values)
    
    locs = []
    sums = df.groupby(['stop_loc'])['read_pairs'].aggregate('sum')
    stop_loc = sums.index
    counts = sums.values
    sums = {'stop_loc': stop_loc, 'read_pairs': counts}
    sums = pd.DataFrame.from_dict(sums)

    for idx in sums.index:
        loc = sums['stop_loc'][idx]
        for _ in range(sums['read_pairs'][idx]):
            locs.append(loc)
    locs = np.array(locs)
    norm_locs = locs/ref_peptide_len

    sns.set(font="Arial")
    sns.set_style("whitegrid")
    sns.set_context("talk", font_scale=.9)
    sns.kdeplot(norm_locs, fill=True)
    plt.title(f'{run_id} stop codons')
    plt.ylabel('Read Density')
    plt.xlabel(f'Rel. Amino Acid Position\nRef Len = {ref_peptide_len}')
    plt.xlim((0,1.05))
    plt.tight_layout()
    plt.savefig(out_path+'_norm'+'.pdf')
    plt.close()



    summary = {
        'ref_peptide_len': ref_peptide_len,
        'indel_rate': indel_rate,
        'frameshift_rate': frameshift_rate,
        'total_read_pairs': total_read_pairs,
        'total_ngs_freq': total_freq,
        'total_variants': len(df),
        'mean_stop_loc': mean_stop_loc,
        'median_stop_loc': median_stop_loc,
        'mode_stop_loc': mode_stop_loc
        }

    for mut in MUT_TYPES:
        mut_rows = df['mut_type'].isin([mut])
        mut_rate = df.loc[mut_rows, 'read_pairs'].sum() / total_read_pairs
        summary[f'{mut}_rate'] = mut_rate

    sum_text = 'Indel rate is {indel_rate}'
    sum_text += ', frameshift rate is {frameshift_rate}'
    sum_text += ', nonsense rate is {nonsense_rate}'
    print(sum_text.format(**summary))

    df.to_csv(out_path+'.csv')
    with open(out_path+'.json', 'w') as f:
        f.write(json.dumps(variants))

    alignment_str = '\n'.join([f">{k}{str(v.aligned[1])}\n{str(v)}" for k,v in alignments.items()])
    with open(out_path+'.pwa', 'w') as f:
        f.write(alignment_str)

    return summary


def main(
    run_table: Path,
    max_codon_shift: int
    ):
    """Main function, takes in run table, performs characterization and
    generates reports

    Parameters
    ----------
    run_table : Path, str
        path to the run table describing the edited pool to be characterized

    Raises
    ------
    ReferenceError
        Raised if the reference CDS and reference amplicon do not align properly
    """
    print(f'Starting run from {run_table}')
    reader = csv.DictReader(open(run_table))
    
    summary = {}
    for row in reader:

        run_params = {}
        for col, val in row.items():
            if col in ['cds_num', 'codons_spanned']:
                val = int(val)
            if col in ['var_fasta_path', 'genbank_path', 'fq_path']:
                val = Path(val)
            if col in ['ref_amplicon', 'cds_fragment', 'sgrna_seq']:
                val = Seq(val).upper()
            
            run_params[col] = val
        
        if 'cds_fragment' not in run_params.keys():
            run_params['cds_fragment'] = None
        else:
            run_params['cds_fragment'] = strand_check(
                run_params['cds_fragment'],
                run_params['cds_strand']
            )

        run_id = f"{run_params['gene_name']}_{run_params['locus']}"
        print(f'Running {run_id}...')
        
        char_kwargs = {k:v for k,v in run_params.items() if k in CHAR_KWARGS}
        ref_kwargs = {k:v for k,v in run_params.items() if k in REF_KWARGS}

        try:
            ref_cds, ref_amplicon = get_ref_data(**ref_kwargs)
            char_kwargs['ref_cds'] = ref_cds
            char_kwargs['ref_amplicon'] = ref_amplicon
        except ReferenceError:
            raise ReferenceError(f'Reference error  {run_id}')

        ref_len = len(ref_cds)
        if ref_len % 3 != 0:
            raise CodonError(f'Ref CDS not divisible by three: {ref_len}')
        ref_peptide = ref_cds.translate()

        variants, alignments = characterize_variants(
            max_codon_shift=max_codon_shift,
            ref_peptide = ref_peptide,
            **char_kwargs
            )

        mod_rate, frameshift_rate = run_crispresso(
            fq_path=run_params['fq_path'],
            ref_amplicon=char_kwargs['ref_amplicon'].seq,
            sgrna_seq=run_params['sgrna_seq'],
            gene_name=run_params['gene_name'],
            locus=run_params['locus'],
            cds_fragment=run_params['cds_fragment']
            )

        out_path = run_params['out_path']
        _summary = summarize(run_id, variants, alignments, len(ref_peptide), out_path)
        _summary['CRISPResso2_pct_modified'] = mod_rate
        _summary['CRISPResso2_pct_frameshift'] = frameshift_rate
        summary[run_id] = _summary
    
    sum_df = pd.DataFrame.from_dict(summary, orient='index')
    sum_df.to_csv(run_table.parents[0]/'summary.csv')




if __name__ == '__main__':
    import sys
    run_table = Path(sys.argv[1])
    _max_codon_shift= sys.argv[2]

    float_err = False
    try:
        max_codon_shift = int(_max_codon_shift)
        flt = float(_max_codon_shift)
        if math.ceil(flt) != max_codon_shift:
            float_err = True
    except ValueError as e:
        float_err = True
    if float_err:
        raise ValueError(f'Max codon shift must be a whole number: passed {_max_codon_shift}')

    main(run_table, max_codon_shift)