# Viral RNA simulation

Here is a very simplistic and idealized (Python) simulation of viral RNA
replication in cells.

As it stands, this is a simulation of a _positive-sense_ RNA viral entering a
cell and replicating through negative (reverse complement (RC)) intermediate
copies.

The aim here was to examine the extent to which, even in a highly idealized
scenario, the actual transcription mutations (in negative strand RNA) in the
replication process could be mistaken for other mutations.

The simulation has _many_ assumptions. Among them:

* RNA molecules do not degrade over time.
* Genetic changes only occur during polymerase transcription (as opposed to
  occurring spontaneously in RNA molecules as a result of some other process,
  e.g., SNPs or APOBEC editing).
* All RNA molecules (positive and negative sense) in all cells are sequenced.
* In sample preparation for sequencing, all RNA molecules are first reverse
  transcribed to make a single-stranded DNA (ssDNA) and a polymerase then
  makes the complementary DNA strand, resulting in a double-stranded DNA
  (dsDNA) that is sequenced.
* Reads from the sequencing are aligned against the (+) RNA reference genome
  of the virus.

Sub-genomic RNAs (sgRNA) are not part of the model. Including them would not
have shed any light on the original question, since changes in sgRNA are not
interpreted any differently from from changes in full-genome copies. Sample
preparation and sequencing of both is identical.

## How should we count "mutations"?

Suppose a (+) RNA molecule is being copied by a polymerase, that a `C`
nucleotide is encountered, and the polymerase accidentally incorporates an
`A` (instead of a `G`) into the negative strand.

There are three ways we might describe this error, all of them questionable.

### 1. Not useful or accurate

It could be classified as a `C->A` change. The polymerase is adding a
nucleotide that should complement a `C` and it mistakenly incorporates an
`A`. Although a `C` is involved, it is not changed and so it feels misleading
label this as a `C` change (but see below).

### 2. Logically attractive?

It could instead be called a `G->A` change, because the base that should have
been added was a `G`.

Naming it this way makes some logical sense because the polymerase might do
its job just fine, but something else (e.g., APOBEC) subsequently changes the
correct `G` into an `A`. That would result in the same (-) RNA as if the
polymerase had made the error, but in this case the change is unequivocally a
`G->A` change.

Because the result of these two different processes is a (-) RNA with an `A`,
we might want to use the same term to describe what happened and call this a
`G->A` change.

### 3. Standard practice

Under ideal conditions, if the (-) RNA with the erroneous `A` is sequenced,
the sample preparation will result in a dsRNA molecule which will have an
`A/T` pair.  The first ssDNA strand will be the RC of the (-) RNA (which is a
RC of the (+) RNA genome) and hence it is a 5' to 3' match of the viral
genome, but with a `T` (RC of the `A` in the (-) RNA).

The second ssDNA synthesized will be a RC of the first ssDNA and have an `A`.
Its RC will match the reference (it is the same as the (-) RNA, which is a RC
of the viral genome).

Both strands of the dsDNA will be sequenced and the resulting reads will be
aligned to the (+) RNA viral reference. 

Alignment software (at least including
[bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and
[bwa](https://bio-bwa.sourceforge.net/)) will put the RC of the second dsDNA
strands (which is a RC of the viral genome) into their resultant
[SAM](https://en.wikipedia.org/wiki/SAM_(file_format)) file, with the
`0x10 REVERSE` bit set in the flag field. I.e., the sequence in the SAM file
for the second DNA strand will match the reference genome, just like the
sequence for the first DNA strand does. I.e., both sequences in the SAM file
will have a `T`.

Then, software like
[pysam](https://pysam.readthedocs.io/en/stable/) will indicate that the read
had a `T` at this location (RNA viral reference genomes use `T` instead of
the `U` that is actually in the RNA and sequencing is of a DNA copy of the
viral RNA, so it makes some sense to here talk of `T` instead of `U`).

In this case, although the polymerase should have incorporated a `G` but
instead incorporated an `A`, this change would be classified as a
`C->U`. However, it is clearly not accurate to think of `C->U` as being the
_actual_ physical mutation. Instead, the result of the change is, supposing
the (-) RNA is subsequently transcribed, that the cell will then contain a
(+) RNA with a `U` at the genome site that in the original (+) RNA there is a
`C`.

People will informally say there has been a `C->U` mutation, but this is not
the case. (It is only the case if a `C` on a (+) RNA is directly mutated to
become a `U`.)

## This is unfortunately all a bit ambiguous

If we want to get a grip on the physical changes occurring in cells with (+)
RNA viral infections, there are at least four basic scenarios that it would
be good to be able to distinguish between. There are 1) misincorporation
errors made in the process of polymerase copying and 2) point mutations
(where one nucleotide is altered in place). These two can result in a change
in a positive or a negative RNA molecule.

Unfortunately, under standard practice these four scenarios cannot be
distinguished.

Here are the main reasons for this:

* Absent other information, we cannot distinguish between changes in RNA
  molecules that are made due to mutations from those due to misincorporation
  errors.  By "other information", I mean for example knowing the relative
  surrounding genome context of changes, aysmmetries in which might implicate
  a template-constrained changes by APOBEC enzymes.
* During sample preparation for sequencing, RNA molecules (both positive and
  negative) are used to produce dsRNA. Any kind of genetic change that
  results in an `A` or a `T`, in either of a (+) or (-) RNA, is going to look
  the same. This is because the sequencing reads from the dsDNA will be
  compared to the viral reference. The same goes for `C` and `G` changes.
* On top of this, we always align sequencing reads against the viral
  reference, so if the reference contains a `C`, then any change at that site
  is inevitably going to be classified as a `C->` whether or not a `C`
  nucleotide was involved (if the change occurs on a (-) RNA, then a `G` is
  involved, or if a (-) RNA has already been changed and that molecule is
  copied and another change happens, completely different nucleotides might
  be involved, with no `C` in the picture at all. But we'd still call it a
  `C->` change because we compare to the reference.

So we can't tell the difference between these pairs (`A/T` or `C/G`), we
can't tell on which RNA strand they occurred, we don't distinguish between
misincorporation errors and mutation, and we always use a naming convention
that makes it look like a reference nucleotide was changed. Pretty dismal.

## What the simulation is trying to assess

The simulation counts the actual changes that occur during polymerization and
compares this underlying reality to what appears to have happened once
standard sample preparation, sequencing, and mutation frequency calling are
done.

The simulation is idealized, with many deliberate simplifications, to
highlight that even in the simplest possible setup it is possible to draw
entirely wrong conclusions.

For example, suppose, as above, that the infecting virus has a `C` and the
polymerase puts an `A` into the (-) RNA RC. That is a single
error. Now suppose that (-) RNA is transcribed 100 times (a ratio of positive
to negatve RNA of 10-100:1 is known for coronaviruses) with no errors. Then
we have one original (+) RNA molecule, one (-) RNA with an `A`, and 100 (+)
RNA with `U`. Now we do sample preparation and sequence everything. The
original molecule results in a dsDNA with a `C/G` pair. The (-) RNA results
in a dsDNA with an `A/T` pair, and the 100 (+) RNAs also each result in a
dsDNA with an `A/T` pair. That's 102 dsDNA molecules. In the sequencing these
get denatured into 204 ssDNA molecules and sequenced. The original (+) RNA
yields reads for two ssDNA that both match the reference (one in RC). Then we have 200 ssDNA which are all aligned to the
reference. Half have a `T` and half have an `A`, and these are all classified
and counted as `C->T` changes.

In this not-unrealistic scenario, a single change (arguably a `C->U`, under
classification #3 above) has resulted in a read depth at that site of 200.

Suppose on the other hand that the polymerase correctly incorporated a `G` in
making the (-) RNA but then while making the 100 (+) RNA from the (-) RNA it
made misincorporation errors five times, each time incorporating a `G` (or an
`A`) instead of the correct `C`. Those errors will show up as ten `C->G` (or
`C->A`) reads when the sequencing data are aligned to the reference. The five
actual polymerase error in making the (+) RNAs will be totally swamped by the
single error in making one (-) RNA (that is then copied many times).

The simulation is designed to show how far the apparent changes (after
sequencing and alignment to the reference) can diverge from the actual
changes that happened. I am writing "changes" as a shorthand for
misincorporation but also for the point mutations that could occur due to
post-transcriptional RNA editing of some kind. We can't distinguish between
them.

## Installing the code

```sh
$ pip install viral-rna-simulation
```

## Running the simulation

Here's an example invocation and output:

```sh
$ viral-rna-simulation --cells 12 --steps 100 --mutation-rate 0.01 --genome-length 100 --ratio 2
RNA molecules: 1701
Total RNA molecule replications: 1689
Actual mutations:
  Total: 1762
    Of which, 1006 in (+) RNA molecules
    Of which, 756 in (-) RNA molecules
  Rate: 0.010432
  From/to: AC:155, AG:174, AT:148, CA:152, CG:139, CT:120, GA:152, GC:151, GT:137, TA:129, TC:146, TG:159
Apparent mutations:
  Total: 7125
  From/to: AC:719, AG:727, AT:622, CA:703, CG:440, CT:334, GA:639, GC:696, GT:658, TA:564, TC:518, TG:505
```

## Command-line options

Use `--help` to see the following:

```sha
$ viral-rna-simulation --help
usage: viral-rna-simulation [-h] [--cells N] (--genome-length N | --genome ACGT...)
                            [--steps N] [--mutation-rate N] [--ratio N]

Simulate the polymerization of viral RNA within cells and the subsequent RNA sample preparation,
sequencing, and alignment to the viral genome.

options:
  -h, --help         show this help message and exit
  --cells N          The number of cells to simulate. This is purely for speed-up so multiple cells
                     can be run in parallel. All cells are seeded with one copy of the same (+) RNA
                     molecule. (default: 1)
  --genome-length N  Specify the viral RNA genome length. A random genome of the given length will be
                     used. Incompatible with --genome. (default: None)
  --genome ACGT...   Specify the infecting viral genome. Incompatible with --genome. (default: None)
  --steps N          The number of replication steps to simulate. In each replication step, a random
                     RNA from each cell will be chosen to replicate. The replication will create the
                     reverse complement sequence with the (+/-) sense flipped. The new molecule will be
                     added to the RNA population within the cell. (default: 1000)
  --mutation-rate N  The per-nucleotide mutation (polymerase misincorporation) rate, applied during
                     RNA molecule replication. (default: 0.001)
  --ratio N          The number of (+) RNA molecules to make from each (-) molecule. In coronaviruses
                     this ratio is often in the range of 10 to 100. (default: 1)
```

<!--
## References

[1] [C→U transition biases in SARS-CoV-2: still rampant 4 years fromthe start of the COVID-19 pandemic](https://pubmed.ncbi.nlm.nih.gov/39475243/), [Peter Simmonds](https://www.ndm.ox.ac.uk/team/peter-simmonds) mBio. 2024 Dec 11;15(12):e0249324. doi: 10.1128/mbio.02493-24. Epub 2024 Oct 30.
-->
