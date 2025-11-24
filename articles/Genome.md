# The Genome

**AlphaSimR** is built upon the basics of genetics for the stochastic
simulation of crop and animal breeding programmes. This vignette will
provide an overview of the genome, it’s inheritance, and the tools used
to read it with supporting examples and code.

## What is a genome?

The genome is the complete set of genetic information within an
organism. It is made up of deoxyribonucleic acid (DNA), consisting of
sequences of four types of nucleotide bases: adenine (A), cytosine (C),
guanine (G), and thymine (T). Within one DNA molecule, there are two
strands joined by base-pairs: A binds with T and C binds with G. This
forms a double helix structure which can be visualised as a spiral
ladder. These DNA molecules are condensed around supporting protein to
form a thread-like structure known as chromosomes. Due to the different
demographic history across species, how the DNA is package can vary
substantial including the number of chromosome and the number of
chromosome copies (ploidy). For example, humans have 23 pairs of
chromosomes (22 autosomnal, one non-autosomnal) for ~three billion
base-pairs. Whereas cattle have 30 pairs of chromosomes for fewer
base-pairs at ~2.7 billion. As these species are both diploid, in each
chromosome pair, one has been inherited from the father (paternal) and
one has been inherited from the mother (maternal).

### Genes, alleles, loci, genotypes, and haplotypes.

Although similar, there are essential differences in the terminology
used. This next section will provide a short glossary for genes, allele,
loci, genotypes, and haplotypes.

Genes are segments of DNA that are translated to proteins. Typically
they will contain a few hundred to more than two billion bases-pairs
with a codon (triplet base-pair) coding for either an amino acid or to
start/stop.

An allele in a gene is a gene variant. More generally, an allele
represents a possible sequence of bases at a particular site on a genome
(known as a locus). When biallelic (two variants per gene) the allele is
either ancestral (denoted 0) or a mutant (denoted 1). Per the name, new
alleles arise from mutation.

A genotype is the combination of alleles across chromosome copies that
an organism has at a specified locus or across loci. For instance, there
are three possible genotypes for a diploid at a biallelic locus with
ancestral allele A and mutant allele a: ancestral homozygous (aa) with
allele dosage 0, heterozygous (aA or Aa) with allele dosage 1, and
mutant homozygous (AA) with allele dosage 2. Genotypes can differ
between individuals due to random segregation in inheritance.

A haplotype (haploid genotype) describes the physical grouping of
alleles along the same chromosome. This means that these group of
alleles tend to be inherited together. New haplotypes are created
through recombination.

What does this look like in **AlphaSimR**?

Let’s start with creating the founding genome using `runMacs`. As
mentioned, cattle have 30 pairs of chromosomes, but for the purpose of
this simulation we are only interested in three. When we inspect the
founderGenome we see that it has a ploidy of two with 12 loci across
three chromosomes (with four segregating sites on each).

``` r
founderGenome = runMacs(nInd = 10, 
                        nChr = 3, 
                        segSites = 4, 
                        species = "CATTLE")
# Set simulation parameters
SP = SimParam$new(founderGenome)
SP$setTrackRec(TRUE)
# Inspect the founderGenome
founderGenome
#> An object of class "MapPop" 
#> Ploidy: 2 
#> Individuals: 10 
#> Chromosomes: 3 
#> Loci: 12
```

Then we can observe the genotypes with `pullSegSiteGeno`. Each row
represents one of the ten individuals and each column represents a
defined loci labelled as the Chromosome_Locus. All will have a value of
zero to two equivalent to the allele dosage. Commonly, this would be
homozygous ancestral (0), heterozygous (1), and homozygous mutant (2).
Since this is a diploid, we are observing these genotypes across the
pair of chromosomes. Where the ancestral allele is not known, the allele
dosage will be randomly assigned.

``` r
pullSegSiteGeno(founderGenome)
#>    1_1 1_2 1_3 1_4 2_1 2_2 2_3 2_4 3_1 3_2 3_3 3_4
#> 1    2   2   2   2   1   1   2   1   0   1   1   2
#> 2    1   1   1   2   0   0   2   0   0   0   2   2
#> 3    0   1   2   2   0   2   2   0   0   0   2   2
#> 4    1   2   2   2   0   1   2   1   0   1   2   1
#> 5    2   2   1   2   1   1   1   0   0   0   2   1
#> 6    1   0   1   2   1   2   1   0   0   1   2   2
#> 7    0   2   1   2   0   0   2   1   1   1   2   2
#> 8    2   2   2   2   1   1   2   2   0   0   2   0
#> 9    1   2   1   1   0   1   2   1   1   1   2   0
#> 10   0   2   2   1   1   2   2   1   0   1   2   0
```

The haplotypes can be observed with `pullSegSiteHaplo`. Each row
represents the sequence of alleles on each chromosome per individual.
Since this is diploid, we observe two for each individual: one maternal
(id_1) and the other paternal (id_2). All will have a value of 0 or 1
for either ancestral or mutant allele respectingly. As before, each
column represents a defined loci.

``` r
pullSegSiteHaplo(founderGenome)
#>      1_1 1_2 1_3 1_4 2_1 2_2 2_3 2_4 3_1 3_2 3_3 3_4
#> 1_1    1   1   1   1   1   1   1   1   0   1   0   1
#> 1_2    1   1   1   1   0   0   1   0   0   0   1   1
#> 2_1    0   0   1   1   0   0   1   0   0   0   1   1
#> 2_2    1   1   0   1   0   0   1   0   0   0   1   1
#> 3_1    0   1   1   1   0   1   1   0   0   0   1   1
#> 3_2    0   0   1   1   0   1   1   0   0   0   1   1
#> 4_1    1   1   1   1   0   1   1   0   0   0   1   0
#> 4_2    0   1   1   1   0   0   1   1   0   1   1   1
#> 5_1    1   1   1   1   0   0   1   0   0   0   1   0
#> 5_2    1   1   0   1   1   1   0   0   0   0   1   1
#> 6_1    1   0   0   1   1   1   0   0   0   1   1   1
#> 6_2    0   0   1   1   0   1   1   0   0   0   1   1
#> 7_1    0   1   1   1   0   0   1   1   0   0   1   1
#> 7_2    0   1   0   1   0   0   1   0   1   1   1   1
#> 8_1    1   1   1   1   1   0   1   1   0   0   1   0
#> 8_2    1   1   1   1   0   1   1   1   0   0   1   0
#> 9_1    1   1   0   0   0   1   1   0   1   0   1   0
#> 9_2    0   1   1   1   0   0   1   1   0   1   1   0
#> 10_1   0   1   1   1   0   1   1   0   0   1   1   0
#> 10_2   0   1   1   0   1   1   1   1   0   0   1   0
```

## How are genomes inherited?

Genetic code is inherited from an individuals mother and father. This
involves cells known as gametes: eggs (from mother) and sperm (from
father). These cells have half the number of chromosomes to other cells
in the species. For example, in diploid species like humans gametes are
haploid cells. These have been developed through interphase (DNA
replication) and meiosis (cell division). The mechanism of interphase
and meiosis are outside the scope of this vignette, however is important
due to their contribution to genetic variation.

### Genetic Variation: the DNA lottery

Broadly speaking, there are three contributions to genetic variation:
mutation, crossing over (or recombination), and independent assortment.

Mutation occurs during DNA replication in interphase. Although the rate
is very small (~2.5x10⁻⁸ mutation per nucleotide), when considering the
whole genome mutation is not negligible.

Both crossing over and independent assortment occurs in the first part
of meiosis known as the reductional division. Crossing over refers to
the swapping of genetic material at random between chromatids resulting
in variation. Independent Assortment refers to the random orientation of
homologous chromosome pairs resulting in gametes with many different
assortments of chromosomes. (recall there are 23 pairs in a human).

Within genetics, albeit quantitative or population, we are interested in
how individual’s genomes differ to the next. Thus, it is this variation
we wish to capture to understand how genetics influence desired (or
undesired) traits in crop and animal breeding.

In **AlphaSimR**:

1.  Random mutations can be added to individuals in populations with
    `mutate`. For the purpose of demonstration, the mutation is at a
    much higher rate than would be expected.

``` r
basePop = newPop(founderGenome)

mutatedBasePop = mutate(basePop, mutRate = 0.1, simParam = SP)

pullSegSiteGeno(basePop)
#>    1_1 1_2 1_3 1_4 2_1 2_2 2_3 2_4 3_1 3_2 3_3 3_4
#> 1    2   2   2   2   1   1   2   1   0   1   1   2
#> 2    1   1   1   2   0   0   2   0   0   0   2   2
#> 3    0   1   2   2   0   2   2   0   0   0   2   2
#> 4    1   2   2   2   0   1   2   1   0   1   2   1
#> 5    2   2   1   2   1   1   1   0   0   0   2   1
#> 6    1   0   1   2   1   2   1   0   0   1   2   2
#> 7    0   2   1   2   0   0   2   1   1   1   2   2
#> 8    2   2   2   2   1   1   2   2   0   0   2   0
#> 9    1   2   1   1   0   1   2   1   1   1   2   0
#> 10   0   2   2   1   1   2   2   1   0   1   2   0

pullSegSiteGeno(mutatedBasePop)
#>    1_1 1_2 1_3 1_4 2_1 2_2 2_3 2_4 3_1 3_2 3_3 3_4
#> 1    2   2   1   1   1   2   2   1   1   2   1   2
#> 2    1   2   1   2   0   0   2   0   0   1   2   1
#> 3    0   1   2   1   0   2   2   0   1   0   2   2
#> 4    2   2   1   1   0   0   2   2   0   1   2   2
#> 5    2   2   1   2   0   1   1   0   0   0   2   1
#> 6    1   1   1   1   2   2   1   0   0   0   2   2
#> 7    0   2   1   2   0   0   2   1   1   1   2   2
#> 8    2   2   2   2   1   1   1   2   0   0   2   1
#> 9    1   2   1   0   0   1   1   1   1   1   2   0
#> 10   2   2   2   1   1   1   2   1   1   1   2   1
```

2.  Independent assortment and recombination can be observed across
    generations.

The first step is to make the cross in the base population to produce
the second generation. For this, `randCross` can be used to create one
cross, producing a pair of full siblings.

``` r
secondGenPop = randCross(pop = basePop, nCrosses = 1, nProgeny = 2, simParam = SP)

str(secondGenPop)
#> Formal class 'Pop' [package "AlphaSimR"] with 18 slots
#>   ..@ id     : chr [1:2] "11" "12"
#>   ..@ iid    : int [1:2] 11 12
#>   ..@ mother : chr [1:2] "4" "4"
#>   ..@ father : chr [1:2] "6" "6"
#>   ..@ sex    : chr [1:2] "H" "H"
#>   ..@ nTraits: int 0
#>   ..@ gv     : num[1:2, 0 ] 
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : NULL
#>   .. .. ..$ : NULL
#>   ..@ pheno  : num[1:2, 0 ] 
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : NULL
#>   .. .. ..$ : NULL
#>   ..@ ebv    : num[1:2, 0 ] 
#>   ..@ gxe    : list()
#>   ..@ fixEff : int [1:2] 1 1
#>   ..@ misc   : list()
#>   ..@ miscPop: list()
#>   ..@ nInd   : int 2
#>   ..@ nChr   : int 3
#>   ..@ ploidy : int 2
#>   ..@ nLoci  : int [1:3] 4 4 4
#>   ..@ geno   :List of 3
#>   .. ..$ : raw [1, 1:2, 1:2] 0e 09 0e 0d
#>   .. ..$ : raw [1, 1:2, 1:2] 06 06 0c 03
#>   .. ..$ : raw [1, 1:2, 1:2] 06 0c 0c 0e

# Collect the iid for mother and father
mother = as.integer(secondGenPop@mother[1])
father = as.integer(secondGenPop@father[1])
```

Then, the haplotypes can be extracted for the mother, father, and two
progeny to observe the inheritance. Recall that id_1 is the maternal
haplotype (inherited from mother) and id_2 is the paternal haplotype
(inherited from father). For clarity, no mutation is included in this
example. All variation observed is due to independent assortment and/or
recombination.

``` r
pullSegSiteHaplo(basePop[basePop@iid == mother,])
#>     1_1 1_2 1_3 1_4 2_1 2_2 2_3 2_4 3_1 3_2 3_3 3_4
#> 4_1   1   1   1   1   0   1   1   0   0   0   1   0
#> 4_2   0   1   1   1   0   0   1   1   0   1   1   1
```

``` r
pullSegSiteHaplo(basePop[basePop@iid == father,])
#>     1_1 1_2 1_3 1_4 2_1 2_2 2_3 2_4 3_1 3_2 3_3 3_4
#> 6_1   1   0   0   1   1   1   0   0   0   1   1   1
#> 6_2   0   0   1   1   0   1   1   0   0   0   1   1
```

``` r
pullSegSiteHaplo(secondGenPop)
#>      1_1 1_2 1_3 1_4 2_1 2_2 2_3 2_4 3_1 3_2 3_3 3_4
#> 11_1   0   1   1   1   0   1   1   0   0   1   1   0
#> 11_2   1   0   0   1   0   1   1   0   0   0   1   1
#> 12_1   0   1   1   1   0   0   1   1   0   0   1   1
#> 12_2   1   0   1   1   1   1   0   0   0   1   1   1
```

In the above, independent assortment will always be observed. This is
clear since despite being full siblings they inherit a different
arrangement of haplotypes from the mother and the father. We may observe
changes to the haplotypes via recombination between the four sites
across a chromosome (either chromosome one, two, and/or three). However,
this depends on the rate of recombination. Nevertheless, even if present
it can be hard or even impossible to observe through the haplotypes
alone, even in this small sample size. Instead we can use the
`pullIbdHaplo` to observe this. For further details on recombination,
please see the vignette on Recombination.

``` r
pullIbdHaplo(basePop[basePop@iid == mother,])
#>     1_1 1_2 1_3 1_4 2_1 2_2 2_3 2_4 3_1 3_2 3_3 3_4
#> 4_1   7   7   7   7   7   7   7   7   7   7   7   7
#> 4_2   8   8   8   8   8   8   8   8   8   8   8   8
pullIbdHaplo(basePop[basePop@iid == father,])
#>     1_1 1_2 1_3 1_4 2_1 2_2 2_3 2_4 3_1 3_2 3_3 3_4
#> 6_1  11  11  11  11  11  11  11  11  11  11  11  11
#> 6_2  12  12  12  12  12  12  12  12  12  12  12  12
pullIbdHaplo(secondGenPop)
#>      1_1 1_2 1_3 1_4 2_1 2_2 2_3 2_4 3_1 3_2 3_3 3_4
#> 11_1   8   8   8   8   7   7   7   7   7   8   7   7
#> 11_2  11  11  11  11  12  12  12  12  12  12  12  12
#> 12_1   8   8   8   8   8   8   8   8   7   7   8   8
#> 12_2  11  12  12  12  11  11  11  11  12  11  11  11
```

## How are genomes read?

For a genome to be read, first it needs to be sequenced. The techniques
for this include Sanger sequencing or next-generation sequencing. These
are then used to identify single nucleotide polymorphisms (SNPs). SNPs
highlight variation at specific loci; equivalent to markers. Other
terminology include traits and QTL: Quantitative Trait Locus. Traits are
characteristics of interest such as milk yield or growth rate that have
a genetic component. QTL are a region of DNA that’s associated with a
specific trait that varies in the population. This could be identifiable
by a marker (eg. SNP) or could only be correlated.

In **AlphaSimR**:

Let’s start a new population.

``` r
founderGenome = runMacs(nInd = 10, 
                        nChr = 3, 
                        segSites = 4, 
                        species = "CATTLE")
# Set simulation parameters
SP = SimParam$new(founderGenome)
```

Then, add an additive trait with three QTLs per chromosome using
`SP$addTraitA` and one SNP per chromosome with `SP$addSnpChip`.

``` r
SP$addTraitA(nQtlPerChr = 3)
SP$addSnpChip(nSnpPerChr = 1)
```

After, create a base population.

``` r
basePop = newPop(founderGenome, simParam = SP)

# Inspect basePop
basePop
#> An object of class "Pop" 
#> Ploidy: 2 
#> Individuals: 10 
#> Chromosomes: 3 
#> Loci: 12 
#> Traits: 1
```

The SNPs can be read with `pullSnpGeno` or, if the location in the
genome is known, `pullMarkerGeno` can be used.

``` r
# From SNP
SNP = pullSnpGeno(pop = basePop, snpChip =1)
print(SNP)
#>    1_3 2_4 3_2
#> 1    0   1   1
#> 2    1   2   1
#> 3    0   1   1
#> 4    0   1   2
#> 5    1   1   1
#> 6    1   1   1
#> 7    2   2   0
#> 8    1   2   1
#> 9    1   2   0
#> 10   1   2   2

# From marker
location = colnames(SNP)
pullMarkerGeno(basePop, markers = location)
#>    1_3 2_4 3_2
#> 1    0   1   1
#> 2    1   2   1
#> 3    0   1   1
#> 4    0   1   2
#> 5    1   1   1
#> 6    1   1   1
#> 7    2   2   0
#> 8    1   2   1
#> 9    1   2   0
#> 10   1   2   2
```

The defined trait (Trait 1) can be inspected via the genetic values
calculated for each individual.

``` r
basePop@gv
#>            Trait1
#>  [1,] -0.80998556
#>  [2,] -0.08194541
#>  [3,] -0.45570651
#>  [4,]  0.48777140
#>  [5,]  0.32859247
#>  [6,]  2.66849231
#>  [7,] -0.72505222
#>  [8,] -1.02403006
#>  [9,] -0.27509763
#> [10,] -0.11303880
```

The genetic map for the QTLs can be retrieved via `getQtlMap`. Here the
site of each QTL is identified per chromosome. The pos notes the
position of the QTL on each chromosome in Morgans.

``` r
getQtlMap(trait = 1, simParam = SP)
#>      id chr site        pos
#> 1_1 1_1   1    1 0.00000000
#> 1_2 1_2   1    2 0.05377956
#> 1_4 1_4   1    4 0.33321374
#> 2_1 2_1   2    1 0.00000000
#> 2_2 2_2   2    2 0.28100082
#> 2_3 2_3   2    3 0.34157652
#> 3_1 3_1   3    1 0.00000000
#> 3_3 3_3   3    3 0.37847421
#> 3_4 3_4   3    4 0.58925124
```

### What is a genetic map?

A genetic map (also known as linkage map) provides the genetic distance
in Morgans between genetic markers, QTLs, and/or segregating sites on a
chromosome. This differs from a physical map which gives the physical
distance in the DNA sequence by the number of nucleotides. Thus, is
based on the actual structure of the genetic material. Although
**AlphaSimR** can import a physical map for all loci through the
`lociMap-class`, it can only output genetic maps of defined segregating
sites, including QTLs, SNPs, and all defined loci. For this, the
observed recombination frequencies between loci are used to show how
genetic information is shuffled in chromosomes. Thus, highlighting how
markers are inherited.

For a genetic map of all defined loci (or segregating sites),
`getGenMap` is used. `getSnpMap` can be used for the genetic map of the
SNP(s) on each chromosome. As hinted previously, `getQtlMap` can be used
for the genetic map of the QTL(s).

``` r
getGenMap(SP)
#>     id chr        pos
#> 1  1_1   1 0.00000000
#> 2  1_2   1 0.05377956
#> 3  1_3   1 0.07121389
#> 4  1_4   1 0.33321374
#> 5  2_1   2 0.00000000
#> 6  2_2   2 0.28100082
#> 7  2_3   2 0.34157652
#> 8  2_4   2 0.60790699
#> 9  3_1   3 0.00000000
#> 10 3_2   3 0.26565846
#> 11 3_3   3 0.37847421
#> 12 3_4   3 0.58925124
getSnpMap(snpChip = 1, simParam = SP)
#>      id chr site        pos
#> 1_3 1_3   1    3 0.07121389
#> 2_4 2_4   2    4 0.60790699
#> 3_2 3_2   3    2 0.26565846
```

## How is this used in plant and animal breeding?

By understanding the genome, breeders can produce crop and livestock
with desirable traits for use in specific production systems. For
example, a breeder may desire to have polled cattle to increase animal
welfare and human safety. When considering the breeding design for a
species, breed, and or system, there are many factors to consider with
potential to introduce unintended consequences. These include balancing
the trade off between short-term genetic gain and long-term loss in
genetic variation as well as antagonistic relationships between desired
traits commonly observed between production and functional traits in
animals. **AlphaSimR** provides a fast and inexpensive way to explore
and test a wide range of breeding designs for long term genetic
improvement across plants and animals through forward-in-time stochastic
simulations.

For a more detailed overview of AlphaSimR in plant and animal breeding
with exercises, please see the edx course: Breeding Programme Modelling
with AlphaSimR. Available at:
<https://learning.edx.org/course/course-v1:EdinburghX+ISB001+1T2021/home>
