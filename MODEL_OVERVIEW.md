# Polyploid Crossing Model: A Plain-Language Overview

## What is this project?

This project helps conserve **slickspot peppergrass** (*Lepidium papilliferum*), a rare plant found only in southwestern Idaho. The species is at risk and has small, fragmented populations. We built a computer model to figure out the best way to cross-pollinate plants in order to maintain genetic diversity.

## The biology in brief

### Every plant carries four copies of a key gene

Most animals (including humans) are *diploid* — they carry two copies of each gene. Slickspot peppergrass is *tetraploid*, meaning it carries **four copies** of each gene. This project focuses on one gene in particular: the **S-locus**, which controls self-incompatibility.

Across all the plants we sampled, we found **94 different versions** (alleles) of this gene. Each individual plant holds exactly 4 of them — for example, one plant might carry alleles 3, 17, 42, and 58.

### Plants reject their own pollen

Slickspot peppergrass has a built-in mechanism to prevent self-fertilization called **self-incompatibility (SI)**. When pollen lands on a flower, the plant checks: *does this pollen carry any of the same S-alleles that I have?* If yes, the pollen is rejected. This forces plants to mate with genetically different partners, which promotes diversity.

### Rare alleles have an advantage

Because of SI, plants with **rare alleles** can mate with more partners — their pollen is less likely to be rejected. This creates a natural balancing force called **negative frequency-dependent selection (NFDS)**. In an ideal world, NFDS would drive all alleles to equal frequency over time.

In practice, small populations lose rare alleles through random chance (genetic drift) faster than NFDS can protect them. Once an allele is gone, it is gone forever.

## What the model does

### Step 1: Figure out who can mate with whom

For every possible pair of plants, the model checks SI compatibility. When two plants cross, each parent contributes 2 of their 4 alleles to the offspring (there are 6 possible combinations of 2-from-4). The model filters out any pollen combinations that would be rejected, then calculates all possible offspring genotypes and their probabilities.

### Step 2: Find the best crosses

Instead of letting plants mate randomly, the model uses **mathematical optimization** to find the ideal set of crosses. The goal is to produce offspring whose allele frequencies are as close to equal as possible — the NFDS equilibrium.

Think of it like a recipe: the optimizer determines how much weight to give each possible cross so that the next generation has the most balanced set of alleles.

### Step 3: Protect rare alleles

Pure optimization can sacrifice rare alleles if doing so reduces overall variance. To prevent this, the model adds three safety mechanisms:

1. **Elite preservation**: The top 10% of individuals carrying the rarest alleles are kept in the next generation unchanged.
2. **Mandatory rare-allele crosses**: For every allele carried by only 1-2 plants, the model forces at least one cross that can produce offspring with that allele.
3. **Optimizer penalty**: The optimizer is penalized for any solution that would reduce rare allele frequencies below their target.

### Step 4: Model population growth

Real populations don't stay the same size forever. The model uses **logistic growth** to simulate how populations change over time:

- Populations grow when they are below their habitat's carrying capacity.
- Growth slows as the population approaches the maximum the habitat can support.
- Random variation (demographic stochasticity) means small populations can fluctuate dramatically from year to year.

The model also tracks **effective population size (Ne)** — a measure of how many plants are actually contributing genetically. Because SI restricts who can mate with whom, Ne is always smaller than the number of plants you can count in the field.

## The four strategies compared

The model runs simulations under four different management strategies:

| Strategy | What it does | Key result |
|----------|-------------|------------|
| **Random mating** | No management — plants mate at random | Alleles are lost quickly; common alleles become even more common |
| **Optimized crossing** | Mathematically optimal crosses, but no allele protection | Good variance reduction, but rare alleles can still be lost |
| **Optimized + Preservation** | Adds elite retention, mandatory rare crosses, and optimizer penalty | Prevents allele loss with only a small cost to convergence speed |
| **Optimized + Preservation + Demography** | Also models population growth toward carrying capacity | Best overall: larger populations provide more genetic "room" and buffer against drift |

## Key results

The model was applied to **124 plants across 24 populations** (4 major populations with 15-31 individuals each).

**Allele preservation works.** Under random mating, populations lose 5-25 alleles within 5 generations. With the preservation strategy, zero alleles are lost.

**Population growth is synergistic.** When preservation is combined with logistic growth, variance reduction improves from roughly 75% to 90% compared to random mating. Growing the population gives more individuals to carry alleles, more compatible mating pairs, and a stronger buffer against random loss.

**All populations are at risk.** Effective population sizes range from about 15 to 28 — well below the conservation threshold of 50 needed to avoid short-term genetic problems. This underscores the importance of active management.

## Practical takeaway

The model provides a ranked list of recommended crosses for each population. A conservation manager can use these recommendations to decide which plants to cross-pollinate each season, prioritizing crosses that:

1. Produce offspring with the most balanced allele frequencies
2. Ensure no rare allele is left without a carrier
3. Support population growth toward a sustainable size

The combination of smart crossing and population growth gives these small, fragmented populations the best chance of maintaining the genetic diversity they need to survive.
