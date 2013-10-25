LDSP - Linear Detection of Selected in Pooled Sequences


We break up the process into two phases.

Phase I:
The Intuition is that data at each SNP are binomial counts, which help estimate the frequency of a SNP in a pool, but they don't
tell you the frequency exactly, they are noisy. But by combining information across multiple corrected SNPs, you can improve the estimated frequency of the test SNP

Phase II:
After we estimate the frequency of the putatively selected SNP in each replicate population, we estimate the group effect using a linear model which also allows us to model genetic drift with a normal error term.
The idea is that in the positively selected population, the group effect will be positive while negative in the negatively selected population, and 0 in the neutrally evolving population.