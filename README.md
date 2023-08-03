Author: Alek Westover
TODOs:
- performance testing
- cache locality
- tune

Description:
parallel k-merge

Say that you have k (a power of 2 for now) sorted lists
and want to merge them into a single sorted list using p (a power of 2)
processors.
We make a tree with p leaves, where each vertex in the tree corresponds to
merging two lists. Parallelism across these vertices is easy.
Within each of these pair merges we can also achieve parallelism via merging odd
/ even halves of the array first and then doing some swaps in parallel.

