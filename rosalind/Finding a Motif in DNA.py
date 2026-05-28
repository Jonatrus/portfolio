"""
Problem: Finding a Motif in DNA

A substring t is contained in string s if t appears as a contiguous
collection of symbols inside s (and t is no longer than s).

Indexing rules:
- String positions are 1-based.
- The position of a symbol is the number of symbols to its left,
  including itself.
- Example:
    s = "AUGCUUCAGAAAGGUCUUACG"
    The positions of 'U' are:
    2, 5, 6, 15, 17, 18

Substring notation:
- A substring of s can be written as s[j:k]
- j = starting position
- k = ending position
- Example:
    s = "AUGCUUCAGAAAGGUCUUACG"
    s[2:5] = "UGCU"

The location of a substring s[j:k] is its starting position j.

If t appears multiple times in s, return all starting positions.

Given:
- Two DNA strings s and t
- Each of length at most 1 kbp

Return:
- All starting positions (1-based) where t appears in s
"""

with open("rosalind_subs.txt") as f:
    s, t = f.read().splitlines()

# Find all occurrences using list comprehension
positions = [i+1 for i in range(len(s) - len(t) + 1) if s.startswith(t, i)]


with open("output.txt", "w") as f:
    f.write(" ".join(map(str, positions)))