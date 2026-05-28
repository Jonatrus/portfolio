"""
Problem 2 — Noise Filtering

Title: NMRP2

Given:
A list of intensities S and a threshold T.

Task:
Return only peaks (from Problem 1) whose intensity ≥ T.

Return:
Filtered peak indices.

Example:

Input:
S = [0, 2, 5, 3, 1, 4, 6, 2]
T = 5

Output:
[2, 6]
"""

S = [0, 2, 5, 3, 1, 4, 6, 2]
T = 5
S_filtered = []

for num in S:
    if num >= T:
        S_filtered.append(num)

print(S_filtered)