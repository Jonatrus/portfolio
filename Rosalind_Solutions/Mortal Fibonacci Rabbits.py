""""
Problem

Recall the definition of the Fibonacci numbers from "Rabbits and Recurrence Relations", which followed 
the recurrence relation Fn=Fn-1+Fn-2 and assumed that each pair of rabbits reaches maturity in one month 
and produces a single pair of offspring (one male, one female) each subsequent month.
Our aim is to somehow modify this recurrence relation to achieve a dynamic programming solution in the case that 
all rabbits die out after a fixed number of months. See Figure 4 for a depiction of a rabbit tree in which rabbits 
live for three months (meaning that they reproduce only twice before dying).

Given: Positive integers n<=100 and m<=20.

Return: The total number of pairs of rabbits that will remain after the n-th month if all rabbits live for m months.

EG. Input: 6 3, output: 4
"""

# F(n) = F(n-1) + F(n-2) - F(n-m)

# Where:
# F(n-1) + F(n-2) is the standard Fibonacci births
# F(n-m) is the deaths — the population that was born m months ago

n = 93 # How many months to iterate.
m = 16 # How many months the rabbits live.

# I will use a sliding window approach.

from collections import deque
# Each index stores rabbit pairs of a certain age
# left = newborns, right = oldest
rabbits = deque([1] + [0] * (m - 1), maxlen=m) # One at index zero because we start with a population of 1.
# Oldest rabbits die automatically using maxlength parameter

for _ in range(1,n):
      # Rabbits older than 0 months can reproduce
        newborns = sum(list(rabbits)[1:])
        # Add newborns to the front and oldest will be popped off.
        rabbits.appendleft(newborns)

print(sum(rabbits))


#------------ Another Solution ---------
""" This is what I first attempted but couldn't make it work. This is similiar to what i coded for 'Rabbits and Reccurence Relations'.
I revisted this code with input from the solutions forum on rosalind for this problem.
"""
n = 93 # How many months to iterate.
m = 16 # How many months the rabbits live.

population=[1]*100 # 100 because the limits the problem states is n<=100, and we will write over the values in this approach. making all the values 1 ensure the populations atarts at one.

for i in range(2,n):
  # Increase population by applying fibonnaci sequence.
  population[i]=population[i-1]+population[i-2]
  # Calculate deaths
  if i>=m:
    population[i]-=population[(i-m)-1]  
print(population[i])