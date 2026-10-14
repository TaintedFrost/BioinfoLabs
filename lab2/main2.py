import math

seq="ACGCGTGCCA"
A=seq.count("A")
C=seq.count("C")
T=seq.count("T")
G=seq.count("G")
Na=0.05
length=len(seq)
GC=(G+C)/length

Tm1=4*(G+C)+2*(A+T)
Tm2=81.5+16.6*math.log10(Na)+41*GC-(600/length)

print("A=",A)
print("T=",T)
print("C=",C)
print("G=",G)
print("Tm1:",Tm1)
print("Tm2:",Tm2)