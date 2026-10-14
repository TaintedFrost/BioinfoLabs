seq="ATCGCGTA"
A=seq.count("A")
C=seq.count("C")
T=seq.count("T")
G=seq.count("G")
Tm=4*(G+C)+2*(A+T)

print("A=",A)
print("T=",T)
print("C=",C)
print("G=",G)
print("Tm:", Tm)